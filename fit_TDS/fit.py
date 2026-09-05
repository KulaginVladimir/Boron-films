import os
import argparse
import csv
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import softmax

from tds_model import T0, n_B, simulate_TDS

for name in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(name, "1")

area = 1.0e-4
Tf_fit = 1100.0

repo = Path(__file__).resolve().parent.parent
experiments_file = repo / "experimental_data" / "experiments_5.csv"

n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))

experiments = None
n_residuals = None
N_traps = None
surface_model = None

starts = {
    3: {
        "fractions": np.array([0.60, 0.10, 0.30]),
        "energies": np.array([0.90, 1.55, 2.30]),
        "Kr0": 1.0e-25,
        "Er": 0.80,
    },
    4: {
        "fractions": np.array([0.62, 0.08, 0.05, 0.25]),
        "energies": np.array([0.85, 1.05, 1.85, 2.35]),
        "Kr0": 1.0e-25,
        "Er": 0.80,
    },
    5: {
        "fractions": np.array(
            [0.19070939, 0.27484574, 0.27484573, 0.16451335, 0.09508579]
        ),
        "energies": np.array([0.92443877, 1.07, 1.19, 1.32677295, 1.62738778]),
        "Kr0": 2.69456812e-25,
        "Er": 0.84338380,
    },
    6: {
        "fractions": np.array(
            [0.12654303, 0.16581731, 0.22846126, 0.22846125, 0.15936622, 0.09135093]
        ),
        "energies": np.array(
            [0.91480892, 1.05522346, 1.15, 1.25, 1.41705683, 1.73039792]
        ),
        "Kr0": 1.39328729e-23,
        "Er": 0.99708420,
    },
}

sink_start = {
    "fractions": np.array([0.30, 0.20, 0.20, 0.15, 0.15]),
    "energies": np.array([1.00, 1.30, 1.60, 2.00, 2.35]),
}


def read_tds(path):
    """Read one experimental TDS file and remove rows with invalid data."""
    data = np.genfromtxt(
        path,
        delimiter=",",
        skip_header=10,
        missing_values="--",
        filling_values=np.nan,
        invalid_raise=False,
    )
    return data[np.all(np.isfinite(data[:, :4]), axis=1)]


def trapezoid_weights(T):
    """Return point weights corresponding to trapezoidal integration over T."""
    q = np.empty_like(T)
    q[0] = 0.5 * (T[1] - T[0])
    q[-1] = 0.5 * (T[-1] - T[-2])
    q[1:-1] = 0.5 * (T[2:] - T[:-2])
    return q


def load_experiments():
    """
    Load the experimental data and prepare fixed inventories and
    normalized least-squares weights.
    """
    out = []

    with open(experiments_file, newline="", encoding="utf-8-sig") as file:
        rows = list(csv.DictReader(file))

    for row in rows:
        experiment_id = row["id"]
        raw = read_tds(experiments_file.parent / f"{experiment_id}.csv")
        fit = raw[(raw[:, 1] >= T0) & (raw[:, 1] <= Tf_fit)]

        thickness = float(row["film_thickness_nm"]) * 1e-9
        beta = float(row["heating_rate_K_per_s"])

        # The initial D inventory is fixed by the experimental TDS integral.
        order_time = np.argsort(fit[:, 0])
        time = fit[order_time, 0]
        flux_time = (fit[order_time, 2] + fit[order_time, 3]) / area
        inventory = np.trapezoid(flux_time, time) / (thickness * n_B)

        order_temperature = np.argsort(fit[:, 1])
        temperature = fit[order_temperature, 1]
        flux = (fit[order_temperature, 2] + fit[order_temperature, 3]) / area

        temperature_span = temperature[-1] - temperature[0]
        l1_norm = np.trapezoid(np.abs(flux), temperature)
        mean_abs_flux = l1_norm / temperature_span

        # Scale the pointwise residual so that its squared norm reproduces
        # the temperature-averaged normalized L2 error for this experiment.
        scale = (
            np.sqrt(trapezoid_weights(temperature) / temperature_span) / mean_abs_flux
        )

        out.append(
            {
                "id": experiment_id,
                "beta": beta,
                "thickness": thickness,
                "inventory": inventory,
                "T": temperature,
                "J": flux,
                "l1_norm": l1_norm,
                "scale": scale,
            }
        )

    factor = np.sqrt(len(out))
    for experiment in out:
        experiment["scale"] /= factor

    return out


def init_worker(
    experiments_local,
    n_residuals_local,
    N_traps_local,
    surface_model_local,
):
    """Initialize shared read-only data in each multiprocessing worker."""
    global experiments, n_residuals, N_traps, surface_model
    experiments = experiments_local
    n_residuals = n_residuals_local
    N_traps = N_traps_local
    surface_model = surface_model_local


def make_x0(start):
    fractions = start["fractions"]
    logits = np.log(fractions[:-1] / fractions[-1])

    if surface_model == "sink":
        return np.r_[logits, start["energies"]]

    return np.r_[
        logits,
        start["energies"],
        np.log10(start["Kr0"]),
        start["Er"],
    ]


def bounds():
    """Return lower and upper bounds in optimization coordinates."""
    n_logits = N_traps - 1

    lower = np.r_[
        np.full(n_logits, -12.0),
        np.full(N_traps, 0.60),
    ]
    upper = np.r_[
        np.full(n_logits, 12.0),
        np.full(N_traps, 2.50),
    ]

    if surface_model == "arrhenius":
        lower = np.r_[lower, -30.0, 0.0]
        upper = np.r_[upper, -16.0, 1.50]

    return lower, upper


def diff_step(x):
    """
    Convert desired absolute finite-difference steps to the relative steps
    expected by scipy.optimize.least_squares.
    """
    absolute = np.r_[
        np.full(N_traps - 1, 2.0e-3),
        np.full(N_traps, 5.0e-4),
    ]

    if surface_model == "arrhenius":
        absolute = np.r_[absolute, 2.0e-3, 5.0e-4]

    return absolute / np.maximum(np.abs(x), 1.0e-2)


def parameters(x):
    """
    Convert optimization coordinates to physical model parameters.

    Trap fractions are represented by a softmax with the last logit fixed
    to zero, which enforces positive fractions summing to unity.
    """
    n_logits = N_traps - 1

    fractions = softmax(np.r_[x[:n_logits], 0.0])
    energies = np.asarray(x[n_logits : n_logits + N_traps])

    if surface_model == "sink":
        return fractions, energies, 0.0, 0.0

    Kr0 = 10.0 ** float(x[n_logits + N_traps])
    Er = float(x[n_logits + N_traps + 1])

    return fractions, energies, Kr0, Er


def residual(x):
    """
    Calculate the concatenated normalized residual vector for all TDS spectra.

    Its squared Euclidean norm is the objective minimized by
    scipy.optimize.least_squares.
    """
    fractions, energies, Kr0, Er = parameters(x)
    parts = []

    try:
        for experiment in experiments:
            temperature, flux = simulate_TDS(
                beta=experiment["beta"],
                thickness=experiment["thickness"],
                inventory=experiment["inventory"],
                fractions=fractions,
                energies=energies,
                Kr0=Kr0,
                E_r=Er,
                Tf=Tf_fit,
                surface_model=surface_model,
            )

            delta = np.interp(experiment["T"], temperature, flux) - experiment["J"]
            parts.append(delta * experiment["scale"])

        return np.concatenate(parts)

    except Exception:
        return np.full(n_residuals, 1.0e4)


def evaluate(x):
    """
    Evaluate the global least-squares objective and the normalized L1 error
    of each experiment for a given parameter vector.
    """
    fractions, energies, Kr0, Er = parameters(x)
    errors = []

    for experiment in experiments:
        temperature, flux = simulate_TDS(
            beta=experiment["beta"],
            thickness=experiment["thickness"],
            inventory=experiment["inventory"],
            fractions=fractions,
            energies=energies,
            Kr0=Kr0,
            E_r=Er,
            Tf=Tf_fit,
            surface_model=surface_model,
        )

        delta = np.interp(experiment["T"], temperature, flux) - experiment["J"]

        errors.append(
            np.trapezoid(np.abs(delta), experiment["T"]) / experiment["l1_norm"]
        )

    r = residual(x)
    return float(r @ r), np.asarray(errors)


def main():
    """Run the simultaneous TDS fit."""
    global experiments, n_residuals, N_traps, surface_model

    parser = argparse.ArgumentParser()
    parser.add_argument("N_traps", type=int)
    parser.add_argument("surface_model", nargs="?", default="arrhenius")
    args = parser.parse_args()

    N_traps = args.N_traps
    surface_model = args.surface_model

    experiments = load_experiments()
    n_residuals = sum(len(experiment["T"]) for experiment in experiments)

    start = sink_start if surface_model == "sink" else starts[N_traps]

    x0 = make_x0(start)
    lower, upper = bounds()

    print(f"Surface model: {surface_model}")
    print(f"N traps: {N_traps}")
    print(f"Free parameters: {len(x0)}")
    print(f"Workers: {n_workers}")
    print()

    context = mp.get_context("spawn")

    with ProcessPoolExecutor(
        max_workers=n_workers,
        mp_context=context,
        initializer=init_worker,
        initargs=(experiments, n_residuals, N_traps, surface_model),
    ) as executor:
        fit = least_squares(
            residual,
            x0,
            bounds=(lower, upper),
            method="trf",
            jac="3-point",
            diff_step=diff_step(x0),
            x_scale="jac",
            loss="linear",
            ftol=2.0e-4,
            xtol=2.0e-4,
            gtol=2.0e-4,
            max_nfev=200,
            verbose=2,
            workers=executor.map,
        )

    objective, errors = evaluate(fit.x)
    fractions, energies, Kr0, Er = parameters(fit.x)

    order = np.argsort(energies)
    fractions = fractions[order]
    energies = energies[order]

    print()
    print("fractions =", np.array2string(fractions, precision=8))
    print("energies  =", np.array2string(energies, precision=8))
    print("dE        =", np.array2string(np.diff(energies), precision=8))

    if surface_model == "sink":
        print("surface   = c_m(0,t) = 0")
    else:
        print(f"Kr0       = {Kr0:.8e} m4/s")
        print(f"Er        = {Er:.8f} eV")

    print(f"objective = {objective:.8e}")
    print(f"SciPy cost= {0.5 * objective:.8e}")
    print(f"mean L1   = {np.mean(errors):.8f}")
    print("L1        =", np.array2string(errors, precision=6))
    print(f"nfev={fit.nfev}, njev={fit.njev}")
    print(fit.message)


if __name__ == "__main__":
    mp.freeze_support()
    main()
