import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

import fit as fit_base
import tds_model


nominal_start = {
    "fractions": np.array([0.19070939, 0.54969147, 0.16451335, 0.09508579]),
    "energies": np.array([0.92443877, 1.13050166, 1.32677295, 1.62738778]),
    "Kr0": 2.69456812e-25,
    "Er": 0.84338380,
}


def init_worker(experiments, n_residuals, k0_fixed, p0_fixed):
    fit_base.init_worker(experiments, n_residuals, 4)
    tds_model.k0 = k0_fixed
    tds_model.p0 = p0_fixed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("k0", type=float)
    parser.add_argument("p0", type=float)
    args = parser.parse_args()

    fit_base.N_traps = 4
    fit_base.experiments = fit_base.load_experiments()
    fit_base.n_residuals = sum(
        len(experiment["T"]) for experiment in fit_base.experiments
    )

    tds_model.k0 = args.k0
    tds_model.p0 = args.p0

    x0 = fit_base.make_x0(nominal_start)
    lower, upper = fit_base.bounds()

    print("Sensitivity fit: four-state model")
    print(f"k0 = {args.k0:.8e} m3/s")
    print(f"p0 = {args.p0:.8e} s^-1")
    print()

    context = mp.get_context("spawn")

    with ProcessPoolExecutor(
        max_workers=fit_base.n_workers,
        mp_context=context,
        initializer=init_worker,
        initargs=(
            fit_base.experiments,
            fit_base.n_residuals,
            args.k0,
            args.p0,
        ),
    ) as executor:
        result = least_squares(
            fit_base.residual,
            x0,
            bounds=(lower, upper),
            method="trf",
            jac="3-point",
            diff_step=fit_base.diff_step(x0),
            x_scale="jac",
            loss="linear",
            ftol=2.0e-4,
            xtol=2.0e-4,
            gtol=2.0e-4,
            max_nfev=200,
            verbose=2,
            workers=executor.map,
        )

    objective, errors = fit_base.evaluate(result.x)
    fractions, energies, Kr0, Er = fit_base.parameters(result.x)

    order = np.argsort(energies)
    fractions = fractions[order]
    energies = energies[order]

    print()
    print("k0        =", f"{args.k0:.8e} m3/s")
    print("p0        =", f"{args.p0:.8e} s^-1")
    print("fractions =", np.array2string(fractions, precision=8))
    print("energies  =", np.array2string(energies, precision=8))
    print("dE        =", np.array2string(np.diff(energies), precision=8))
    print(f"Kr0       = {Kr0:.8e} m4/s")
    print(f"Er        = {Er:.8f} eV")
    print(f"objective = {objective:.8e}")
    print(f"SciPy cost= {0.5 * objective:.8e}")
    print(f"mean L1   = {np.mean(errors):.8f}")
    print("L1        =", np.array2string(errors, precision=6))
    print(f"nfev={result.nfev}, njev={result.njev}")
    print(result.message)


if __name__ == "__main__":
    mp.freeze_support()
    main()
