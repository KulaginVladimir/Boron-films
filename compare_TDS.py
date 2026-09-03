import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fit_TDS.tds_model import T0, n_B, simulate_TDS


area = 1.0e-4
Tf_fit = 1100.0
Tf_plot = 1300.0

fractions = np.array([0.19070939, 0.54969147, 0.16451335, 0.09508579])
energies = np.array([0.92443877, 1.13050166, 1.32677295, 1.62738778])
Kr0 = 2.69456812e-25
E_r = 0.84338380

base = Path(__file__).resolve().parent
experimental_data = base / "experimental_data"
experiments_file = experimental_data / "experiments_5.csv"
calculated_spectra = base / "calculated_spectra"

parser = argparse.ArgumentParser()
parser.add_argument("--recalculate", action="store_true")
args = parser.parse_args()


def read_tds(path):
    data = np.genfromtxt(
        path,
        delimiter=",",
        skip_header=10,
        missing_values="--",
        filling_values=np.nan,
        invalid_raise=False,
    )
    return data[np.all(np.isfinite(data[:, :4]), axis=1)]


def get_inventory(data, thickness):
    fit = data[(data[:, 1] >= T0) & (data[:, 1] <= Tf_fit)]

    order = np.argsort(fit[:, 0])
    time = fit[order, 0]
    flux = (fit[order, 2] + fit[order, 3]) / area

    return np.trapezoid(flux, time) / (thickness * n_B)


plt.rcParams.update(
    {
        "font.size": 8.5,
        "axes.titlesize": 9.0,
        "axes.labelsize": 9.0,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
    }
)

with open(experiments_file, newline="") as file:
    experiments = list(csv.DictReader(file))

calculated_spectra.mkdir(exist_ok=True)

fig, axes = plt.subplots(
    2,
    3,
    figsize=(11.2, 6.5),
    sharex=True,
)
axes = axes.ravel()

for i, (ax, experiment) in enumerate(zip(axes[:5], experiments)):
    experiment_id = experiment["id"]
    beta = float(experiment["heating_rate_K_per_s"])
    thickness = float(experiment["film_thickness_nm"]) * 1e-9

    data = read_tds(experimental_data / f"{experiment_id}.csv")

    plot = data[(data[:, 1] >= T0) & (data[:, 1] <= Tf_plot)]
    order = np.argsort(plot[:, 1])

    T_exp = plot[order, 1]
    J_exp = (plot[order, 2] + plot[order, 3]) / area

    spectrum_file = calculated_spectra / f"{experiment_id}.csv"

    if spectrum_file.exists() and not args.recalculate:
        calculated = np.genfromtxt(
            spectrum_file,
            delimiter=",",
            names=True,
        )

        T_calc = calculated["T_K"]
        J_calc = calculated["J_calc"]
        release = [calculated[f"J_trap{j + 1}"] for j in range(4)]
    else:
        inventory = get_inventory(data, thickness)

        T_model, J_model, release_model = simulate_TDS(
            beta=beta,
            thickness=thickness,
            inventory=inventory,
            fractions=fractions,
            energies=energies,
            Kr0=Kr0,
            E_r=E_r,
            Tf=Tf_plot,
            return_trap_contribution=True,
        )

        T_calc = T_exp
        J_calc = np.interp(T_exp, T_model, J_model)

        release = [
            np.interp(T_exp, T_model[1:], contribution)
            for contribution in release_model
        ]

        np.savetxt(
            spectrum_file,
            np.column_stack((T_exp, J_exp, J_calc, *release)),
            delimiter=",",
            header="T_K,J_exp,J_calc,J_trap1,J_trap2,J_trap3,J_trap4",
            comments="",
        )

    mask_fit = (T_exp >= T0) & (T_exp <= Tf_fit)

    error = np.trapezoid(
        np.abs(J_calc[mask_fit] - J_exp[mask_fit]), T_exp[mask_fit]
    ) / np.trapezoid(np.abs(J_exp[mask_fit]), T_exp[mask_fit])

    ax.plot(
        T_calc,
        J_calc,
        linewidth=1.9,
        label="calculated",
        zorder=4,
    )
    ax.plot(
        T_exp,
        J_exp,
        linewidth=1.6,
        label="experiment",
        zorder=5,
    )

    for j, contribution in enumerate(release):
        ax.fill_between(
            T_calc,
            0.0,
            contribution,
            alpha=0.22,
            linewidth=0.0,
            label=f"Trap {j + 1}",
            zorder=1,
        )

    ax.set_title(
        f"{experiment_id}\n"
        rf"$\beta={beta:g}$ K s$^{{-1}}$, "
        rf"$L={thickness * 1e9:.0f}$ nm",
        pad=6,
    )

    ax.text(
        0.97,
        0.95,
        rf"$\mathrm{{err}}={error:.3f}$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox=dict(
            boxstyle="round,pad=0.18",
            facecolor="white",
            alpha=0.75,
            edgecolor="none",
        ),
    )

    ax.set_xlim(T0, Tf_plot)
    ax.tick_params(direction="in")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(alpha=0.20, linewidth=0.6)

for ax in axes[5:]:
    ax.axis("off")

for ax in axes[3:5]:
    ax.set_xlabel("Temperature, K")

for ax in axes[::3]:
    if ax.axison:
        ax.set_ylabel(r"Flux, D m$^{-2}$ s$^{-1}$")

handles, labels = axes[0].get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.985),
    ncol=6,
    frameon=False,
    columnspacing=1.3,
    handlelength=2.0,
)

fig.subplots_adjust(
    left=0.075,
    right=0.985,
    bottom=0.10,
    top=0.84,
    wspace=0.28,
    hspace=0.38,
)

fig.savefig(
    base / "tds_comparison.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.04,
)
plt.show()
