import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fit_TDS.tds_model import T0, n_B, simulate_TDS


area = 1.0e-4
Tf_fit = 1100.0
Tf_plot = 1300.0

fractions = np.array([0.12654304, 0.16581731, 0.45692251, 0.15936622, 0.09135092])
energies = np.array([0.91480893, 1.05522347, 1.20260982, 1.41705684, 1.73039792])
Kr0 = 1.39328746e-23
E_r = 0.99708420

sink_fractions = np.array([0.13816642, 0.22944871, 0.36937497, 0.15742860, 0.10558130])
sink_energies = np.array([1.28625148, 1.50171091, 1.69711525, 1.99016510, 2.43786185])

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
        J_sink = calculated["J_sink"]
        release = [calculated[f"J_trap{j + 1}"] for j in range(5)]
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
            surface_model="arrhenius",
        )

        T_sink_model, J_sink_model = simulate_TDS(
            beta=beta,
            thickness=thickness,
            inventory=inventory,
            fractions=sink_fractions,
            energies=sink_energies,
            Kr0=0.0,
            E_r=0.0,
            Tf=Tf_plot,
            surface_model="sink",
        )

        T_calc = T_exp
        J_calc = np.interp(T_exp, T_model, J_model)
        J_sink = np.interp(T_exp, T_sink_model, J_sink_model)

        release = [
            np.interp(T_exp, T_model[1:], contribution)
            for contribution in release_model
        ]

        np.savetxt(
            spectrum_file,
            np.column_stack((T_exp, J_exp, J_calc, J_sink, *release)),
            delimiter=",",
            header=("T_K,J_exp,J_calc,J_sink,J_trap1,J_trap2,J_trap3,J_trap4,J_trap5"),
            comments="",
        )

    mask_fit = (T_exp >= T0) & (T_exp <= Tf_fit)

    error = np.trapezoid(
        np.abs(J_calc[mask_fit] - J_exp[mask_fit]), T_exp[mask_fit]
    ) / np.trapezoid(np.abs(J_exp[mask_fit]), T_exp[mask_fit])

    sink_error = np.trapezoid(
        np.abs(J_sink[mask_fit] - J_exp[mask_fit]), T_exp[mask_fit]
    ) / np.trapezoid(np.abs(J_exp[mask_fit]), T_exp[mask_fit])

    ax.plot(
        T_calc,
        J_calc,
        linewidth=1.5,
        label="Finite-rate desorption",
        zorder=4,
    )
    ax.plot(
        T_exp,
        J_exp,
        linewidth=3,
        label="Experiment",
        color="tab:red",
        alpha=0.65,
        zorder=0,
    )
    ax.plot(
        T_calc,
        J_sink,
        linewidth=1.5,
        linestyle="--",
        label="Instantaneous desorption",
        color="black",
        zorder=3,
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
        rf"$\mathrm{{err}}={error:.3f}$"
        "\n"
        rf"$\mathrm{{err}}_{{\rm sink}}={sink_error:.3f}$",
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
    ncol=4,
    frameon=False,
    columnspacing=1.3,
    handlelength=2.0,
)

fig.subplots_adjust(
    left=0.075,
    right=0.985,
    bottom=0.10,
    top=0.80,
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
