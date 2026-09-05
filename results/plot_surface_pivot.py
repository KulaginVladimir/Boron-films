from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


k_B = 8.617333262e-5

results = Path(__file__).resolve().parent
sensitivity_file = results / "sensitivity.csv"

data = np.genfromtxt(
    sensitivity_file,
    delimiter=",",
    names=True,
)

T = np.linspace(350.0, 900.0, 2000)

ln_Kr = np.array(
    [np.log(Kr0) - Er / (k_B * T) for Kr0, Er in zip(data["Kr0"], data["Er"])]
)

spread = np.std(ln_Kr, axis=0)

i_min = np.argmin(spread)
T_min = T[i_min]
T_p = 10.0 * np.round(T_min / 10.0)

print(f"Minimum-spread temperature: {T_min:.1f} K")
print(f"Reported pivot temperature: {T_p:.0f} K")
print(f"std(ln K_r) at pivot: {np.interp(T_p, T, spread):.4f}")

fig, ax = plt.subplots(figsize=(7.0, 5.0))

x = 1000.0 / T

for Kr0, Er in zip(data["Kr0"], data["Er"]):
    ax.plot(
        x,
        np.log(Kr0) - Er / (k_B * T),
        linewidth=1.2,
    )

ax.axvline(
    1000.0 / T_p,
    linestyle="--",
    linewidth=1.2,
    label=rf"$T_\mathrm{{p}}={T_p:.0f}$ K",
)

ax.set_xlabel(r"$1000/T$, K$^{-1}$")
ax.set_ylabel(r"$\ln K_\mathrm{r}(T)$")
ax.legend(frameon=False)
ax.grid(alpha=0.2)

fig.tight_layout()
fig.savefig(
    results / "surface_pivot.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
