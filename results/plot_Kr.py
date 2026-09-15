from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


k_B = 8.617333262e-5

results = Path(__file__).resolve().parent
data = np.genfromtxt(
    results / "sensitivity.csv",
    delimiter=",",
    names=True,
)

T = np.linspace(300.0, 1100.0, 2000)

fig, ax = plt.subplots(figsize=(7.0, 5.0))

for row in data:
    Kr = row["Kr0"] * np.exp(-row["Er"] / (k_B * T))

    label = (
        rf"$D_0/D_{{0,\mathrm{{nom}}}}={row['D0_factor']:g}$, "
        rf"$\Delta E_D={row['ED_shift_eV']:+.1f}$ eV, "
        rf"$p_0/p_{{0,\mathrm{{nom}}}}={row['p0_factor']:g}$"
    )

    ax.plot(T, Kr, linewidth=1.2, label=label)

ax.set_yscale("log")
ax.set_xlim(300.0, 1100.0)
ax.set_xlabel("Temperature, K")
ax.set_ylabel(r"$K_\mathrm{r}$, m$^4$ s$^{-1}$")
ax.legend(frameon=False, fontsize=7)

fig.tight_layout()
fig.savefig(
    results / "Kr_dependencies.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
