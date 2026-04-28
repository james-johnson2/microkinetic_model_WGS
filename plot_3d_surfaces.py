from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from model_data import ads_energies, metals
from volcano_plot import build_surface_grid, redox_surface, carboxyl_surface, formate_surface

RESULTS_DIR = Path(__file__).resolve().parent / "results" / "3d_plots"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

mechanisms = [
    ("redox", redox_surface),
    ("carboxyl", carboxyl_surface),
    ("formate", formate_surface),
]


def plot_3d_surface(name, surface_fn):
    E_CO_mesh, E_O_mesh, Z = build_surface_grid(surface_fn)

    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        E_CO_mesh,
        E_O_mesh,
        Z,
        cmap="viridis",
        linewidth=0,
        antialiased=True,
        alpha=0.9,
    )

    for metal in metals:
        x = ads_energies[metal]["CO"]
        y = ads_energies[metal]["O"]
        z = surface_fn(x, y)
        ax.scatter(x, y, z, color="red", s=30)
        ax.text(x, y, z + 0.2, metal, color="red", fontsize=9)

    ax.set_xlabel("E_CO / eV")
    ax.set_ylabel("E_O / eV")
    ax.set_zlabel("log10(TOF)")
    ax.set_title(f"{name.capitalize()} 3D Volcano Surface")
    fig.colorbar(surf, ax=ax, shrink=0.6, label="log10(TOF)")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"{name}_3d_surface.png", dpi=300, bbox_inches="tight")
    return fig


def main():
    for name, fn in mechanisms:
        plot_3d_surface(name, fn)
    plt.show()


if __name__ == "__main__":
    main()
