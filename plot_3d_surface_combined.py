from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from volcano_plot import build_surface_grid, redox_surface, carboxyl_surface, formate_surface

RESULTS_DIR = Path(__file__).resolve().parent / "results" / "3d_plots"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def build_combined_surface():
    E_CO_mesh, E_O_mesh, Z_redox = build_surface_grid(redox_surface)
    _, _, Z_carboxyl = build_surface_grid(carboxyl_surface)
    _, _, Z_formate = build_surface_grid(formate_surface)

    stacked = np.stack([Z_redox, Z_carboxyl, Z_formate], axis=0)
    Z_best = np.nanmax(stacked, axis=0)
    return E_CO_mesh, E_O_mesh, Z_best


def plot_combined_surface():
    E_CO_mesh, E_O_mesh, Z_best = build_combined_surface()

    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        E_CO_mesh,
        E_O_mesh,
        Z_best,
        cmap="viridis",
        linewidth=0,
        antialiased=True,
        alpha=0.9,
    )

    ax.set_xlabel("E_CO / eV")
    ax.set_ylabel("E_O / eV")
    ax.set_zlabel("log10(TOF)")
    ax.set_title("Combined Best Mechanism 3D Volcano Surface")
    fig.colorbar(surf, ax=ax, shrink=0.6, label="log10(TOF)")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "combined_best_3d_surface.png", dpi=300, bbox_inches="tight")
    return fig


def main():
    plot_combined_surface()
    plt.show()


if __name__ == "__main__":
    main()
