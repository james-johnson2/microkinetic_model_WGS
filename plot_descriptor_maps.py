from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


ROOT_DIR = Path(__file__).resolve().parent
TOF_CSV = ROOT_DIR / "tof_results.csv"
if not TOF_CSV.exists():
    TOF_CSV = ROOT_DIR / "results" / "tof_results.csv"

OUTPUT_DIR = ROOT_DIR / "results" / "descriptor_maps"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MECHANISM_COLORS = {
    "redox": "black",
    "carboxyl": "#1f77b4",
    "formate": "#2ca02c",
}


def read_tof_results(csv_path: Path) -> list[dict[str, float | str]]:
    rows = []
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "metal": row["metal"],
                    "mechanism": row["mechanism"],
                    "dE_CO": float(row["dE_CO"]),
                    "dE_O": float(row["dE_O"]),
                    "log10_tof": float(row["log10_tof"]),
                }
            )
    return rows


def plot_mechanism_maps(rows: list[dict[str, float | str]]) -> None:
    mechanisms = sorted({str(row["mechanism"]) for row in rows})
    vmin = min(float(row["log10_tof"]) for row in rows)
    vmax = max(float(row["log10_tof"]) for row in rows)

    for mechanism in mechanisms:
        mechanism_rows = [row for row in rows if row["mechanism"] == mechanism]

        fig, ax = plt.subplots(figsize=(7, 5.5))
        scatter = ax.scatter(
            [float(row["dE_CO"]) for row in mechanism_rows],
            [float(row["dE_O"]) for row in mechanism_rows],
            c=[float(row["log10_tof"]) for row in mechanism_rows],
            s=120,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            edgecolor="black",
            linewidth=0.7,
        )

        for row in mechanism_rows:
            ax.annotate(
                str(row["metal"]),
                (float(row["dE_CO"]), float(row["dE_O"])),
                textcoords="offset points",
                xytext=(6, 5),
                fontsize=9,
            )

        colorbar = fig.colorbar(scatter, ax=ax)
        colorbar.set_label(r"$\log_{10}(\mathrm{TOF})$")

        ax.set_xlabel(r"$\Delta E_{\mathrm{CO}}\ \mathrm{(eV)}$")
        ax.set_ylabel(r"$\Delta E_{\mathrm{O}}\ \mathrm{(eV)}$")
        ax.set_title(f"{mechanism.capitalize()} Descriptor Map")
        fig.tight_layout()
        fig.savefig(OUTPUT_DIR / f"{mechanism}_descriptor_map.png", dpi=300)
        plt.close(fig)


def scale_point_sizes(values: list[float], min_size: float = 90.0, max_size: float = 360.0) -> list[float]:
    min_value = min(values)
    max_value = max(values)

    if max_value == min_value:
        return [(min_size + max_size) / 2 for _ in values]

    return [
        min_size + (value - min_value) * (max_size - min_size) / (max_value - min_value)
        for value in values
    ]


def plot_dominant_mechanism_map(rows: list[dict[str, float | str]]) -> None:
    rows_by_metal: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        rows_by_metal[str(row["metal"])].append(row)

    dominant_rows = []
    for metal, metal_rows in rows_by_metal.items():
        dominant_row = max(metal_rows, key=lambda row: float(row["log10_tof"]))
        total_tof = sum(10 ** float(row["log10_tof"]) for row in metal_rows)
        dominant_rows.append(
            {
                "metal": metal,
                "dominant_mechanism": str(dominant_row["mechanism"]),
                "dE_CO": float(dominant_row["dE_CO"]),
                "dE_O": float(dominant_row["dE_O"]),
                "log10_total_tof": math.log10(total_tof),
            }
        )

    sizes = scale_point_sizes([row["log10_total_tof"] for row in dominant_rows])

    fig, ax = plt.subplots(figsize=(7, 5.5))
    for row, size in zip(dominant_rows, sizes):
        mechanism = row["dominant_mechanism"]
        ax.scatter(
            row["dE_CO"],
            row["dE_O"],
            s=size,
            color=MECHANISM_COLORS.get(mechanism, "gray"),
            edgecolor="black",
            linewidth=0.7,
            alpha=0.9,
        )
        ax.annotate(
            row["metal"],
            (row["dE_CO"], row["dE_O"]),
            textcoords="offset points",
            xytext=(6, 5),
            fontsize=9,
        )

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgecolor="black",
            markersize=9,
            label=mechanism.capitalize(),
        )
        for mechanism, color in MECHANISM_COLORS.items()
    ]
    ax.legend(handles=handles, title="Dominant mechanism")

    ax.set_xlabel(r"$\Delta E_{\mathrm{CO}}\ \mathrm{(eV)}$")
    ax.set_ylabel(r"$\Delta E_{\mathrm{O}}\ \mathrm{(eV)}$")
    ax.set_title("Dominant Mechanism Descriptor Map")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "dominant_mechanism_map.png", dpi=300)
    plt.close(fig)


def main() -> None:
    rows = read_tof_results(TOF_CSV)
    plot_mechanism_maps(rows)
    plot_dominant_mechanism_map(rows)
    print(f"Saved descriptor maps to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
