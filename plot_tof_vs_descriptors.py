from pathlib import Path
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

ROOT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = ROOT_DIR / "results"
TOF_RESULTS_DIR = RESULTS_DIR / "tof_comparisons"
RESULTS_DIR.mkdir(exist_ok=True)
TOF_RESULTS_DIR.mkdir(exist_ok=True)
TOF_CSV = ROOT_DIR / "tof_results.csv"
if not TOF_CSV.exists():
    TOF_CSV = RESULTS_DIR / "tof_results.csv"

from model_data import metals, exp_tof, ads_energies

markers = {"redox": "o", "carboxyl": "s", "formate": "^"}
mechanism_colors = {"redox": "black", "carboxyl": "#1f77b4", "formate": "#2ca02c"}
colors = cm.tab10(np.linspace(0, 1, len(metals)))
metal_colors = dict(zip(metals, colors))
descriptors = {
    "E_CO": {"species": "CO", "label": "E_CO / eV"},
    "E_O": {"species": "O", "label": "E_O / eV"},
    "E_H": {"species": "H", "label": "E_H / eV"},
    "E_OH": {"species": "OH", "label": "E_OH / eV"},
}


def descriptor_output_dir(descriptor):
    output_dir = TOF_RESULTS_DIR / descriptor
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def descriptor_value(metal, descriptor):
    return ads_energies[metal][descriptors[descriptor]["species"]]


def read_tof_results(csv_path):
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return [
            {
                "metal": row["metal"],
                "mechanism": row["mechanism"],
                "log10_tof": float(row["log10_tof"]),
                "tof": float(row["tof"]),
            }
            for row in reader
        ]


def plot_tof_vs_descriptor(results, descriptor, output_name, mechanism_filter=None, show_segmentation=True):
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot calculated TOF
    for metal in metals:
        color = metal_colors[metal]
        for mechanism in ["redox", "carboxyl", "formate"]:
            if mechanism_filter and mechanism != mechanism_filter:
                continue
            row = next((r for r in results if r["metal"] == metal and r["mechanism"] == mechanism), None)
            if row:
                ax.scatter(
                    descriptor_value(metal, descriptor),
                    row["tof"],
                    marker=markers[mechanism],
                    color=color,
                    edgecolor="black",
                    s=70,
                    alpha=0.85,
                )
    
    # Plot experimental data
    for metal in metals:
        if metal in exp_tof:
            ax.scatter(
                descriptor_value(metal, descriptor),
                exp_tof[metal],
                marker='x',
                color=metal_colors[metal],
                s=70,
                alpha=0.85,
            )
    
    # Annotate metals above model data clusters
    for metal in metals:
        relevant_results = [r for r in results if r["metal"] == metal and (not mechanism_filter or r["mechanism"] == mechanism_filter)]
        if relevant_results:
            max_tof = max(r["tof"] for r in relevant_results)
            x_pos = descriptor_value(metal, descriptor)
            ax.annotate(
                metal,
                (x_pos, max_tof),
                textcoords="offset points",
                xytext=(0, 10),
                fontsize=8,
                color=metal_colors[metal],
                ha='center'
            )
    
    # Highlight regimes for combined E_CO/E_O plots using the current axis window.
    if show_segmentation and not mechanism_filter and descriptor in {"E_CO", "E_O"}:
        x_min, x_max = ax.get_xlim()
        if descriptor == "E_CO":
            regimes = [(x_min, -1.75, "Strong E_CO"), (-1.75, -0.5, "Intermediate E_CO"), (-0.5, x_max, "Weak E_CO")]
        else:
            regimes = [(x_min, -4.5, "Strong E_O"), (-4.5, -3.5, "Intermediate E_O"), (-3.5, x_max, "Weak E_O")]

        for left, right, label in regimes:
            ax.axvspan(left, right, alpha=0.12, color='blue' if "Strong" in label else 'green' if "Intermediate" in label else 'red')
            x_text = left + (right - left) * 0.5
            ax.text(
                x_text,
                0.95,
                label,
                transform=ax.get_xaxis_transform(),
                fontsize=9,
                color='black',
                ha='center',
                va='top',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='none')
            )

    # Create legend without metals
    handles = []
    for mechanism in ["redox", "carboxyl", "formate"]:
        if not mechanism_filter or mechanism == mechanism_filter:
            handles.append(plt.Line2D([0], [0], marker=markers[mechanism], color='black', linestyle='None', markersize=8, label=mechanism.capitalize()))
    handles.append(plt.Line2D([0], [0], marker='x', color='black', linestyle='None', markersize=8, label='Experimental'))
    
    ax.legend(handles=handles, title="Mechanism / Exp")
    
    ax.set_xlabel(descriptors[descriptor]["label"])
    ax.set_ylabel("TOF / s^-1")
    ax.set_yscale("log")
    title = f"TOF vs {descriptor}"
    if mechanism_filter:
        title += f" ({mechanism_filter.capitalize()})"
    ax.set_title(title)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
    fig.tight_layout()
    fig.savefig(descriptor_output_dir(descriptor) / output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_exp_vs_descriptor(descriptor, output_name):
    fig, ax = plt.subplots(figsize=(8, 6))
    
    for metal in metals:
        if metal in exp_tof:
            ax.scatter(
                descriptor_value(metal, descriptor),
                exp_tof[metal],
                marker='o',
                color='black',
                s=70,
                alpha=0.85,
            )
            ax.annotate(
                metal,
                (descriptor_value(metal, descriptor), exp_tof[metal]),
                textcoords="offset points",
                xytext=(0, 10),
                fontsize=8,
                color='black',
                ha='center'
            )
    
    ax.set_xlabel(descriptors[descriptor]["label"])
    ax.set_ylabel("TOF / s^-1")
    ax.set_yscale("log")
    ax.set_title(f"Experimental TOF vs {descriptor}")
    fig.tight_layout()
    fig.savefig(descriptor_output_dir(descriptor) / output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_connected_mechanisms(results, descriptor, output_name):
    fig, ax = plt.subplots(figsize=(8, 6))

    for mechanism in ["redox", "carboxyl", "formate"]:
        mechanism_points = []

        for metal in metals:
            row = next((r for r in results if r["metal"] == metal and r["mechanism"] == mechanism), None)
            if row:
                mechanism_points.append(
                    {
                        "metal": metal,
                        "x": descriptor_value(metal, descriptor),
                        "tof": row["tof"],
                    }
                )

        mechanism_points.sort(key=lambda point: point["x"])

        ax.plot(
            [point["x"] for point in mechanism_points],
            [point["tof"] for point in mechanism_points],
            marker=markers[mechanism],
            color=mechanism_colors[mechanism],
            markeredgecolor="black",
            linewidth=1.6,
            markersize=7,
            label=mechanism.capitalize(),
        )

        for point in mechanism_points:
            ax.annotate(
                point["metal"],
                (point["x"], point["tof"]),
                textcoords="offset points",
                xytext=(0, 8),
                fontsize=8,
                color=mechanism_colors[mechanism],
                ha="center",
            )

    ax.set_xlabel(descriptors[descriptor]["label"])
    ax.set_ylabel("TOF / s^-1")
    ax.set_yscale("log")
    ax.set_title(f"Connected Mechanism TOF vs {descriptor}")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.legend(title="Mechanism")
    fig.tight_layout()
    fig.savefig(descriptor_output_dir(descriptor) / output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig


def main():
    results = read_tof_results(TOF_CSV)

    for descriptor in descriptors:
        plot_tof_vs_descriptor(results, descriptor, f"tof_vs_{descriptor}.png")

        if descriptor in {"E_CO", "E_O"}:
            plot_tof_vs_descriptor(
                results,
                descriptor,
                f"tof_vs_{descriptor}_no_segmentation.png",
                show_segmentation=False,
            )

        for mechanism in ["redox", "carboxyl", "formate"]:
            plot_tof_vs_descriptor(
                results,
                descriptor,
                f"tof_vs_{descriptor}_{mechanism}.png",
                mechanism_filter=mechanism,
            )

        plot_exp_vs_descriptor(descriptor, f"exp_tof_vs_{descriptor}.png")
        plot_connected_mechanisms(results, descriptor, f"tof_vs_{descriptor}_connected.png")

    plt.show()


if __name__ == "__main__":
    main()
