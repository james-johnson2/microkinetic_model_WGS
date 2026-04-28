from pathlib import Path
import csv
import matplotlib.pyplot as plt

ROOT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = ROOT_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)
TOF_CSV = ROOT_DIR / "tof_results.csv"


def read_tof_results(csv_path):
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return [
            {
                "metal": row["metal"],
                "mechanism": row["mechanism"],
                "E_O": float(row["dE_O"]),
                "E_CO": float(row["dE_CO"]),
                "log10_tof": float(row["log10_tof"]),
                "tof": float(row["tof"]),
            }
            for row in reader
        ]


def plot_tof_vs_descriptor(results, descriptor, output_name):
    fig, ax = plt.subplots(figsize=(8, 6))
    for mechanism, marker, color in [
        ("redox", "o", "tab:blue"),
        ("carboxyl", "s", "tab:orange"),
        ("formate", "^", "tab:green"),
    ]:
        subset = [row for row in results if row["mechanism"] == mechanism]
        ax.scatter(
            [row[descriptor] for row in subset],
            [row["tof"] for row in subset],
            label=mechanism.capitalize(),
            marker=marker,
            color=color,
            edgecolor="black",
            s=70,
            alpha=0.85,
        )
        for row in subset:
            ax.annotate(
                row["metal"],
                (row[descriptor], row["tof"]),
                textcoords="offset points",
                xytext=(4, 4),
                fontsize=8,
                color=color,
                alpha=0.9,
            )

    ax.set_xlabel("E_O / eV" if descriptor == "E_O" else "E_CO / eV")
    ax.set_ylabel("TOF / s^-1")
    ax.set_yscale("log")
    ax.set_title("TOF vs " + ("E_O" if descriptor == "E_O" else "E_CO"))
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.legend(title="Mechanism")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / output_name, dpi=300, bbox_inches="tight")
    return fig


def main():
    results = read_tof_results(TOF_CSV)
    plot_tof_vs_descriptor(results, "E_O", "tof_vs_E_O.png")
    plot_tof_vs_descriptor(results, "E_CO", "tof_vs_E_CO.png")
    plt.show()


if __name__ == "__main__":
    main()
