from pathlib import Path
import matplotlib.pyplot as plt

from model_data import activation_energies, enthalpies, metals

RESULTS_DIR = Path(__file__).resolve().parent / "results" / "reaction_coordinates"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MECHANISM_TITLES = {
    "redox": "Redox",
    "carboxyl": "Carboxyl",
    "formate": "Formate",
}

MECHANISM_STYLES = {
    "redox": {"color": "black", "label": "Redox mechanism"},
    "carboxyl": {"color": "#00A9D6", "label": "Carboxyl mechanism"},
    "formate": {"color": "#2E8B57", "label": "Formate mechanism"},
}


def add_step(labels, energies, metal, ts_label, product_label, ea_key, enthalpy_key):
    current_energy = energies[-1]
    labels.append(ts_label)
    energies.append(current_energy + activation_energies[ea_key][metal])
    labels.append(product_label)
    energies.append(current_energy + enthalpies[enthalpy_key][metal])


def build_redox_profile(metal):
    labels = ["CO*+H2O*"]
    energies = [0.0]

    add_step(
        labels,
        energies,
        metal,
        "TS1",
        "CO*+OH*+H*",
        "Ea1: CO* + H2O* -> CO* + OH* + H*",
        "ΔH1: CO* + H2O* -> CO* + OH* + H*",
    )

    add_step(
        labels,
        energies,
        metal,
        "TS2",
        "CO*+O*+2H*",
        "Ea2: CO* + OH* + H* -> CO* + O* + 2H*",
        "ΔH2: CO* + OH* + H* -> CO* + O* + 2H*",
    )

    add_step(
        labels,
        energies,
        metal,
        "TS3",
        "CO2+2H*",
        "Ea3: CO* + O* + 2H* -> CO2(g) + 2H*",
        "ΔH3: CO* + O* + 2H* -> CO2(g) + 2H*",
    )

    return labels, energies


def build_carboxyl_profile(metal):
    labels = ["CO*+H2O*"]
    energies = [0.0]

    add_step(
        labels,
        energies,
        metal,
        "TS1",
        "CO*+OH*+H*",
        "Ea1: CO* + H2O* -> CO* + OH* + H*",
        "ΔH1: CO* + H2O* -> CO* + OH* + H*",
    )

    add_step(
        labels,
        energies,
        metal,
        "TS4",
        "COOH*+H*",
        "Ea4: CO* + OH* + H* -> COOH* + H*",
        "ΔH4: CO* + OH* + H* -> COOH* + H*",
    )

    add_step(
        labels,
        energies,
        metal,
        "TS5",
        "CO2+2H*",
        "Ea5: COOH* + H* -> CO2(g) + 2H*",
        "ΔH5: COOH* + H* -> CO2(g) + 2H*",
    )

    return labels, energies


def build_formate_profile(metal):
    labels = ["CO*+H2O*"]
    energies = [0.0]

    add_step(
        labels,
        energies,
        metal,
        "TS1",
        "CO*+OH*+H*",
        "Ea1: CO* + H2O* -> CO* + OH* + H*",
        "ΔH1: CO* + H2O* -> CO* + OH* + H*",
    )

    add_step(
        labels,
        energies,
        metal,
        "TS2",
        "CO*+O*+2H*",
        "Ea2: CO* + OH* + H* -> CO* + O* + 2H*",
        "ΔH2: CO* + OH* + H* -> CO* + O* + 2H*",
    )

    add_step(
        labels,
        energies,
        metal,
        "TS6",
        "CHO*+O*+H*",
        "Ea6: CO* + O* + 2H* -> CHO* + O* + H*",
        "ΔH6: CO* + O* + 2H* -> CHO* + O* + H*",
    )

    add_step(
        labels,
        energies,
        metal,
        "TS7",
        "HCOO*+H*",
        "Ea7: CHO* + O* + H* -> HCOO* + H*",
        "ΔH7: CHO* + O* + H* -> HCOO* + H*",
    )

    add_step(
        labels,
        energies,
        metal,
        "TS8",
        "CO2+2H*",
        "Ea8: HCOO* + H* -> CO2(g) + 2H*",
        "ΔH8: HCOO* + H* -> CO2(g) + 2H*",
    )

    return labels, energies


PROFILE_BUILDERS = {
    "redox": build_redox_profile,
    "carboxyl": build_carboxyl_profile,
    "formate": build_formate_profile,
}


def draw_plateau_profile(ax, labels, energies, color, label, y_offset=0.0):
    bar_half_width = 0.24
    xs = list(range(len(labels)))

    for i, (x, state_label, energy) in enumerate(zip(xs, labels, energies)):
        display_energy = energy + y_offset
        ax.hlines(
            display_energy,
            x - bar_half_width,
            x + bar_half_width,
            color=color,
            linewidth=3.2,
            label=label if i == 0 else None,
        )

        if i < len(labels) - 1:
            next_energy = energies[i + 1] + y_offset
            ax.plot(
                [x + bar_half_width, xs[i + 1] - bar_half_width],
                [display_energy, next_energy],
                color=color,
                linestyle="--",
                linewidth=1.1,
                alpha=0.85,
            )

        text_color = color if state_label.startswith("TS") else "black"
        vertical_offset = 0.08 if state_label.startswith("TS") else -0.12
        vertical_alignment = "bottom" if state_label.startswith("TS") else "top"
        ax.text(
            x,
            display_energy + vertical_offset,
            state_label,
            ha="center",
            va=vertical_alignment,
            fontsize=8,
            color=text_color,
            fontweight="bold" if state_label.startswith("TS") else "normal",
        )

        if i > 0:
            previous_energy = energies[i - 1]
            delta_energy = energy - previous_energy
            midpoint_x = x - 0.5
            midpoint_y = (energy + previous_energy) / 2 + y_offset
            ax.text(
                midpoint_x,
                midpoint_y,
                f"{delta_energy:+.2f}",
                color=color,
                fontsize=8,
                ha="center",
                va="center",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.0),
            )


def plot_metal_profiles(metal):
    fig, ax = plt.subplots(figsize=(13, 7))
    all_energies = []

    for mechanism_name, build_profile_fn in PROFILE_BUILDERS.items():
        style = MECHANISM_STYLES[mechanism_name]
        labels, energies = build_profile_fn(metal)
        all_energies.extend(energies)
        draw_plateau_profile(
            ax,
            labels,
            energies,
            color=style["color"],
            label=style["label"],
        )

    min_energy = min(all_energies)
    max_energy = max(all_energies)
    y_padding = max(0.6, 0.12 * (max_energy - min_energy))

    ax.axhline(0, color="black", linestyle=(0, (3, 3)), linewidth=1.0, alpha=0.75)
    ax.set_xlim(-0.6, 10.6)
    ax.set_ylim(min_energy - y_padding, max_energy + y_padding)
    ax.set_xticks([])
    ax.set_xlabel("Reaction coordinate", fontsize=12)
    ax.set_ylabel("Relative energy from enthalpy and Ea (eV)")
    ax.set_title(f"Activation Energy Reaction Coordinate ({metal})")
    ax.grid(False)
    ax.legend(loc="upper right", frameon=False, fontsize=11)

    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.annotate(
        "",
        xy=(1.01, 0),
        xycoords=("axes fraction", "axes fraction"),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="-|>", linewidth=1.5, color="black"),
        annotation_clip=False,
    )
    ax.annotate(
        "",
        xy=(0, 1.02),
        xycoords=("axes fraction", "axes fraction"),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="-|>", linewidth=1.5, color="black"),
        annotation_clip=False,
    )

    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"reaction_coordinate_{metal}.png", dpi=300, bbox_inches="tight")
    return fig


def main():
    for plot_file in RESULTS_DIR.glob("reaction_coordinate_*.png"):
        plot_file.unlink()

    for metal in metals:
        plot_metal_profiles(metal)

    plt.show()


if __name__ == "__main__":
    main()
