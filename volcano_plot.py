from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from model_data import (
    VIB_DATA,
    activation_energies,
    ads_energies,
    enthalpies,
    h_planck,
    k_B,
    metals,
    p_CO,
    p_CO2,
    p_H2,
    p_H2O,
    T,
)
from TOF_calc import K, TOF_formate, carb_TOF, dG, wgsr_redox_tof


RESULTS_DIR = Path(__file__).resolve().parent / "results" / "volcano_heatmaps"

# Use the same metal set as the model data for fitting descriptor scaling lines
fit_metals = metals

# The volcano surfaces use CO and O adsorption energies as the two descriptors
X_desc = np.array([
    [ads_energies[m]["CO"], ads_energies[m]["O"]]
    for m in fit_metals
])


def fit_descriptor_model(values_by_metal):
    #Fit a linear scaling relation against E_CO and E_O
    y = np.array([values_by_metal[m] for m in fit_metals])
    model = LinearRegression()
    model.fit(X_desc, y)
    return model


def predict(model, E_CO, E_O):
    #Predict a scaled thermodynamic or kinetic value at one descriptor point
    return float(model.predict(np.array([[E_CO, E_O]]))[0])


# H2O adsorption is not a descriptor, so it is predicted from E_CO and E_O
ads_models = {
    "H2O": fit_descriptor_model({
        m: ads_energies[m]["H2O"] for m in fit_metals
    })
}

# Scale every tabulated reaction enthalpy and activation energy over descriptor space
enthalpy_models = {
    key: fit_descriptor_model({m: vals[m] for m in fit_metals})
    for key, vals in enthalpies.items()
}

ea_models = {
    key: fit_descriptor_model({m: vals[m] for m in fit_metals})
    for key, vals in activation_energies.items()
}


def scaled_h(reaction, E_CO, E_O):
    #Predict one reaction enthalpy at a descriptor point
    return predict(enthalpy_models[reaction], E_CO, E_O)


def scaled_ea(reaction, E_CO, E_O):
    #Predict one activation energy at a descriptor point
    return predict(ea_models[reaction], E_CO, E_O)


def safe_log10_tof(tof):
    #Return log10(TOF), using NaN for invalid surface points
    if tof <= 0 or not np.isfinite(tof):
        return np.nan
    return np.log10(tof)


def redox_surface(E_CO, E_O):
    #Continuous redox-mechanism log10(TOF) over E_CO/E_O descriptor space
    E_ads_H2O = predict(ads_models["H2O"], E_CO, E_O)

    K1 = K(dG(0, [VIB_DATA["CO_NIST"]], [VIB_DATA["CO"]], E_CO))
    K2 = K(dG(0, [VIB_DATA["H2O_NIST"]], [VIB_DATA["H2O"]], E_ads_H2O))
    K3 = K(dG(
        scaled_h("ΔH1: CO* + H2O* -> CO* + OH* + H*", E_CO, E_O),
        [VIB_DATA["H2O"], VIB_DATA["CO"]],
        [VIB_DATA["OH"], VIB_DATA["H"], VIB_DATA["CO"]]
    ))
    K4 = K(dG(
        scaled_h("ΔH2: CO* + OH* + H* -> CO* + O* + 2H*", E_CO, E_O),
        [VIB_DATA["OH"], VIB_DATA["H"], VIB_DATA["CO"], VIB_DATA["H"]],
        [VIB_DATA["O"], VIB_DATA["H"], VIB_DATA["CO"], VIB_DATA["H"]]
    ))
    K6 = K(dG(
        scaled_h("ΔH13: CO2(g) + 2H* -> CO2(g) + H2(g)", E_CO, E_O),
        [VIB_DATA["H"], VIB_DATA["H"]],
        [VIB_DATA["H2_NIST"]]
    ))

    k5 = (k_B * T / h_planck) * np.exp(
        -scaled_ea("Ea3: CO* + O* + 2H* -> CO2(g) + 2H*", E_CO, E_O)
        / (k_B * T)
    )

    denominator = (
        1
        + K1 * p_CO
        + K2 * p_H2O
        + K2 * K3 * p_H2O * (K6 / p_H2)**0.5
        + (K2 * K3 * K4 * K6 * p_H2O) / p_H2
        + (p_H2 / K6)**0.5
    )

    tof = k5 * K1 * K2 * K3 * K4 * K6 * (p_CO * p_H2O / p_H2) / denominator**2
    return safe_log10_tof(tof)


def carboxyl_surface(E_CO, E_O):
    #Continuous carboxyl-mechanism log10(TOF) over E_CO/E_O descriptor space
    E_ads_H2O = predict(ads_models["H2O"], E_CO, E_O)

    dG1 = dG(0, [VIB_DATA["H2O_NIST"]], [VIB_DATA["H2O"]], E_ads_H2O)
    dG2 = dG(0, [VIB_DATA["CO_NIST"]], [VIB_DATA["CO"]], E_CO)
    dG3 = dG(
        scaled_h("ΔH1: CO* + H2O* -> CO* + OH* + H*", E_CO, E_O),
        [VIB_DATA["H2O"]],
        [VIB_DATA["OH"], VIB_DATA["H"]]
    )
    dG5 = dG(
        scaled_h("ΔH5: COOH* + H* -> CO2(g) + 2H*", E_CO, E_O),
        [VIB_DATA["COOH"]],
        [VIB_DATA["CO2_NIST"], VIB_DATA["H"]]
    )
    dG6 = dG(
        scaled_h("ΔH13: CO2(g) + 2H* -> CO2(g) + H2(g)", E_CO, E_O),
        [VIB_DATA["H"]],
        [VIB_DATA["H2_NIST"]]
    )

    K1 = K(dG1)
    K2 = K(dG2)
    K3 = K(dG3)
    K5 = K(dG5)
    K6 = K(dG6)

    k4 = (k_B * T / h_planck) * np.exp(
        -scaled_ea("Ea4: CO* + OH* + H* -> COOH* + H*", E_CO, E_O)
        / (k_B * T)
    )

    theta_star = (
        1
        + K1 * p_H2O
        + K2 * p_CO
        + K1 * K3 * p_H2O * np.sqrt(K6 / p_H2)
        + (p_CO2 / K5) * np.sqrt(p_H2 / K6)
        + np.sqrt(p_H2 / K6)
    ) ** -1

    tof = k4 * K1 * K2 * K3 * p_CO * p_H2O * np.sqrt(K6 / p_H2) * theta_star**2
    return safe_log10_tof(tof)


def formate_surface(E_CO, E_O):
    #Continuous formate-mechanism log10(TOF) over E_CO/E_O descriptor space
    E_ads_H2O = predict(ads_models["H2O"], E_CO, E_O)

    dG1 = dG(0, [VIB_DATA["CO_NIST"]], [VIB_DATA["CO"]], E_CO)
    dG2 = dG(0, [VIB_DATA["H2O_NIST"]], [VIB_DATA["H2O"]], E_ads_H2O)
    dG3 = dG(
        scaled_h("ΔH1: CO* + H2O* -> CO* + OH* + H*", E_CO, E_O),
        [VIB_DATA["H2O"]],
        [VIB_DATA["OH"], VIB_DATA["H"]]
    )
    dG4 = dG(
        scaled_h("ΔH2: CO* + OH* + H* -> CO* + O* + 2H*", E_CO, E_O),
        [VIB_DATA["OH"]],
        [VIB_DATA["O"], VIB_DATA["H"]]
    )
    dG5 = dG(
        scaled_h("ΔH6: CO* + O* + 2H* -> CHO* + O* + H*", E_CO, E_O),
        [VIB_DATA["CO"], VIB_DATA["H"]],
        [VIB_DATA["CHO"]]
    )
    dG6 = dG(
        scaled_h("ΔH7: CHO* + O* + H* -> HCOO* + H*", E_CO, E_O),
        [VIB_DATA["CHO"], VIB_DATA["O"]],
        [VIB_DATA["HCOO"]]
    )
    dG7 = dG(
        scaled_h("ΔH9: CHO* + OH* -> HCOOH*", E_CO, E_O),
        [VIB_DATA["CHO"], VIB_DATA["OH"]],
        [VIB_DATA["HCOOH"]]
    )
    dG8 = dG(
        scaled_h("ΔH10: HCOOH* -> HCOO* + H*", E_CO, E_O),
        [VIB_DATA["HCOOH"]],
        [VIB_DATA["HCOO"], VIB_DATA["H"]]
    )
    dG10 = dG(
        scaled_h("ΔH13: CO2(g) + 2H* -> CO2(g) + H2(g)", E_CO, E_O),
        [VIB_DATA["H"], VIB_DATA["H"]],
        [VIB_DATA["H2_NIST"]]
    )

    K1 = K(dG1)
    K2 = K(dG2)
    K3 = K(dG3)
    K4 = K(dG4)
    K5 = K(dG5)
    K6 = K(dG6)
    K7 = K(dG7)
    K8 = K(dG8)
    K10 = K(dG10)

    k9 = (k_B * T / h_planck) * np.exp(
        -scaled_ea("Ea8: HCOO* + H* -> CO2(g) + 2H*", E_CO, E_O)
        / (k_B * T)
    )

    alpha = np.sqrt(p_H2 / K10)
    A = (
        1
        + K1 * p_CO
        + K2 * p_H2O
        + alpha
        + (K2 * K3 * p_H2O) / alpha
        + (K2 * K3 * K4 * p_H2O) / alpha**2
        + K1 * K5 * p_CO * alpha
    )
    B = K1 * K2 * K3 * K4 * K5 * K6 * p_CO * p_H2O * np.sqrt(K10 / p_H2)
    C = K1 * K2 * K3 * K5 * K7 * p_CO * p_H2O
    D = 2 * (B + C)

    theta = 2 / (A + np.sqrt(A**2 + 4 * D))
    tof = k9 * B * theta**2
    return safe_log10_tof(tof)


# Each entry defines how to calculate the continuous surface and actual metal points.
mechanisms = {
    "redox": {
        "title": "Redox Mechanism",
        "surface_fn": redox_surface,
        "metal_tof_fn": wgsr_redox_tof,
    },
    "carboxyl": {
        "title": "Carboxyl Mechanism",
        "surface_fn": carboxyl_surface,
        "metal_tof_fn": carb_TOF,
    },
    "formate": {
        "title": "Formate Mechanism",
        "surface_fn": formate_surface,
        "metal_tof_fn": TOF_formate,
    },
}


def build_surface_grid(surface_fn):
    #Evaluate one mechanism over a regular E_CO/E_O descriptor grid
    E_CO_vals = np.linspace(-2.4, 0.2, 100)
    E_O_vals = np.linspace(-7.0, -2.0, 100)

    E_CO_mesh, E_O_mesh = np.meshgrid(E_CO_vals, E_O_vals)
    Z = np.full_like(E_CO_mesh, np.nan, dtype=float)

    for i in range(E_O_mesh.shape[0]):
        for j in range(E_CO_mesh.shape[1]):
            Z[i, j] = surface_fn(E_CO_mesh[i, j], E_O_mesh[i, j])

    # Clipping keeps extremely low/high values from washing out the color scale.
    return E_CO_mesh, E_O_mesh, np.clip(Z, -20, 5)


def plot_volcano_surface(mechanism_name, E_CO_mesh, E_O_mesh, Z_plot):
    #Create one 3D volcano surface and overlay the tabulated metal points
    mechanism = mechanisms[mechanism_name]
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="3d")

    surf = ax.plot_surface(
        E_CO_mesh,
        E_O_mesh,
        Z_plot,
        cmap="jet",
        linewidth=0,
        antialiased=True,
        alpha=0.85
    )

    for metal in fit_metals:
        x = ads_energies[metal]["CO"]
        y = ads_energies[metal]["O"]
        z = mechanism["metal_tof_fn"](metal)

        ax.scatter(x, y, z, color="red", s=40)
        ax.text(x, y, z + 0.35, metal, color="red", fontsize=10)

    ax.set_xlabel("E_CO / eV")
    ax.set_ylabel("E_O / eV")
    ax.set_zlabel("log10(TOF)")
    ax.set_title(f"{mechanism['title']} Volcano Surface")

    fig.colorbar(surf, ax=ax, shrink=0.6, label="log10(TOF)")
    plt.tight_layout()


def plot_volcano_heatmap(mechanism_name, E_CO_mesh, E_O_mesh, Z_plot):
    """Create one 2D contour volcano map and overlay the tabulated metal points."""
    mechanism = mechanisms[mechanism_name]
    fig, ax = plt.subplots(figsize=(8, 6))
    contour = ax.contourf(E_CO_mesh, E_O_mesh, Z_plot, levels=40, cmap="jet")
    fig.colorbar(contour, ax=ax, label="log10(TOF)")

    for metal in fit_metals:
        x = ads_energies[metal]["CO"]
        y = ads_energies[metal]["O"]

        ax.scatter(x, y, color="black", s=35)
        ax.text(x + 0.03, y + 0.03, metal, color="white", fontsize=9)

    ax.set_xlabel("E_CO / eV")
    ax.set_ylabel("E_O / eV")
    ax.set_title(f"{mechanism['title']} Volcano Heatmap")
    fig.tight_layout()
    return fig


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    for mechanism_name, mechanism in mechanisms.items():
        E_CO_mesh, E_O_mesh, Z_plot = build_surface_grid(mechanism["surface_fn"])
        plot_volcano_surface(mechanism_name, E_CO_mesh, E_O_mesh, Z_plot)
        heatmap_fig = plot_volcano_heatmap(mechanism_name, E_CO_mesh, E_O_mesh, Z_plot)
        heatmap_fig.savefig(
            RESULTS_DIR / f"{mechanism_name}_volcano_heatmap.png",
            dpi=300,
            bbox_inches="tight",
        )

    plt.show()


if __name__ == "__main__":
    main()
