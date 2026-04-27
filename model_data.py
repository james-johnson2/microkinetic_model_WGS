from __future__ import annotations

import csv
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"


def _read_constants(filename: str) -> dict[str, float]:
    with (DATA_DIR / filename).open(newline="", encoding="utf-8") as f:
        return {row["name"]: float(row["value"]) for row in csv.DictReader(f)}


def _read_metals(filename: str) -> list[str]:
    with (DATA_DIR / filename).open(newline="", encoding="utf-8") as f:
        return [row["metal"] for row in csv.DictReader(f)]


def _read_vibrations(filename: str) -> dict[str, list[float]]:
    with (DATA_DIR / filename).open(newline="", encoding="utf-8") as f:
        return {
            row["species"]: [float(value) for value in row["frequencies_cm-1"].split(";")]
            for row in csv.DictReader(f)
        }


def _read_reaction_table(filename: str) -> dict[str, dict[str, float]]:
    with (DATA_DIR / filename).open(newline="", encoding="utf-8") as f:
        rows = csv.DictReader(f)
        return {
            row["reaction"]: {
                metal: float(value)
                for metal, value in row.items()
                if metal != "reaction" and value != ""
            }
            for row in rows
        }


def _read_species_by_metal(filename: str) -> dict[str, dict[str, float]]:
    with (DATA_DIR / filename).open(newline="", encoding="utf-8") as f:
        rows = csv.DictReader(f)
        return {
            row["metal"]: {
                species: float(value)
                for species, value in row.items()
                if species != "metal" and value != ""
            }
            for row in rows
        }


def _read_value_by_metal(filename: str, value_column: str) -> dict[str, float]:
    with (DATA_DIR / filename).open(newline="", encoding="utf-8") as f:
        return {row["metal"]: float(row[value_column]) for row in csv.DictReader(f)}


constants = _read_constants("constants.csv")

k_B = constants["k_B"]
T = constants["T"]
h = constants["h"]
h_planck = constants["h_planck"]
p_CO = constants["p_CO"]
p_H2O = constants["p_H2O"]
p_H2 = constants["p_H2"]
p_CO2 = constants["p_CO2"]

metals = _read_metals("metals.csv")
VIB_DATA = _read_vibrations("vibrational_frequencies.csv")
enthalpies = _read_reaction_table("enthalpies.csv")
activation_energies = _read_reaction_table("activation_energies.csv")
ads_energies = _read_species_by_metal("adsorption_energies.csv")
exp_tof = _read_value_by_metal("experimental_tof.csv", "tof")
