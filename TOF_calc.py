import numpy as np
import csv
from model_data import (
    VIB_DATA,
    activation_energies,
    ads_energies,
    enthalpies,
    h,
    h_planck,
    k_B,
    metals,
    p_CO,
    p_CO2,
    p_H2,
    p_H2O,
    T,
)

# ----------------------
# Block 2: Equations For general equilibrium constant derivations
# ----------------------

# Entropy Calculation
def S(vibs):
  entropy = 0
  for i in vibs:
    entropy += k_B * np.log(1 - np.exp((-h * i) / (k_B * T)))
  return -entropy

def zpe(vibs):
  sum = 0
  for i in vibs:
    sum += h * i
  return 0.5 * sum

# Gibbs Free Energy Calculation
def dG(enthalpy, vibs_reacts, vibs_prods, E2=0):

  dG = 0
  s_react = 0
  s_prod = 0
  zpe_react = 0
  zpe_prod = 0

  for i in vibs_reacts:
    zpe_react += zpe(i)
    s_react += S(i)

  for j in vibs_prods:
    zpe_prod += zpe(j)
    s_prod += S(j)

  if enthalpy != 0:
    dG = enthalpy - (s_prod-s_react)*T

  else:
    G1 = E2 + zpe_prod - s_prod*T
    G2 = zpe_react - s_react*T
    dG = G1-G2

  return dG

# Equilibrium Constant Calculation
def K(delta_G):
  return np.exp(-delta_G / (k_B * T))

# ----------------------
# Block 3: Mechanism Specific Calculations
# ----------------------
# Redox Mechanism
#
# Steps
# (1) CO(g) + * <=> CO*
# (2) H2O(g) + * <=> H2O*
# (3) CO* + H2O* + * <=> OH* + H* + CO*
# (4) OH* + H* + CO* + * <=> O* + 2H* + CO*
# (5) CO* + O* + 2H* => CO2(g) + 2* + 2H*   [RDS]
# (6) CO2(g) + 2H* <=> CO2(g) + 2* + H2(g)

def wgsr_redox_tof(metal):
    # Step 1: CO(g) + * <=> CO*
    K1 = K(dG(0,
        [VIB_DATA["CO_NIST"]],
        [VIB_DATA["CO"]],
        ads_energies[metal]["CO"]
    ))

    # Step 2: H2O(g) + * <=> H2O*
    K2 = K(dG(
        0,
        [VIB_DATA["H2O_NIST"]],
        [VIB_DATA["H2O"]],
        ads_energies[metal]["H2O"]
    ))

    # Step 3: CO* + H2O* + * <=> OH* + H* + CO*
    K3 = K(dG(
        enthalpies["ΔH1: CO* + H2O* -> CO* + OH* + H*"][metal],
        [VIB_DATA["H2O"], VIB_DATA["CO"]],
        [VIB_DATA["OH"], VIB_DATA["H"], VIB_DATA["CO"]]
    ))

    # Step 4: OH* + * <=> O* + H*
    K4 = K(dG(
        enthalpies["ΔH2: CO* + OH* + H* -> CO* + O* + 2H*"][metal],
        [VIB_DATA["OH"], VIB_DATA["H"], VIB_DATA["CO"], VIB_DATA["H"]],
        [VIB_DATA["O"], VIB_DATA["H"], VIB_DATA["CO"], VIB_DATA["H"]]
    ))

    # Step 6: 2H* <=> H2(g) + 2*
    K6 = K(dG(
        enthalpies["ΔH13: CO2(g) + 2H* -> CO2(g) + H2(g)"][metal],
        [VIB_DATA["H"], VIB_DATA["H"]],
        [VIB_DATA["H2_NIST"]]
    ))

    # RDS: Step 5
    k5 = (k_B * T / h_planck) * np.exp(
        -activation_energies["Ea3: CO* + O* + 2H* -> CO2(g) + 2H*"][metal] / (k_B * T)
    )

    D = (
        1
        + K1 * p_CO
        + K2 * p_H2O
        + K2 * K3 * p_H2O * (K6 / p_H2)**0.5
        + (K2 * K3 * K4 * K6 * p_H2O) / p_H2
        + (p_H2 / K6)**0.5
    )

    TOF = k5 * K1 * K2 * K3 * K4 * K6 * (p_CO * p_H2O / p_H2) / D**2

    return np.log10(TOF)
# ----------------------
# Carboxyl Mechanism
#
# Steps:
# 1) H2O + * <=> H2O*
# 2) CO + * <=> CO*
# 3) H2O* + * <=> OH* + H*
# 4) CO* + OH* --> COOH + * [RATE DETERMINING STEP]
# 5) COOH <=> CO2 + H*
# 6) 2H* <=> H2 + 2*
#
# TOF = k_4[CO*][OH*]
#
# \theta_{star} = (1 + K1[H2O] + k2[CO] + k1K3[H2O]\sqrt{\frac{K6}{[H2]}}
# + \frac{[CO2]}{K5}\sqrt{\frac{[H2]}{K6}} + \sqrt{\frac{[H2]}{K6}})^-1

def carb_TOF(metal):
  #Detla G calculations
  dG1_carb = dG(0, [VIB_DATA["H2O_NIST"]], [VIB_DATA["H2O"]], ads_energies[metal]["H2O"])
  dG2_carb = dG(0, [VIB_DATA['CO_NIST']], [VIB_DATA["CO"]], ads_energies[metal]["CO"])
  dG3_carb = dG(enthalpies['ΔH1: CO* + H2O* -> CO* + OH* + H*'][metal], [VIB_DATA['H2O']],
                [VIB_DATA["OH"], VIB_DATA["H"]])
  dG4_carb = dG(enthalpies['ΔH4: CO* + OH* + H* -> COOH* + H*'][metal], [VIB_DATA["CO"], VIB_DATA["OH"]],
                [VIB_DATA["COOH"]])
  dG5_carb = dG(enthalpies['ΔH5: COOH* + H* -> CO2(g) + 2H*'][metal], [VIB_DATA["COOH"]], [VIB_DATA['CO2_NIST'], VIB_DATA['H']])
  dG6_carb = dG(enthalpies['ΔH13: CO2(g) + 2H* -> CO2(g) + H2(g)'][metal], [VIB_DATA['H']], [VIB_DATA['H2_NIST']])

  #rds k calculation
  k4_carb = ((k_B * T)/ h_planck) * np.exp(-activation_energies['Ea4: CO* + OH* + H* -> COOH* + H*'][metal] / (k_B * T))

  #K calculations
  K1_carb = K(dG1_carb)
  K2_carb = K(dG2_carb)
  K3_carb = K(dG3_carb)
  K5_carb = K(dG5_carb)
  K6_carb = K(dG6_carb)
  #Overal TOF calculations
  theta_star_carb = (1 + K1_carb * p_H2O + K2_carb * p_CO + K1_carb * K3_carb * p_H2O * np.sqrt(K6_carb / p_H2) +
                   (p_CO2 / K5_carb) * np.sqrt(p_H2 / K6_carb) + np.sqrt(p_H2 / K6_carb)) ** (-1)
  TOF_carb = k4_carb * K1_carb * K2_carb * K3_carb * p_CO * p_H2O * np.sqrt(K6_carb / p_H2) * theta_star_carb ** (2)
  return np.log10(TOF_carb)


# Steps:
# 1) CO + * <=> CO*
# 2) H2O + * <=> H2O*
# 3) H2O* <=> OH* + H*
# 4) OH* + H* <=> O* + 2H*
# 5) CO* + H* <=> CHO* + *
# 6) CHO* + O* <=> HCOO**
# 7) CHO* + OH* <=> HCOOH**
# 8) HCOOH** <=> HCOO** + H*
# 9) HCOO** --> CO2 + H* [RATE DETERMINING STEP]
# 10) 2H* <=> H2 + 2*
#
# TOF = k_9[HCOO**]

def TOF_formate(metal):

  dG1f = dG(0, [VIB_DATA['CO_NIST']], [VIB_DATA['CO']], ads_energies[metal]["CO"])
  dG2f = dG(0, [VIB_DATA['H2O_NIST']], [VIB_DATA['H2O']], ads_energies[metal]["H2O"])
  dG3f = dG(enthalpies["ΔH1: CO* + H2O* -> CO* + OH* + H*"][metal], [VIB_DATA['H2O']], [VIB_DATA['OH'], VIB_DATA['H']])
  dG4f = dG(enthalpies["ΔH2: CO* + OH* + H* -> CO* + O* + 2H*"][metal], [VIB_DATA['OH']], [VIB_DATA['O'], VIB_DATA['H']])
  dG5f = dG(enthalpies["ΔH6: CO* + O* + 2H* -> CHO* + O* + H*"][metal], [VIB_DATA['CO'], VIB_DATA['H']], [VIB_DATA['CHO']])
  dG6f = dG(enthalpies["ΔH7: CHO* + O* + H* -> HCOO* + H*"][metal], [VIB_DATA['CHO'], VIB_DATA['O']], [VIB_DATA['HCOO']])
  dG7f = dG(enthalpies["ΔH9: CHO* + OH* -> HCOOH*"][metal], [VIB_DATA['CHO'], VIB_DATA['OH']], [VIB_DATA['HCOOH']])
  dG8f = dG(enthalpies["ΔH10: HCOOH* -> HCOO* + H*"][metal], [VIB_DATA['HCOOH']], [VIB_DATA['HCOO'], VIB_DATA['H']])
  dG10f = dG(enthalpies["ΔH13: CO2(g) + 2H* -> CO2(g) + H2(g)"][metal], [VIB_DATA['H'], VIB_DATA['H']], [VIB_DATA['H2_NIST']])

  K1_f = K(dG1f)
  K2_f = K(dG2f)
  K3_f = K(dG3f)
  K4_f = K(dG4f)
  K5_f = K(dG5f)
  K6_f = K(dG6f)
  K7_f = K(dG7f)
  K8_f = K(dG8f)
  K10_f = K(dG10f)

  k9_f = ((k_B*T)/h_planck) * np.exp((-activation_energies['Ea8: HCOO* + H* -> CO2(g) + 2H*'][metal])/(k_B*T))

  A = 1 + K1_f * p_CO + K2_f * p_H2O + np.sqrt(p_H2/K10_f) + (K2_f * K3_f * p_H2O)/np.sqrt(p_H2/K10_f) + (K2_f * K3_f * K4_f * p_H2O)/(p_H2/K10_f) + K1_f * K5_f * p_CO * np.sqrt(p_H2/K10_f)
  B = K1_f * K2_f * K3_f * K4_f *K5_f * K6_f * p_CO * p_H2O * np.sqrt(K10_f/p_H2)
  C = K1_f * K2_f * K3_f * K5_f *K7_f * p_CO * p_H2O
  D = 2 * (B + C)

  theta = 2 / (A + np.sqrt(A**2 + 4*D))

  TOF = k9_f * B * theta**2
  #TOF = k9_f * K1_f * K2_f * K3_f * K4_f *K5_f * K6_f * p_CO * p_H2O * np.sqrt(K10_f/p_H2)

  return np.log10(TOF)

def calculate_tof_rows():
    rows = []

    for metal in metals:
        values = {
            "redox": wgsr_redox_tof(metal),
            "carboxyl": carb_TOF(metal),
            "formate": TOF_formate(metal),
        }

        for mechanism, log10_tof in values.items():
            rows.append(
                {
                    "metal": metal,
                    "mechanism": mechanism,
                    "dE_O": ads_energies[metal]["O"],
                    "dE_CO": ads_energies[metal]["CO"],
                    "log10_tof": log10_tof,
                    "tof": 10 ** log10_tof,
                }
            )

    return rows


def write_tof_results(output_csv="tof_results.csv"):
    fieldnames = ["metal", "mechanism", "dE_O", "dE_CO", "log10_tof", "tof"]
    rows = calculate_tof_rows()

    with open(output_csv, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote TOF results to {output_csv}")


def main():
    write_tof_results()


if __name__ == "__main__":
    main()
