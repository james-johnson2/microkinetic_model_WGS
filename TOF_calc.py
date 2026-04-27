import numpy as np
# ----------------------
# Block 1: All Data
# ----------------------

k_B = 8.617333262 * (10 ** (-5))
T = 573.15 #Kelvin
h = 1.23984 * (10 ** (-4))
h_planck = 4.135667696 * (10 **(-15))
p_CO = 0.243
p_H2O = 0.314
p_H2 = 0.001
p_CO2 = 0.001
metals = ["Co", "Ni", "Cu", "Rh", "Pd", "Ag", "Ir", "Pt", "Au"]

# ----------------------
# Vibrational Data
# Source for Vibrational Frequencies: https://pubs.acs.org/doi/full/10.1021/jp7099702
# ASSUMPTION: All for Pt metal; however, the vibrational frequencies are assumed to be the same
# for all the metals. This is common in literature.
VIB_DATA = {
    "H":    [1091, 556, 548],
    "O":    [452, 371, 369],
    "OH":   [3744, 894, 502],
    "H2O":  [3794, 3690, 1526, 533, 452],
    "CO":   [1867, 332, 298, 298],
    "CO2":  [2438, 1348, 620, 617],
    "COOH": [3541, 1702, 1208, 1126, 659, 615, 505, 281, 267],   # top H down
    "HCOO": [2982, 1557, 1307, 1299, 954, 739, 348, 322, 299],
    "CHO":  [2875, 1784, 1141, 791, 530, 282],
    "HCOOH": [3541, 1702, 1208, 1126, 659, 615, 505, 281, 267],
    "H2_NIST": [4161], # NIST for the gas species
    "CO2_NIST": [1333, 2349, 667],
    "H2O_NIST": [3657, 1595, 3756],
    "CO_NIST": [2143]
}
# ----------------------

# ----------------------
# Enthalpy Data
# Source for Enthalpy Data: https://pubs.acs.org/doi/10.1021/jp2034467
enthalpies = {
    "ΔH1: CO* + H2O* -> CO* + OH* + H*": {
        "Co": -0.77, "Ni": -0.41, "Cu": 0.21, "Rh": 0.26, "Pd": 1.65,
        "Ag": 1.93, "Ir": 1.18, "Pt": 1.63, "Au": 1.43
    },
    "ΔH2: CO* + OH* + H* -> CO* + O* + 2H*": {
        "Co": -0.07, "Ni": 0.03, "Cu": 0.75, "Rh": -0.09, "Pd": 0.36,
        "Ag": 1.73, "Ir": -0.18, "Pt": 0.14, "Au": 1.49
    },
    "ΔH3: CO* + O* + 2H* -> CO2(g) + 2H*": {
        "Co": 0.89, "Ni": 0.73, "Cu": -0.90, "Rh": 0.57, "Pd": -0.51,
        "Ag": -2.69, "Ir": 0.19, "Pt": -0.59, "Au": -2.93
    },
    "ΔH4: CO* + OH* + H* -> COOH* + H*": {
        "Co": 1.35, "Ni": 1.22, "Cu": 0.44, "Rh": 0.84, "Pd": 0.37,
        "Ag": -0.40, "Ir": 0.17, "Pt": -0.25, "Au": -0.95
    },
    "ΔH5: COOH* + H* -> CO2(g) + 2H*": {
        "Co": -0.56, "Ni": -0.52, "Cu": -0.59, "Rh": -0.56, "Pd": -0.52,
        "Ag": -0.56, "Ir": -0.56, "Pt": -0.59, "Au": -0.49
    },
    "ΔH6: CO* + O* + 2H* -> CHO* + O* + H*": {
        "Co": 1.14, "Ni": 1.18, "Cu": 0.64, "Rh": 1.19, "Pd": 1.11,
        "Ag": 0.71, "Ir": 0.73, "Pt": 0.61, "Au": 0.52
    },
    "ΔH7: CHO* + O* + H* -> HCOO* + H*": {
        "Co": -0.32, "Ni": -0.39, "Cu": -2.96, "Rh": -0.66, "Pd": -1.22,
        "Ag": -3.20, "Ir": -0.89, "Pt": -1.52, "Au": -2.54
    },
    "ΔH8: HCOO* + H* -> CO2(g) + 2H*": {
        "Co": -0.06, "Ni": -0.06, "Cu": 0.18, "Rh": -0.10, "Pd": -0.35,
        "Ag": -0.26, "Ir": -0.07, "Pt": -0.44, "Au": -0.83
    },
    "ΔH9: CHO* + OH* -> HCOOH*": {
        "Co": 0.29, "Ni": 0.14, "Cu": -0.74, "Rh": -0.07, "Pd": -0.64,
        "Ag": -1.59, "Ir": -0.29, "Pt": -0.51, "Au": -1.72
    },
    "ΔH10: HCOOH* -> HCOO* + H*": {
        "Co": -0.88, "Ni": -0.80, "Cu": -0.47, "Rh": -0.79, "Pd": -0.41,
        "Ag": -0.20, "Ir": -0.85, "Pt": -0.42, "Au": -0.60
    },
    "ΔH11: HCOOH* -> COOH* + H*": {
        "Co": 0.84, "Ni": 0.86, "Cu": 1.22, "Rh": 0.84, "Pd": 0.99,
        "Ag": 1.80, "Ir": 0.94, "Pt": 1.06, "Au": 1.94
    },
    "ΔH12: COOH* + H* -> CO2(g) + 2H*": {
        "Co": -1.68, "Ni": -1.75, "Cu": -1.52, "Rh": -1.73, "Pd": -1.95,
        "Ag": -1.35, "Ir": -1.72, "Pt": -1.94, "Au": -1.71
    },
    "ΔH13: CO2(g) + 2H* -> CO2(g) + H2(g)": {
        "Co": 1.19, "Ni": 1.18, "Cu": 0.66, "Rh": 1.16, "Pd": 1.15,
        "Ag": 0.31, "Ir": 0.98, "Pt": 1.03, "Au": 0.45
    },
}
# ----------------------

# ----------------------
# Activation Energy Data
# Source for Activation Energy Data: https://pubs.acs.org/doi/10.1021/jp2034467
activation_energies = {
    "Ea1: CO* + H2O* -> CO* + OH* + H*": {
        "Co": 0.45, "Ni": 0.79, "Cu": 1.13, "Rh": 1.00, "Pd": 1.24,
        "Ag": 1.62, "Ir": 0.80, "Pt": 1.02, "Au": 2.03
    },
    "Ea2: CO* + OH* + H* -> CO* + O* + 2H*": {
        "Co": 0.98, "Ni": 1.09, "Cu": 1.74, "Rh": 0.81, "Pd": 1.35,
        "Ag": 2.31, "Ir": 0.76, "Pt": 0.95, "Au": 2.12
    },
    "Ea3: CO* + O* + 2H* -> CO2(g) + 2H*": {
        "Co": 1.45, "Ni": 1.44, "Cu": 0.62, "Rh": 1.34, "Pd": 0.84,
        "Ag": 0.90, "Ir": 1.11, "Pt": 0.74, "Au": 0.03
    },
    "Ea4: CO* + OH* + H* -> COOH* + H*": {
        "Co": 1.56, "Ni": 1.45, "Cu": 0.88, "Rh": 1.06, "Pd": 0.97,
        "Ag": 0.43, "Ir": 0.79, "Pt": 0.59, "Au": 0.30
    },
    "Ea5: COOH* + H* -> CO2(g) + 2H*": {
        "Co": 0.55, "Ni": 0.60, "Cu": 0.60, "Rh": 0.68, "Pd": 0.55,
        "Ag": 0.61, "Ir": 0.59, "Pt": 0.66, "Au": 0.58
    },
    "Ea6: CO* + O* + 2H* -> CHO* + O* + H*": {
        "Co": 1.40, "Ni": 1.45, "Cu": 0.88, "Rh": 1.43, "Pd": 1.40,
        "Ag": 0.39, "Ir": 1.12, "Pt": 0.93, "Au": 0.37
    },
    "Ea7: CHO* + O* + H* -> HCOO* + H*": {
        "Co": 0.97, "Ni": 0.94, "Cu": 0.19, "Rh": 0.74, "Pd": 0.66,
        "Ag": 0.12, "Ir": 0.91, "Pt": 0.34, "Au": 0.18
    },
    "Ea8: HCOO* + H* -> CO2(g) + 2H*": {
        "Co": 1.31, "Ni": 1.36, "Cu": 1.26, "Rh": 1.35, "Pd": 1.22,
        "Ag": 1.25, "Ir": 1.38, "Pt": 1.26, "Au": 1.24
    },
    "Ea9: CHO* + OH* -> HCOOH*": {
        "Co": 0.97, "Ni": 0.82, "Cu": 0.49, "Rh": 0.69, "Pd": 0.60,
        "Ag": 0.18, "Ir": 0.52, "Pt": 0.59, "Au": 0.27
    },
    "Ea10: HCOOH* -> HCOO* + H*": {
        "Co": 0.22, "Ni": 0.80, "Cu": 0.34, "Rh": 0.38, "Pd": 0.88,
        "Ag": 0.96, "Ir": 0.37, "Pt": 0.75, "Au": 1.14
    },
    "Ea11: HCOOH* -> COOH* + H*": {
        "Co": 1.07, "Ni": 1.43, "Cu": 1.82, "Rh": 1.41, "Pd": 1.63,
        "Ag": 2.39, "Ir": 1.27, "Pt": 1.47, "Au": 2.41
    },
    "Ea12: COOH* + H* -> CO2(g) + 2H*": {
        "Co": 0.15, "Ni": 0.17, "Cu": 0.21, "Rh": 0.20, "Pd": 0.26,
        "Ag": 0.24, "Ir": 0.18, "Pt": 0.15, "Au": 0.23
    },
    "Ea13: CO2(g) + 2H* -> CO2(g) + H2(g)": {
        "Co": 1.15, "Ni": 1.18, "Cu": 0.43, "Rh": 1.14, "Pd": 1.15,
        "Ag": 0.37, "Ir": 0.93, "Pt": 0.99, "Au": 0.38
    },
}

# Adsorption energy data
# Source for data: https://pubs.acs.org/doi/10.1021/cs400664z
ads_energies = {
    "Au":  {"HCOOH": -0.16, "HCOO": -1.78, "COOH": -1.35, "CO": -0.26, "OH": -1.56, "O": -2.47, "C": -3.68, "H": -2.04, "H2O": -0.5},
    "Ag":  {"HCOOH": -0.18, "HCOO": -2.31, "COOH": -1.19, "CO": -0.10, "OH": -2.28, "O": -3.16, "C": -3.15, "H": -2.07, "H2O": -0.5},
    "Cu":  {"HCOOH": -0.23, "HCOO": -2.76, "COOH": -1.52, "CO": -0.72, "OH": -2.69, "O": -4.14, "C": -4.30, "H": -2.38, "H2O": -0.5},
    "Pt":  {"HCOOH": -0.37, "HCOO": -2.35, "COOH": -2.40, "CO": -1.74, "OH": -2.11, "O": -3.73, "C": -6.57, "H": -2.70, "H2O": -0.5},
    "Pd":  {"HCOOH": -0.39, "HCOO": -2.34, "COOH": -2.19, "CO": -1.95, "OH": -2.22, "O": -3.73, "C": -6.46, "H": -2.83, "H2O": -0.5},
    "Ni":  {"HCOOH": -0.32, "HCOO": -2.80, "COOH": -2.25, "CO": -1.90, "OH": -2.98, "O": -4.94, "C": -6.38, "H": -2.81, "H2O": -0.5},
    "Ir":  {"HCOOH": -0.43, "HCOO": -2.94, "COOH": -2.63, "CO": -1.83, "OH": -2.60, "O": -4.69, "C": -6.88, "H": -2.73, "H2O": -0.5},
    "Rh":  {"HCOOH": -0.49, "HCOO": -2.94, "COOH": -2.58, "CO": -1.92, "OH": -2.72, "O": -4.74, "C": -7.05, "H": -2.82, "H2O": -0.5},
    "Co":  {"HCOOH": -0.33, "HCOO": -2.95, "COOH": -2.25, "CO": -1.83, "OH": -3.13, "O": -5.23, "C": -6.59, "H": -2.98, "H2O": -0.5},
    "Os":  {"HCOOH": -0.69, "HCOO": -3.62, "COOH": -3.05, "CO": -2.10, "OH": -3.05, "O": -5.47, "C": -7.37, "H": -2.81, "H2O": -0.5},
    "Ru":  {"HCOOH": -0.81, "HCOO": -3.50, "COOH": -3.00, "CO": -1.98, "OH": -3.22, "O": -5.29, "C": -6.99, "H": -2.90, "H2O": -0.5},
    "Re":  {"HCOOH": -0.64, "HCOO": -3.55, "COOH": -2.83, "CO": -1.92, "OH": -3.72, "O": -6.69, "C": -7.47, "H": -3.04, "H2O": -0.5}
}


exp_tof = {
    "Ru": 0.1929,
    "Rh": 0.0086,
    "Pd": 0.0135,
    "Os": 0.0615,
    "Ir": 0.0032,
    "Pt": 0.0635,
    "Co": 0.2472,
    "Ni": 0.1029,
    "Re": 0.3839,
    "Cu": 12.1900,
    "Au": 0.0356,
}
# ----------------------

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

def wgsr_redox_tof(metal, p_CO=0.243, p_H2O=0.314, p_H2=0.001, vib_data=VIB_DATA, enthalpy_db=enthalpies):
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

TOF_formate = []


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

import numpy as np
import matplotlib.pyplot as plt

# ----------------------
# Compute model TOFs
# ----------------------
redox_tofs = []
carboxyl_tofs = []
formate_tofs = []
CO_energies = []

for metal in metals:
    redox_tof = wgsr_redox_tof(metal)
    redox_tofs.append(redox_tof)

    carboxyl_tof = carb_TOF(metal)
    carboxyl_tofs.append(carboxyl_tof)

    formate_tof = TOF_formate(metal)
    formate_tofs.append(formate_tof)

    # FIXED: correct indexing
    CO_energies.append(ads_energies[metal]["CO"])

# ----------------------
# Extract overlapping metals
# ----------------------
common_metals = [m for m in metals if m in exp_tof]

# Align everything to common metals
exp_vals = [exp_tof[m] for m in common_metals]
redox_common = [redox_tofs[metals.index(m)] for m in common_metals]
carboxyl_common = [carboxyl_tofs[metals.index(m)] for m in common_metals]
formate_common = [formate_tofs[metals.index(m)] for m in common_metals]
CO_common = [ads_energies[m]["CO"] for m in common_metals]

# ----------------------
# Plot (CO adsorption energy on x-axis)
# ----------------------
plt.figure(figsize=(10,6))

plt.scatter(CO_common, exp_vals, label='Experimental')
plt.scatter(CO_common, redox_common, label='Redox')
plt.scatter(CO_common, carboxyl_common, label='Carboxyl')
plt.scatter(CO_common, formate_common, label='Formate')

plt.xlabel("CO Adsorption Energy (eV)")
plt.ylabel("TOF")
plt.title("TOF vs CO Adsorption Energy")
plt.legend()
plt.tight_layout()
plt.show()
