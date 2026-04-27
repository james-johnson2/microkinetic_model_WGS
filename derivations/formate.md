# Formate Mechanism for WGS

## Mechanism

1) CO + * ⇌ CO*  
2) H2O + * ⇌ H2O*  
3) H2O* + * ⇌ OH* + H*  
4) OH* + * ⇌ O* + H*  
5) CO* + H* ⇌ CHO* + *  
6) CHO* + O* ⇌ HCOO**  
7) CHO* + OH* ⇌ HCOOH**  
8) HCOOH** ⇌ HCOO** + H*  
9) HCOO** → CO2 + H* (RDS)  
10) 2H* ⇌ H2 + 2*  

Steps 1–8 and 10 are quasi-equilibrated; step 9 is rate-determining.

---

## Rate

$$
\mathrm{TOF} = k_9 [\mathrm{HCOO**}]
$$

---

## Equilibrium Relations

From steps 1, 2, and 10:

$$
[\mathrm{CO*}] = K_1 P_{CO} \theta
$$

$$
[\mathrm{H2O*}] = K_2 P_{H2O} \theta
$$

$$
K_{10} = \frac{P_{H2} \theta^2}{[\mathrm{H*}]^2}
\Rightarrow
[\mathrm{H*}] = \sqrt{\frac{P_{H2}}{K_{10}}} \theta
$$

Define:

$$
\alpha = \sqrt{\frac{P_{H2}}{K_{10}}}, \quad [\mathrm{H*}] = \alpha \theta
$$

---

## From Step 3

$$
K_3 = \frac{[\mathrm{OH*}][\mathrm{H*}]}{[\mathrm{H2O*}] \theta}
$$

$$
[\mathrm{OH*}] = \frac{K_3 [\mathrm{H2O*}] \theta}{[\mathrm{H*}]}
= \frac{K_2 K_3 P_{H2O}}{\alpha} \theta
$$

---

## From Step 4

$$
K_4 = \frac{[\mathrm{O*}][\mathrm{H*}]}{[\mathrm{OH*}] \theta}
$$

$$
[\mathrm{O*}] = \frac{K_4 [\mathrm{OH*}] \theta}{[\mathrm{H*}]}
= \frac{K_2 K_3 K_4 P_{H2O}}{\alpha^2} \theta
$$

---

## From Step 5

$$
K_5 = \frac{[\mathrm{CHO*}] \theta}{[\mathrm{CO*}][\mathrm{H*}]}
$$

$$
[\mathrm{CHO*}] = \frac{K_5 [\mathrm{CO*}][\mathrm{H*}]}{\theta}
= K_1 K_5 P_{CO} \alpha \theta
$$

---

## Formate via Route 1 (CHO + O)

$$
K_6 = \frac{[\mathrm{HCOO**}]}{[\mathrm{CHO*}][\mathrm{O*}]}
$$

$$
[\mathrm{HCOO**}] =
K_1 K_2 K_3 K_4 K_5 K_6
P_{CO} P_{H2O}
\frac{1}{\alpha}
\theta^2
$$

Define:

$$
B = K_1 K_2 K_3 K_4 K_5 K_6 P_{CO} P_{H2O} \sqrt{\frac{K_{10}}{P_{H2}}}
$$

$$
[\mathrm{HCOO**}] = B \theta^2
$$

---

## Formate via Route 2 (CHO + OH)

$$
K_7 = \frac{[\mathrm{HCOOH**}]}{[\mathrm{CHO*}][\mathrm{OH*}]}
$$

$$
[\mathrm{HCOOH**}] =
K_1 K_2 K_3 K_5 K_7
P_{CO} P_{H2O}
\theta^2
$$

Define:

$$
C = K_1 K_2 K_3 K_5 K_7 P_{CO} P_{H2O}
$$

---

## From Step 8

$$
K_8 = \frac{[\mathrm{HCOO**}][\mathrm{H*}]}{[\mathrm{HCOOH**}] \theta}
$$

$$
[\mathrm{HCOO**}] =
\frac{K_8 C}{\alpha} \theta^2
$$

---

## Thermodynamic Consistency

$$
K_4 K_6 = K_7 K_8
$$

---

## Final Rate Expression

$$
\mathrm{TOF} =
k_9 K_1 K_2 K_3 K_4 K_5 K_6
P_{CO} P_{H2O}
\sqrt{\frac{K_{10}}{P_{H2}}}
\theta^2
$$

---

## Site Balance

$$
\begin{aligned}
1 =\;& \theta + [\mathrm{CO*}] + [\mathrm{H2O*}] + [\mathrm{H*}] \\
&+ [\mathrm{OH*}] + [\mathrm{O*}] + [\mathrm{CHO*}] \\
&+ 2[\mathrm{HCOO**}] + 2[\mathrm{HCOOH**}]
\end{aligned}
$$

Substitute:

$$
1 = A \theta + D \theta^2
$$

where:

$$
\begin{aligned}
A = 1 + K_1 P_{CO} + K_2 P_{H2O} + \alpha + \frac{K_2 K_3 P_{H2O}}{\alpha} + \frac{K_2 K_3 K_4 P_{H2O}}{\alpha^2} + K_1 K_5 P_{CO} \alpha
\end{aligned}
$$

$$
D = 2(B + C)
$$

---

## Solve for θ

$$
D \theta^2 + A \theta - 1 = 0
$$

$$
\begin{aligned}
\theta = \frac{-A + \sqrt{A^2 + 4D}}{2D} = \frac{2}{A + \sqrt{A^2 + 4D}}
\end{aligned}
$$
