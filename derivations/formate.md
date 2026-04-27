# Formate Mechanism for WGS

## Mechanism

1) $\mathrm{CO} + {*} \rightleftharpoons \mathrm{CO}^{*}$
2) $\mathrm{H_2O} + {*} \rightleftharpoons \mathrm{H_2O}^{*}$
3) $\mathrm{H_2O}^{*} + {*} \rightleftharpoons \mathrm{OH}^{*} + \mathrm{H}^{*}$
4) $\mathrm{OH}^{*} + {*} \rightleftharpoons \mathrm{O}^{*} + \mathrm{H}^{*}$
5) $\mathrm{CO}^{*} + \mathrm{H}^{*} \rightleftharpoons \mathrm{CHO}^{*} + {*}$
6) $\mathrm{CHO}^{*} + \mathrm{O}^{*} \rightleftharpoons \mathrm{HCOO}^{**}$
7) $\mathrm{CHO}^{*} + \mathrm{OH}^{*} \rightleftharpoons \mathrm{HCOOH}^{**}$
8) $\mathrm{HCOOH}^{**} \rightleftharpoons \mathrm{HCOO}^{**} + \mathrm{H}^{*}$
9) $\mathrm{HCOO}^{**} \to \mathrm{CO_2} + \mathrm{H}^{*}$ (RDS)
10) $2\mathrm{H}^{*} \rightleftharpoons \mathrm{H_2} + 2{*}$

Steps 1–8 and 10 are quasi-equilibrated; step 9 is rate-determining.

## Rate

$$\mathrm{TOF} = k_9 [\mathrm{HCOO}^{**}]$$

## Equilibrium Relations

From steps 1, 2, and 10:

$$[\mathrm{CO}^{*}] = K_1 P_{\mathrm{CO}}\, \theta$$

$$[\mathrm{H_2O}^{*}] = K_2 P_{\mathrm{H_2O}}\, \theta$$

$$K_{10} = \frac{P_{\mathrm{H_2}} \theta^2}{[\mathrm{H}^{*}]^2} \;\Rightarrow\; [\mathrm{H}^{*}] = \sqrt{\frac{P_{\mathrm{H_2}}}{K_{10}}}\, \theta$$

Define $\alpha = \sqrt{P_{\mathrm{H_2}} / K_{10}}$, so $[\mathrm{H}^{*}] = \alpha\theta$.

From step 3:

$$K_3 = \frac{[\mathrm{OH}^{*}][\mathrm{H}^{*}]}{[\mathrm{H_2O}^{*}]\,\theta}
\;\Rightarrow\;
[\mathrm{OH}^{*}] = \frac{K_3 [\mathrm{H_2O}^{*}]\, \theta}{[\mathrm{H}^{*}]}
= \frac{K_2 K_3 P_{\mathrm{H_2O}}\, \theta}{\alpha}$$

From step 4 (corrected stoichiometry: $\mathrm{OH}^{*} + {*} \rightleftharpoons \mathrm{O}^{*} + \mathrm{H}^{*}$):

$$K_4 = \frac{[\mathrm{O}^{*}][\mathrm{H}^{*}]}{[\mathrm{OH}^{*}]\,\theta}
\;\Rightarrow\;
[\mathrm{O}^{*}] = \frac{K_4 [\mathrm{OH}^{*}]\, \theta}{[\mathrm{H}^{*}]}
= \frac{K_2 K_3 K_4 P_{\mathrm{H_2O}}}{\alpha^2}\, \theta$$

Note: $[\mathrm{O}^{*}]$ is linear in $\theta$, not independent of it. This is the key correction.

From step 5:

$$K_5 = \frac{[\mathrm{CHO}^{*}]\, \theta}{[\mathrm{CO}^{*}][\mathrm{H}^{*}]}
\;\Rightarrow\;
[\mathrm{CHO}^{*}] = \frac{K_5 [\mathrm{CO}^{*}][\mathrm{H}^{*}]}{\theta}
= K_1 K_5 P_{\mathrm{CO}}\, \alpha\, \theta$$

## Formate via Route 1 (CHO + O)

$$K_6 = \frac{[\mathrm{HCOO}^{**}]}{[\mathrm{CHO}^{*}][\mathrm{O}^{*}]}
\;\Rightarrow\;
[\mathrm{HCOO}^{**}] = K_6 [\mathrm{CHO}^{*}][\mathrm{O}^{*}]$$

$$[\mathrm{HCOO}^{**}] = K_1 K_2 K_3 K_4 K_5 K_6\, P_{\mathrm{CO}} P_{\mathrm{H_2O}} \cdot \frac{1}{\alpha} \cdot \theta^2
= K_1 K_2 K_3 K_4 K_5 K_6\, P_{\mathrm{CO}} P_{\mathrm{H_2O}} \sqrt{\frac{K_{10}}{P_{\mathrm{H_2}}}}\, \theta^2$$

Define

$$B \equiv K_1 K_2 K_3 K_4 K_5 K_6\, P_{\mathrm{CO}} P_{\mathrm{H_2O}} \sqrt{\frac{K_{10}}{P_{\mathrm{H_2}}}}, \qquad [\mathrm{HCOO}^{**}] = B\, \theta^2$$

## Formate via Route 2 (CHO + OH, then deprotonation)

$$K_7 = \frac{[\mathrm{HCOOH}^{**}]}{[\mathrm{CHO}^{*}][\mathrm{OH}^{*}]}
\;\Rightarrow\;
[\mathrm{HCOOH}^{**}] = K_7 \cdot (K_1 K_5 P_{\mathrm{CO}}\, \alpha\, \theta) \cdot \left(\frac{K_2 K_3 P_{\mathrm{H_2O}}}{\alpha}\, \theta\right)
= K_1 K_2 K_3 K_5 K_7\, P_{\mathrm{CO}} P_{\mathrm{H_2O}}\, \theta^2$$

Define $C \equiv K_1 K_2 K_3 K_5 K_7\, P_{\mathrm{CO}} P_{\mathrm{H_2O}}$, so $[\mathrm{HCOOH}^{**}] = C\theta^2$.

From step 8:

$$K_8 = \frac{[\mathrm{HCOO}^{**}][\mathrm{H}^{*}]}{[\mathrm{HCOOH}^{**}]\,\theta}
\;\Rightarrow\;
[\mathrm{HCOO}^{**}] = \frac{K_8 [\mathrm{HCOOH}^{**}]\, \theta}{[\mathrm{H}^{*}]}
= \frac{K_8 \cdot C\theta^2 \cdot \theta}{\alpha\,\theta}
= \frac{K_8 C}{\alpha}\, \theta^2$$

Expanding $C$ and $\alpha$:

$$[\mathrm{HCOO}^{**}] = K_1 K_2 K_3 K_5 K_7 K_8\, P_{\mathrm{CO}} P_{\mathrm{H_2O}} \sqrt{\frac{K_{10}}{P_{\mathrm{H_2}}}}\, \theta^2$$

## Thermodynamic Consistency

Both routes must give identical $[\mathrm{HCOO}^{**}]$ at equilibrium. Equating the two expressions for $B\theta^2$:

$$\boxed{K_4 K_6 = K_7 K_8}$$

This is the cycle condition for the closed loop $\mathrm{CHO}^{*} + \mathrm{OH}^{*} \to \mathrm{HCOO}^{**} + \mathrm{H}^{*}$ traversed two ways.

## Final Rate Expression

$$\mathrm{TOF} = k_9\, B\, \theta^2
= k_9 K_1 K_2 K_3 K_4 K_5 K_6\, P_{\mathrm{CO}} P_{\mathrm{H_2O}} \sqrt{\frac{K_{10}}{P_{\mathrm{H_2}}}}\, \theta^2$$

The $\theta^2$ factor reflects the dual-site nature of $\mathrm{HCOO}^{**}$ and must be retained — it is not generally close to unity at WGS conditions.

## Site Balance

Counting each surface species by its site occupancy (1 for monodentate, 2 for bidentate $\mathrm{HCOO}^{**}$ and $\mathrm{HCOOH}^{**}$):

$$1 = \theta + [\mathrm{CO}^{*}] + [\mathrm{H_2O}^{*}] + [\mathrm{H}^{*}] + [\mathrm{OH}^{*}] + [\mathrm{O}^{*}] + [\mathrm{CHO}^{*}] + 2[\mathrm{HCOO}^{**}] + 2[\mathrm{HCOOH}^{**}]$$

Substituting:

$$\begin{aligned}
1 &= \theta + K_1 P_{\mathrm{CO}}\, \theta + K_2 P_{\mathrm{H_2O}}\, \theta + \alpha \theta + \frac{K_2 K_3 P_{\mathrm{H_2O}}}{\alpha}\theta + \frac{K_2 K_3 K_4 P_{\mathrm{H_2O}}}{\alpha^2}\theta \\
  &\quad + K_1 K_5 P_{\mathrm{CO}}\, \alpha\, \theta + 2B \theta^2 + 2C \theta^2
\end{aligned}$$

Collecting:

$$1 = A\, \theta + D\, \theta^2$$

where

$$A = 1 + K_1 P_{\mathrm{CO}} + K_2 P_{\mathrm{H_2O}} + \alpha + \frac{K_2 K_3 P_{\mathrm{H_2O}}}{\alpha} + \frac{K_2 K_3 K_4 P_{\mathrm{H_2O}}}{\alpha^2} + K_1 K_5 P_{\mathrm{CO}}\, \alpha$$

$$D = 2(B + C)$$

The site balance is a quadratic in $\theta$:

$$D\, \theta^2 + A\, \theta - 1 = 0$$

$$\theta = \frac{-A + \sqrt{A^2 + 4D}}{2D} = \frac{2}{A + \sqrt{A^2 + 4D}}$$

(positive root chosen since $0 \leq \theta \leq 1$; the second form is numerically stable when $D \ll A^2$).

## Summary of Corrections

1. Step 4 rebalanced: $\mathrm{OH}^{*} + {*} \rightleftharpoons \mathrm{O}^{*} + \mathrm{H}^{*}$ (the original had two $\mathrm{H}^{*}$ on the product side, which is unbalanced).
2. $[\mathrm{O}^{*}]$ now scales linearly in $\theta$, which propagates to give $[\mathrm{HCOO}^{**}] \propto \theta^2$ as expected for a bidentate species formed from two singly-bound surface intermediates.
3. TOF carries an explicit $\theta^2$ factor; the empty-site coverage cannot be set to unity a priori.
4. Site balance is quadratic, not the previous form $a\theta^2 + (b-1)\theta + c = 0$ that came from the spurious $1/\theta$ in $[\mathrm{O}^{*}]$.
5. Consistency condition $K_4 K_6 = K_7 K_8$ unchanged — it was correct in the original.
