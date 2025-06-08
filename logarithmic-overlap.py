import numpy as np
import matplotlib.pyplot as plt

"""
# Equation writing:

import sympy as s
s.init_printing()
k, yp, lam, up, Cp, A = s.symbols('κ y⁺ λ u⁺ C⁺ A')
e = s.E
dudy_eq = (1/(k*yp))*(1 - e**(-yp/A))**2
dudy = s.Derivative(up, yp)
eq = s.Eq(dudy, dudy_eq)
cp_log = (s.log(A) - 1.5)*(1/k)
eq2 = s.Eq(Cp, cp_log)
s.pprint(eq)
s.pprint(eq2)
"""

# Universal variables:
kappa = 0.41

"""
Gersten and Herwig (1992) p.406

Using ansatz:

du⁺   λ       1         1  
─── = ─ + ────────── + ────
dy⁺   κ   κ⋅(2 - y⁺)   y⁺⋅κ

Gives:

_   λ + log(2)
C = ──────────
        κ  
"""

# Variables
lam = np.linspace(-1, 1, 10)
C_bar_GH_exp = 2.1

# Solutions
C_bar_GH = (lam + np.log(2))/kappa
lam_exp = kappa*C_bar_GH_exp - np.log(2)


# Plotting
plt.figure()
plt.title(r"$\bar{C}$ vs $\lambda$")
plt.plot(lam, C_bar_GH)
plt.plot(lam_exp, C_bar_GH_exp, 'ro', label=fr'Experimental $\bar{{C}}$ ≈ {C_bar_GH_exp:.1f} at λ = {lam_exp:.4f}')
plt.xlabel('λ')
plt.ylabel(r"$\bar{C}$")
plt.legend()
plt.grid(True)
plt.show(block=False)

"""
van Driest

Using ansatz:
                  2
       ⎛     -y⁺ ⎞ 
       ⎜     ────⎟ 
       ⎜      A  ⎟ 
du⁺    ⎝1 - ℯ    ⎠ 
─── = ────────────
dy⁺       y⁺⋅κ    

Gives:
     log(A) - 1.5
C⁺ = ────────────
          κ      
"""

# Variables
A = np.linspace(1, 17, 100)
C_plus_D_exp = 5

# Solutions
C_plus_D = np.log(A)/kappa - 1.5
A_exp = np.exp(kappa*(C_plus_D_exp + 1.5))


# Plotting
plt.figure()
plt.title(r"$C^+$ vs A")
plt.plot(A, C_plus_D)
plt.plot(A_exp, C_plus_D_exp, 'ro', label=fr'Experimental $C^+$ ≈ {C_plus_D_exp:.2f} at A = {A_exp:.4f}')
plt.xlabel('A')
plt.ylabel('C+')
plt.legend()
plt.grid(True)
plt.show()





