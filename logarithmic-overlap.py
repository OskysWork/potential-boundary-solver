import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.integrate import quad

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


"""
Spalding
### Scrappy code
"""

# Constants
log_const = 5
E = np.exp(-kappa*log_const)

def residual_function(u_plus, y_plus):
    """
	Defines residual func
    """
    # Example: Spalding's equation
    f_u = u_plus + (np.exp(kappa * u_plus) - 1 - kappa * u_plus - (kappa * u_plus)**2 / 2 - (kappa * u_plus)**3 / 6)*E
    return f_u - y_plus  # Residual = LHS - RHS


def dudy_plus_ansatz(u_plus):
    """
	Defines rearranged Spalding's gradient ansatz
    """
    numerator = kappa * np.exp(kappa * u_plus) - kappa - kappa**2 * u_plus - 0.5 * kappa**3 * u_plus**2
    derivative = 1 + numerator*E
    return 1 / derivative


def solve_u_plus(y_plus):
    sol = root_scalar(residual_function, args=(y_plus,), bracket=[0, 30], method='brentq')
    return sol.root


def integrand_1(y_plus):
    u_plus = solve_u_plus(y_plus)
    return dudy_plus_ansatz(u_plus)


def integrand_2(y_plus):
    u_plus = solve_u_plus(y_plus)
    return dudy_plus_ansatz(u_plus) - 1 / (kappa * y_plus)


y_plus_max = 10000


integral_1, _ = quad(integrand_1, 0, 1)
integral_2, _ = quad(integrand_2, 1, y_plus_max)


C_plus = integral_1 + integral_2
print(C_plus)





