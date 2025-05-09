import os
import sys
import matplotlib.pyplot as plt
import numpy as np

current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)

from src.utilities import read_aerofoil
from src.panels import discretise_panels, solve_panels, compute_cp

x, y = read_aerofoil("../data/NACA6409.dat")
#x, y = read_aerofoil("../data/1%-Thickness.dat")
#x, y = read_aerofoil("../data/NACA0006.dat")
alpha = np.radians(5)

panels = discretise_panels(x, y)
panels = solve_panels(panels, U_free=1.0, alpha=alpha)

gamma_0 = [panel.gamma_0 for panel in panels]
gamma_1 = [panel.gamma_1 for panel in panels]

xp, Cp = compute_cp(panels, U_free=1.0, alpha=alpha)

#print(gamma_0)
#print(gamma_1)
N = np.linspace(1, len(x)-1, len(x)-1)

plt.plot(N, gamma_0, label='Gamma 0')
plt.plot(N, gamma_1, label='Gamma 1')
plt.legend()

plt.show()

plt.figure()
plt.plot(xp, Cp, '-o')
plt.gca().invert_yaxis()
plt.xlabel('x/c')
plt.ylabel('$C_p$')
plt.title('Pressure Coefficient')
plt.show()