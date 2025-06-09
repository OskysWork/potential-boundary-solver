import numpy as np
import matplotlib.pyplot as plt

"""
Troubleshoot notes:

- Identified +-inf in i-1 u velocity profiles (start of range(1, Ny-1))
- From about count 19 there is inf error, values leading to are very large
"""

Nx = 200
Ny = 200
delta0 = 0.001
L = 0.05

x = np.linspace(0, L, Nx)
y = np.linspace(0, delta0, Ny)
dx = L / (Nx-1)         # i.e: dx = x[1] - x[0]
dy = delta0 / (Ny-1)    # i.e: dy = y[1] - y[0]

rho = 1
mu = 1.95*1e-5
dp_dx = 0

u = np.zeros((Nx, Ny))
v = np.zeros((Nx, Ny))
nut = np.zeros((Nx, Ny))
mut = np.zeros((Nx, Ny))

n1 = 0
n1s = True
n2 = 0
n2s = True
n3 = 0
n3s = True
n4 = 0
n4s = True
count = 0

for j in range(Ny):
    eta = y[j]/delta0
    u[0, j] = 1.5 * eta - 0.5 * eta**3 if eta<1 else 1


for i in range(1, Nx):
    delta = delta0 + 0.005 * x[i]     # Placeholder layer growth relation
    count+=1
    print(count)

    for j in range(Ny):
        eta = y[j] / delta
        nut[i-1, j] = 0.01 * (eta*(1-eta)) if eta<1 else 0
        mut[i-1, j] = rho * nut[i-1, j]

    for j in range(1, Ny-1):
        du_dy = (u[i-1, j+1] - u[i-1, j-1]) / (2*dy)

        d2u_dy2 = (u[i-1, j+1] - 2*u[i-1, j] + u[i-1, j-1]) / (dy)**2

        visc_term = (mu + mut[i-1, j])*d2u_dy2 + (
            (mut[i-1, j+1] - mut[i-1, j-1])/(2*dy)
            ) * du_dy
        
        du_dx = (
            (visc_term - dp_dx)/rho - v[i-1, j]*du_dy
            )/u[i-1, j] if abs(u[i-1, j]) > 1e-8 else 0
        
        u[i, j] = u[i-1, j] + dx * du_dx

    u[i, 0] = 0             # No slip condition
    u[i, -1] = 1

    v[i, -1] = 0
    for j in reversed(range(1, Ny)):
        du_dx = (u[i, j] - u[i-1, j]) / dx
        v[i, j-1] = v[i, j] - dy * du_dx


plt.plot(y, u[0, :])
plt.show()