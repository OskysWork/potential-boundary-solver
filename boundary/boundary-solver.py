import numpy as np
import matplotlib.pyplot as plt

"""
Troubleshoot notes:

- Identified +-inf in i-1 u velocity profiles (start of range(1, Ny-1))
- From about count 19 there is inf error, values leading to are very large
"""

Nx = 160000
Ny = 80
delta0 = 0.05
L = 0.4

x = np.linspace(0, L, Nx)
y = delta0 * (np.tanh(2 * np.linspace(0, 1, Ny)) / np.tanh(2))
#y = np.linspace(0, delta0, Ny)
dx = x[1] - x[0]
print(dx)
dy = y[1] - y[0]
print(dy**2)

rho = 1
mu = 1.95*1e-5
dp_dx = 0

u = np.zeros((Nx, Ny))
v = np.zeros((Nx, Ny))
nut = np.zeros((Nx, Ny))
mut = np.zeros((Nx, Ny))

count = 0

for j in range(Ny):
    eta = y[j]/delta0
    u[0, j] = 1.5 * eta - 0.5 * eta**3 if eta<1 else 1


for i in range(1, Nx):
    delta = delta0 + 0.005 * x[i]     # Placeholder layer growth relation
    count+=1
    #print(count)

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
plt.show(block=False)
plt.plot(y, u[10, :])
plt.show(block=False)
plt.plot(y, u[100, :])
plt.show(block=False)
plt.plot(y, u[1000, :])
plt.show(block=False)
plt.plot(y, u[10000, :], label='1')
plt.show(block=False)
plt.plot(y, u[20000, :], label='2')
plt.show(block=False)
plt.plot(y, u[30000, :], label='3')
plt.xlabel('y position')
plt.ylabel('Velocity distribution')
plt.legend()
plt.show()