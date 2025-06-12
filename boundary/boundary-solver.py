import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


Nx = 160000
Ny = 80
delta0 = 0.05
L = 0.02 

x = np.linspace(0, L, Nx)
#y = np.linspace(0, delta0, Ny)
y = np.linspace(0, 1, 80)**1.3
y = y/y[-1] * delta0
dx = x[1] - x[0]
#dy = y[1] - y[0]
dy = [y[1] - y[0]] + [(y[i+1] - y[i-1]) / 2 for i in range(1, len(y)-1)] + [y[-1] - y[-2]]

print(dx)
print(dy[0]**2)
print(dy[1]**2)

rho = 1
mu = 1.95*1e-5
nu = mu / rho
dp_dx = 0.01#0.01*1e5

u = np.zeros((Nx, Ny))
v = np.zeros((Nx, Ny))
nut = np.zeros((Nx, Ny))
mut = np.zeros((Nx, Ny))
tau_v = np.zeros((Nx, Ny))
tau_t = np.zeros((Nx, Ny))
dudy = np.zeros((Nx, Ny))
kappa = 0.41

count = 0

for j in range(Ny):
    eta = y[j]/delta0
    u[0, j] = 1.5 * eta - 0.5 * eta**3 if eta<1 else 1

"""
for j in range(Ny-1):
    dudy[0, j] = (u[0, j+1] - u[0, j-1]) / (2*dy[j])

dudy[0, 0] = (u[0, 1] - u[0, 0]) / dy[0]
dudy[0, -1] = (u[0, -1] - u[0, -2]) / dy[-1]
"""

for i in tqdm(range(1, Nx)):
    #delta = delta0 + 0.005 * x[i]     # Placeholder layer growth relation
    delta = delta0 + np.sqrt(x[i]*mu*rho)
    count+=1
    #print(count)


    eta = y / delta
    nut[i-1, :] = np.where(eta < 1, 0.01 * (eta*(1-eta)), 0)
    #nut[i-1, :] = ((kappa*y)**2)*abs(dudy[i, :])
    mut[i-1, :] = rho * nut[i-1, :]

    for j in range(1, Ny-1):
        du_dy = (u[i-1, j+1] - u[i-1, j-1]) / (2*dy[j])
        dudy[i, j] = du_dy

        d2u_dy2 = (u[i-1, j+1] - 2*u[i-1, j] + u[i-1, j-1]) / (dy[j])**2

        visc_term = (mu + mut[i-1, j])*d2u_dy2 + (
            (mut[i-1, j+1] - mut[i-1, j-1])/(2*dy[j])
            ) * du_dy
        tau_v[i-1, j] = mu*d2u_dy2
        tau_t[i-1, j] = mut[i-1, j]*d2u_dy2 + (
            (mut[i-1, j+1] - mut[i-1, j-1])/(2*dy[j])
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
        v[i, j-1] = v[i, j] - dy[j] * du_dx

dudy_w = (u[100, 1] - u[100, 0]) / dy[0]
tau_w = mu*dudy_w
u_tau = np.sqrt(tau_w / rho)

y_plus = y*u_tau / nu
x_plus = u[100, :] / u_tau


plt.figure()
plt.plot(u[0, :], y)
plt.show(block=False)
plt.plot(u[10, :], y)
plt.show(block=False)
plt.plot(u[100, :], y)
plt.show(block=False)
plt.plot(u[1000, :], y)
plt.show(block=False)
plt.plot(u[10000, :], y, label='1')
plt.show(block=False)
plt.plot(u[40000, :], y, label='4')
plt.show(block=False)
plt.plot(u[80000-1, :], y, label='8')
plt.show(block=False)
plt.plot(u[160000-1, :], y, label='16')
plt.title("Boundary Layer 40y (with y and delta functions + 0.01Pa)")
plt.ylabel('y position')
plt.xlabel('Velocity distribution')
plt.legend()
plt.show(block=False)

plt.figure()
plt.plot(y_plus, x_plus)
plt.xscale('log')
plt.title('Plot for x=100')
plt.xlabel('y+')
plt.ylabel('u+')
plt.show(block=False)

plt.figure()
plt.plot(y_plus, tau_v[40000], label=fr'$\tau_{{\nu}}$')
plt.plot(y_plus, tau_t[40000], label=fr'$\tau_{{t}}$')
plt.title(fr'$\tau_{{\nu}}$ and $\tau_{{t}}$ with $y^+$')
plt.ylabel('Shear stresses')
plt.xlabel(fr'$y^+$')
plt.legend()
plt.show()

X, Y = np.meshgrid(x, y)
plt.figure(figsize=(8, 4))
plt.contourf(X, Y, u.T, levels=50, cmap="viridis")
plt.colorbar(label="u [m/s]")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Velocity field u(x, y)")
plt.tight_layout()
plt.show()

