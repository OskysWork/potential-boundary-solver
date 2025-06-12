import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


Nx = 160000
Ny = 80
delta0 = 0.05
L = 0.4 

x = np.linspace(0, L, Nx)
#y = delta0 * (np.tanh(2 * np.linspace(0, 1, Ny)) / np.tanh(2))
y = np.linspace(0, delta0, Ny)
dx = x[1] - x[0]
print(dx)
dy = y[1] - y[0]
print(dy**2)

rho = 1
mu = 1.95*1e-5
nu = mu / rho
dp_dx = 0#0.01*1e5

u = np.zeros((Nx, Ny))
v = np.zeros((Nx, Ny))
nut = np.zeros((Nx, Ny))
mut = np.zeros((Nx, Ny))
kappa = 0.41

count = 0

for j in range(Ny):
    eta = y[j]/delta0
    u[0, j] = 1.5 * eta - 0.5 * eta**3 if eta<1 else 1


for i in tqdm(range(1, Nx)):
    #delta = delta0 + 0.005 * x[i]     # Placeholder layer growth relation
    delta = delta0 + np.sqrt(x[i]*mu*rho)
    count+=1
    #print(count)


    eta = y / delta
    nut[i-1, :] = np.where(eta < 1, 0.01 * (eta*(1-eta)), 0)
    mut[i-1, :] = rho * nut[i-1, :]

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

dudy_w = (u[100, 1] - u[100, 0]) / dy
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
plt.title("Boundary Layer 40y (with y and delta functions + 0.01bar)")
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

