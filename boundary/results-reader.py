import pandas as pd
import os
import matplotlib.pyplot as plt

file_path = os.path.join('.', 'boundary-layer.csv')

file = pd.read_csv(file_path)

y = file['y'].tolist()
u_0 = file['u[0]'].tolist()
u_10 = file['u[10000]'].tolist()
u_80 = file['u[80000]'].tolist()
u_160 = file['u[160000]'].tolist()
u_320 = file['u[320000]'].tolist()
u_640 = file['u[640000]'].tolist()
u_plus = file['u+'].tolist()
y_plus = file['y+'].tolist()
tau_v = file['tau_v[40000]'].tolist()
tau_t = file['tau_t[40000]'].tolist()
tau_v_h = file['tau_v[640000]']
tau_t_h = file['tau_t[640000]']

print(y)

print(file.columns.tolist())

plt.figure()
plt.plot(u_0, y, label='0')
plt.plot(u_10, y, label='10')
plt.plot(u_80, y, label='80')
plt.plot(u_160, y, label='160')
plt.plot(u_320, y, label='320')
plt.plot(u_640, y, label='640')
plt.ylabel('Velocity u [m/s]')
plt.xlabel('Wall distance y [m]')
plt.legend()
plt.show(block=False)

plt.figure()
plt.plot(tau_v, y_plus, label='Tau_v')
plt.plot(tau_t, y_plus, label='Tau_t')
plt.plot(tau_v_h, y_plus, label='Tau_v_h')
plt.plot(tau_t_h, y_plus, label='Tau_t_h')
plt.legend()
plt.show()