import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

angles = np.array([0, 1, 2, 3, 4, 5, 6])   # degrees
velocities = np.array([5,5.46,5.85,6.24,6.63,7,7.4,7.8,8.2,8.6,9,9.3,9.75,10.15,10.5,11,11.3,11.7,12,12.5,12.9,13.3,13.65,14,14.4,14.8,15.2,15.6,16,16.4,16.8,17.15,17.5,18,18.3,18.7,19.1,19.5,19.9,20.3,20.6,21,21.45,21.84,22.23])    # m/s

# case matrix: shape (n_angles, n_velocities)
# 0 = no tone, 1 = tone, 2 = hump+tone
cases = np.array([
    [0,0,0,0,0,1,1,1,1,1,2,2,2,1,1,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,0],
    [0,0,0,0,0,0,0,0,0,0,1,2,2,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,0],
    [0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
])

data_c=pd.read_excel('Desquesnes.xlsx', engine='openpyxl',header=None)

rho = 1.225       # kg/m^3
mu = 1.81e-5      # kg/(m s)
c = 0.135          # chord [m] ← CHANGE IF NEEDED

# Compute Reynolds numbers
Re = rho * velocities * c / mu

# Marker styles
styles = {
    0: dict(marker='^', facecolors='none', edgecolors='lightgray', label='No tone'),
    1: dict(marker='o', facecolors='none', edgecolors='black', label='Tone'),
    2: dict(marker='o', facecolors='black', edgecolors='black', label='Hump + tone')
}

plt.figure(figsize=(7,6))

# Plot each point
for i, alpha in enumerate(angles):
    for j, Re_val in enumerate(Re):
        case = cases[i, j]
        plt.scatter(
            Re_val,
            alpha,
            s=60,
            linewidths=1,
            **styles[case]
        )

# Axis formatting
plt.plot(10**data_c[0],data_c[1],'k')
plt.xscale('log')
plt.xlabel('Reynolds number')
plt.ylabel('Angle of attack ($^\circ$)')
plt.ylim(-0.5, 14)

# Custom legend (avoid duplicates)
handles = []
labels = []
for k, v in styles.items():
    h = plt.scatter([], [], s=60, linewidths=1, **v)
    handles.append(h)
    labels.append(v['label'])

plt.legend(handles, labels, loc='upper left', frameon=False,ncol=3,)

plt.tight_layout()
plt.savefig('discrete_tones.png',dpi=600)
