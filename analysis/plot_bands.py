import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Load JSON configuration
with open("config.json", "r") as f:
    config = json.load(f)

dos_data = config["dos_data"]
bands_data = config["bands_data"]
sym_points_index = config["sym_points_index"]
sym_labels = config["sym_labels"]
num_kpoints = config["num_kpoints"]
num_bands = config["num_bands"]
homo = config["homo"]

# Load band structure data
data = np.loadtxt(bands_data)
k = data[:num_kpoints, 0]
bands = np.reshape(data[:, 1], (-1, len(k)))

fig, (ax1, ax2) = plt.subplots(1, 2, width_ratios=[3, 1], sharey=True)

for i in range(num_bands):
    ax1.scatter(k, bands[i, :], linewidth=1, alpha=0.5, color='k')
ax1.set_xlim(np.min(k), np.max(k))
ax1.set_ylim(homo-5, homo+8)
ax1.set_ylabel("Energy (eV)", fontsize=25)

# Add symmetry points
for sym_point_i in sym_points_index:
    ax1.axvline(k[sym_point_i], linewidth=0.75, color='k', alpha=0.5)
sym_points = k[sym_points_index]
ax1.set_xticks(sym_points)
ax1.set_xticklabels(sym_labels)
ax1.axhline(homo, linewidth=2, label="E Fermi")

data = np.loadtxt(dos_data)
ax2.plot(data[:, 1], data[:, 0], linewidth=3, alpha=0.5, color='k')
ax2.set_ylim(homo-5, homo+8)
ax2.set_xlabel("dos(E)", fontsize=25)
ax2.axhline(homo, linewidth=2, label="E Fermi")


plt.show()

