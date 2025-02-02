import json
import numpy as np
import matplotlib.pyplot as plt

# Load JSON configuration
with open("config_ph.json", "r") as f:
    config = json.load(f)

phonons_data = config["phonons_data"]
sym_points_index = config["sym_points_index"]
sym_labels = config["sym_labels"]
num_kpoints = config["num_kpoints"]
num_bands = config["num_bands"]

# Load band structure data
data = np.loadtxt(phonons_data)
k = data[:, 0]
phs = []
for band in range(num_bands):
    phs.append(data[:, 1+band])

plt.figure()

for band, ph in enumerate(phs):
    plt.scatter(k, ph, linewidth=1, alpha=0.5, color='k')
plt.xlim(np.min(k), np.max(k))
plt.ylabel(f"Phonon Frequency [$cm^{-1}$]", fontsize=25)

# Add symmetry points
for sym_point_i in sym_points_index:
    plt.axvline(k[sym_point_i], linewidth=0.75, color='k', alpha=0.5)
sym_points = k[sym_points_index]
plt.yticks(fontsize=15)
plt.xticks(ticks=sym_points, labels=sym_labels, fontsize=15)
plt.show()

