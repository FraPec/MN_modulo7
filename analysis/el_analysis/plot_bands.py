import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Load JSON configuration
with open("config_el.json", "r") as f:
    config = json.load(f)

dos_data = config["dos_data"]
bands_data = config["bands_data"]
sym_points_index = config["sym_points_index"]
sym_labels = config["sym_labels"]
num_bands = config["num_bands"]
homo = config["homo"]
lumo = config["lumo"]

# Load band structure data
data = np.loadtxt(bands_data)
num_kpoints = data.shape[0]//num_bands
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
ax1.tick_params(axis='y', labelsize=20)
ax1.set_xticklabels(sym_labels, fontsize=20)

data = np.loadtxt(dos_data)
ax2.plot(data[:, 1], data[:, 0], linewidth=3, alpha=0.5, color='k')
ax2.set_ylim(homo-5, homo+8)
ax2.set_xlabel("dos(E)", fontsize=25)
ax2.tick_params(axis='x', labelsize=20)

# Find the index of the valence band maximum (VBM) closest to HOMO
vbm_idx = np.argmin(np.abs(bands - homo))
# Unravel col and row in the bands matrix
row, col = np.unravel_index(vbm_idx, bands.shape)
# Value
vbm = bands[row, col]
# Add the vbm line to the plot to show the gap
ax1.axhline(vbm, linewidth=2, label="Valence band maximum", color="red")
ax2.axhline(vbm, linewidth=2, label="Valence band maximum", color="red")

# Find the index of the conduction band minimum (CBM) closest to LUMO
cbm_idx = np.argmin(np.abs(bands - lumo))
# Unravel col and row in the bands matrix
row, col = np.unravel_index(cbm_idx, bands.shape)
# Value
cbm = bands[row, col]
# Add the vbm line to the plot to show the gap
ax1.axhline(cbm, linewidth=2, label="Conduction band minimum", color="blue")
ax2.axhline(cbm, linewidth=2, label="Conduction band minimum", color="blue")

# Compute the band gap
band_gap = cbm - vbm

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.65, 0.5, f"Gap = {band_gap:.3f} eV",
         fontsize=15, transform=ax1.transAxes, verticalalignment='top', bbox=props)

ax1.legend(fontsize=15)

plt.show()

