import json
import numpy as np
import matplotlib.pyplot as plt

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
num_kpoints = data.shape[0] // num_bands
k = data[:num_kpoints, 0]
bands = np.reshape(data[:, 1], (-1, len(k)))

# Find the index of the valence band maximum (VBM) closest to HOMO
vbm_idx = np.argmin(np.abs(bands - homo))
row, col = np.unravel_index(vbm_idx, bands.shape)
vbm = bands[row, col]

# Subtract VBM from the band energies
bands -= vbm

# Load DOS data and subtract VBM
dos_data = np.loadtxt(dos_data)
dos_data[:, 0] -= vbm

fig, (ax1, ax2) = plt.subplots(1, 2, width_ratios=[3, 1], sharey=True)

# Plot bands with adjusted energies
for i in range(num_bands):
    ax1.scatter(k, bands[i, :], linewidth=1, alpha=0.5, color='k')
ax1.set_xlim(np.min(k), np.max(k))
ax1.set_ylim(np.min(dos_data[:, 0]), np.max(dos_data[:, 0]))
ax1.set_ylabel(r"$E-E_{vbm} [eV]$", fontsize=25)

# Add symmetry points
for sym_point_i in sym_points_index:
    ax1.axvline(k[sym_point_i], linewidth=0.75, color='k', alpha=0.5)
ax1.set_xticks(k[sym_points_index])
ax1.tick_params(axis='y', labelsize=20)
ax1.set_xticklabels(sym_labels, fontsize=20)

# Plot DOS with adjusted energies
ax2.plot(dos_data[:, 1], dos_data[:, 0], linewidth=3, alpha=0.5, color='k')
ax2.set_xlabel("dos(E)", fontsize=25)
ax2.tick_params(axis='x', labelsize=20)
dos_max = np.max(dos_data[:, 1])
ax2.set_xlim(0, dos_max + 0.1 * dos_max)

# Add the VBM and CBM lines to the plot
ax1.axhline(0, linewidth=2, label="Valence band maximum", color="red")
ax2.axhline(0, linewidth=2, color="red")

# Find the index of the conduction band minimum (CBM) closest to LUMO
cbm_idx = np.argmin(np.abs(bands - lumo))
row, col = np.unravel_index(cbm_idx, bands.shape)
cbm = bands[row, col]

# Add the CBM line to the plot
ax1.axhline(cbm - vbm, linewidth=2, label="Conduction band minimum", color="blue")
ax2.axhline(cbm - vbm, linewidth=2, color="blue")

# Compute the band gap
band_gap = cbm - vbm

# Add a box showing the band gap value
props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
ax1.text(np.min(k) + 0.1 * (np.max(k) - np.min(k)), (cbm - vbm) / 2, 
         f"Band gap\n{band_gap:.3f} eV", fontsize=15, bbox=props)

ax1.legend(fontsize=15, framealpha=0.9)
plt.show()

