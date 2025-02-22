import os
import numpy as np
import matplotlib.pyplot as plt 

# Collect file paths
file_paths = []
for file in os.listdir("data_RINGS/angles/"):
        if ("angle" in file) and (".dat" in file):
                file_paths.append(os.path.join("data_RINGS/angles/", file))

# Compute ADFs
weight_TeTe = 0.170627 
weight_GeTe = 2.589158 
weight_GeGe = 1.321388 

adfTeTeTe = np.loadtxt("data_RINGS/angles/angle_Te_Te_Te.dat")[:, 1] * weight_TeTe**2
adfTeTeGe = np.loadtxt("data_RINGS/angles/angle_Te_Te_Ge.dat")[:, 1] * (weight_TeTe * weight_GeTe)
adfGeTeTe = np.loadtxt("data_RINGS/angles/angle_Ge_Te_Te.dat")[:, 1] * (weight_TeTe * weight_GeTe)
adfGeTeGe = np.loadtxt("data_RINGS/angles/angle_Ge_Te_Ge.dat")[:, 1] * weight_GeTe**2
adfTe = adfTeTeTe + adfTeTeGe + adfGeTeTe + adfGeTeGe 

adfGeGeGe = np.loadtxt("data_RINGS/angles/angle_Ge_Ge_Ge.dat")[:, 1] * weight_GeGe**2
adfGeGeTe = np.loadtxt("data_RINGS/angles/angle_Ge_Ge_Te.dat")[:, 1] * (weight_GeGe * weight_GeTe)
adfTeGeGe = np.loadtxt("data_RINGS/angles/angle_Te_Ge_Ge.dat")[:, 1] * (weight_GeGe * weight_GeTe)
adfTeGeTe = np.loadtxt("data_RINGS/angles/angle_Te_Ge_Te.dat")[:, 1] * weight_GeTe**2
adfGe = adfGeGeGe + adfGeGeTe + adfTeGeGe + adfTeGeTe

angles = np.loadtxt("data_RINGS/angles/angle_Ge_Te_Te.dat")[:, 0]

# Normalize computed ADFs
adfTe /= np.max(adfTe)
adfGe /= np.max(adfGe)

# Load reference data
angleGe_ref, adfGe_ref = np.loadtxt("adfGe_ref.dat", unpack=True, delimiter=',')
angleTe_ref, adfTe_ref = np.loadtxt("adfTe_ref.dat", unpack=True, delimiter=',')

# Normalize reference ADFs
adfTe_ref /= np.max(adfTe_ref)
adfGe_ref /= np.max(adfGe_ref)

# Plot comparison
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

axes[1].scatter(angles, adfTe, color='red', marker='x', label='Computed ADF Te')
axes[1].scatter(angleTe_ref, adfTe_ref, color='blue', marker='o', label='Reference ADF Te')
axes[1].set_ylabel("ADF Te [a.u.]", fontsize=20)
axes[1].set_xlabel("Angle (degrees)", fontsize=20)
axes[1].legend(fontsize=20)
axes[1].grid(True)
axes[1].xaxis.set_tick_params(labelsize=15)
axes[1].yaxis.set_tick_params(labelsize=15)

axes[0].scatter(angles, adfGe, color='red', marker='x', label='Computed ADF Ge')
axes[0].scatter(angleGe_ref, adfGe_ref, color='blue', marker='o', label='Reference ADF Ge')
axes[0].set_ylabel("ADF Ge [a.u.]", fontsize=20)
axes[0].legend(fontsize=20)
axes[0].grid(True)
axes[0].yaxis.set_tick_params(labelsize=15)

plt.tight_layout()
plt.show()

