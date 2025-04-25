import numpy as np

# Load FES data
data = np.loadtxt('fes.dat10.dat')

phi = data[:, 0]
psi = data[:, 1]
fes = data[:, 2]

# Get unique phi and psi grid points
phi_vals = np.unique(phi)
psi_vals = np.unique(psi)

# Reshape FES into 2D grid
fes_grid = fes.reshape(len(phi_vals), len(psi_vals))

# Project onto psi by averaging over phi (axis=0)
fes_psi = np.mean(fes_grid, axis=0)

# Save 1D projected data
np.savetxt('fes_psi.dat', np.column_stack([psi_vals, fes_psi]), fmt="%.3f")

print("Saved 1D projected FES to fes_psi.dat")
