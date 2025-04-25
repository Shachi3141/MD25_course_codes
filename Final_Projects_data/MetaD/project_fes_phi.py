import numpy as np

# Load FES data
data = np.loadtxt('fes.dat10.dat')

phi = data[:,0]
psi = data[:,1]
fes = data[:,2]

# Get unique phi and psi grid points
phi_vals = np.unique(phi)
psi_vals = np.unique(psi)

# Reshape FES into 2D grid
fes_grid = fes.reshape(len(phi_vals), len(psi_vals))

# Project onto phi by averaging over psi (axis=1)
fes_phi = np.mean(fes_grid, axis=1)

# Save 1D projected data
np.savetxt('fes_phi.dat', np.column_stack([phi_vals, fes_phi]), fmt="%.3f")

print("Saved 1D projected FES to fes_phi.dat")

