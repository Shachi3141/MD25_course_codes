#!/usr/bin/env python3
"""
Calculate Free Energy Surface from OPES COLVAR data
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Load data
data = np.loadtxt('COLVAR', comments='#', usecols=(1,2))  # Load phi and psi columns
phi, psi = data[:,0], data[:,1]

# Create grid
xgrid = np.linspace(-np.pi, np.pi, 100)
ygrid = np.linspace(-np.pi, np.pi, 100)
xx, yy = np.meshgrid(xgrid, ygrid)

# Calculate 2D kernel density estimate
kde = gaussian_kde(np.vstack([phi, psi]))
positions = np.vstack([xx.ravel(), yy.ravel()])
f = np.reshape(kde(positions).T, xx.shape)

# Convert to free energy (kT units)
f = -np.log(f)
f -= np.min(f)  # Set minimum to zero

# Save FES data
np.savetxt('fes.dat', np.column_stack((xx.ravel(), yy.ravel(), f.ravel())))

print("Free energy surface saved to fes.dat")
