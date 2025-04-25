#!/usr/bin/env python3
"""
Ramachandran Plot Generator for Alanine Dipeptide (OPES Simulation)
Creates a 2D free energy surface with contour lines from OPES data
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata

# Configuration
CV_FILE = "COLVAR"          # OPES collective variables output
FES_FILE = "fes.dat"        # Free energy surface data
OUTPUT_FILE = "opes_ramachandran.png"

# Load data
data = np.loadtxt(CV_FILE, skiprows=1)
phi = data[:,1]  # Phi angles (radians)
psi = data[:,2]  # Psi angles (radians)

try:
    # Try to load pre-computed FES if available
    fes_data = np.loadtxt(FES_FILE)
    xx, yy, f = fes_data[:,0], fes_data[:,1], fes_data[:,2]
    grid_x, grid_y = np.mgrid[-np.pi:np.pi:100j, -np.pi:np.pi:100j]
    f = griddata((xx, yy), f, (grid_x, grid_y), method='cubic')
    xx, yy = grid_x, grid_y
    energy_label = "Free Energy (kJ/mol)"
except:
    # Fallback to kernel density estimate
    print("No FES file found, using kernel density estimate")
    xx, yy = np.mgrid[-np.pi:np.pi:100j, -np.pi:np.pi:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    kernel = gaussian_kde(np.vstack([phi, psi]))
    f = np.reshape(kernel(positions).T, xx.shape)
    energy_label = "Probability Density"

# Create figure
plt.figure(figsize=(10, 8), dpi=300)
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.labelsize': 14
})

# Enhanced colormap (similar to PLUMED)
cmap = plt.cm.viridis.reversed()

# Plot free energy surface
heatmap = plt.pcolormesh(xx, yy, f, shading='auto', cmap=cmap,
                        norm=colors.PowerNorm(gamma=0.4))

# Add contour lines (every 5 kJ/mol if FES available)
contour_levels = np.arange(0, np.nanmax(f), 5) if 'fes_data' in locals() else 10
contour = plt.contour(xx, yy, f, levels=contour_levels,
                     colors='white', linewidths=0.8, alpha=0.7)
plt.clabel(contour, inline=True, fontsize=9, fmt='%d')

# Colorbar and labels
cbar = plt.colorbar(heatmap)
cbar.set_label(energy_label, labelpad=15)
plt.xlabel(r'$\phi$ Dihedral Angle (rad)', fontsize=14)
plt.ylabel(r'$\psi$ Dihedral Angle (rad)', fontsize=14)
plt.title('Alanine Dipeptide Free Energy Surface\nOPES Simulation', pad=20)

# Secondary structure markers
structures = {
    r'$\alpha_R$-helix': (-0.9, -0.7, 'magenta', 'P'),
    r'$\beta$-sheet': (-1.4, 1.4, 'red', '*'),
    r'$\alpha_L$-helix': (0.9, 0.7, 'cyan', 'o'),
    'PPII helix': (-1.3, 0.7, 'yellow', 's')
}

for label, (x, y, color, marker) in structures.items():
    plt.scatter(x, y, c=color, marker=marker, s=120, edgecolor='black', linewidth=0.5, label=label)

plt.legend(loc='upper right', fontsize=10)

# Adjust axes
plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
          [r'$-\pi$', r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$'])
plt.yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
          [r'$-\pi$', r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$'])
plt.xlim(-np.pi, np.pi)
plt.ylim(-np.pi, np.pi)

# Save high-quality output
plt.tight_layout()
plt.savefig(OUTPUT_FILE, bbox_inches='tight', dpi=300, transparent=True)
print(f"OPES Ramachandran plot saved as '{OUTPUT_FILE}'")

# Optional: Show plot if running interactively
# plt.show()
