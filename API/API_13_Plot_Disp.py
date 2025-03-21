import numpy as np
import matplotlib.pyplot as plt
# Geometry Parameters
geoX = 9
geoY = 9
nodsx = 41
nodsy = 41
# Loading Data
dis_type = "2D Exponential"
U = np.load("..\\" + dis_type + "\\U.npy")
X, Y = np.meshgrid(np.linspace(0, geoX, nodsx), np.linspace(0, geoY, nodsy))
fig, axes = plt.subplots(1, 2)
# Plot the first displacement component
pc1 = axes[0].pcolormesh(X, Y, U[0,:,:,0], cmap='RdYlBu')
fig.colorbar(pc1, ax=axes[0])
axes[0].set_title("U1 Distribution")
axes[0].set_xlabel("X")
axes[0].set_ylabel("Y")
axes[0].tick_params(axis='both', which='both', length=0)
# Plot the second displacement component
pc2 = axes[1].pcolormesh(X, Y, U[0,:,:,1], cmap='RdYlBu')
fig.colorbar(pc2, ax=axes[1])
axes[1].set_title("U2 Distribution")
axes[1].set_xlabel("X")
axes[1].set_ylabel("Y")
axes[1].tick_params(axis='both', which='both', length=0)
plt.tight_layout()
plt.show()