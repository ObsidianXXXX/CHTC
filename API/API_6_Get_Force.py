from GFE.predata import *
from scipy import sparse
import os
# Geometry Parameter
geoX = 9
geoY = 9
nelx = 40
nely = 40
# Get Force
P = 0.1
Sample = PlFGMInfo(geoX, geoY, nelx, nely)
F = Sample.GetForce(P)
print(F)
# Save Force Vector
dis_type = "2D Exponential"
savepath = os.path.abspath("..\\" + dis_type)
os.chdir(savepath)
sparse.save_npz('Force', F, compressed=True)