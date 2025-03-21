from GFE.predata import *
import os
# Geometry Parameter
geoX = 9
geoY = 9
nelx = 40
nely = 40
# Read Parameter
dis_type = "2D Exponential"
alpha = np.load(os.path.join("..\\" + dis_type,'alpha.npy'))
beta = np.load(os.path.join("..\\" + dis_type,'beta.npy'))
dis_type = '2D Exponential'
# Batch Info
Batch_Num = 20
BatchDir = os.path.abspath("..\\" + dis_type + "\\Batch\\")
# Get Modulus Field & Stiffness Matrix
BatchMat(geoX, geoY, nelx, nely, BatchDir, Batch_Num, alpha, beta, dis_type)