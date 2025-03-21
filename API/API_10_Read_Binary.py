from GFE.predata import *
import numpy as np
import os
# Geometry Parameter
geoX = 9
geoY = 9
nelx = 40
nely = 40
dis_type = "2D Exponential"
BatchDir = os.path.abspath("..\\" + dis_type + "\\Batch\\")
savepath = os.path.abspath("..\\" + dis_type)
# Read Expotential Parameter
alpha = np.load(os.path.join("..\\" + dis_type,'alpha.npy'))
beta = np.load(os.path.join("..\\" + dis_type,'beta.npy'))
# Batch Parameter
Batch_Num = 20
# Change to Batch Directory
os.chdir(BatchDir)
# Read Data
Batch_dict, Stress, Strain, Modulus = \
    BatchData(geoX, geoY, nelx, nely, BatchDir, Batch_Num, alpha, beta, dis_type)
# Save Data
os.chdir(savepath)
np.save('Stress', Stress)
np.save('Strain', Strain)
np.save('Modulus', Modulus)
DictToBin(Batch_dict, 'Batch_dict')