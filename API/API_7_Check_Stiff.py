from scipy import sparse
import os

dis_type = "2D Exponential"
BatchDir = os.path.abspath("..\\" + dis_type + "\\Batch\\")
os.chdir(BatchDir)
Stiff = sparse.load_npz('1Stiff.npz')
print(Stiff)