# from abaqus import *
# from odbAccess import *
# from textRepr import *
# from abaqusConstants import *
import os
import sys
import numpy as np
from inspect import getsourcefile
strPythonPath, strFileName = os.path.split(os.path.abspath(getsourcefile(lambda:0)))
sys.path.append('%s' %(strPythonPath))
dis_type = "2D Exponential"
BatchDir = os.path.abspath("..\\" + dis_type + "\\Batch\\")
savepath = os.path.abspath("..\\" + dis_type)

from GFE.odbinfo import *

nodsx = 41
nodsy = 41
Batch_Num = 20
U, disp = BatchDisp(nodsx,nodsy,BatchDir,Batch_Num)
os.chdir(savepath)
np.save('U',U)
np.save('disp', disp)

# execfile('..\\CHTC\\AIP\\API_10_Read_Odb.py')
