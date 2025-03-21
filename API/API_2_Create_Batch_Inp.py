from GFE.abqinp import *
import os   
# Parameter Settings
## Batch Parameters
EleNum = 1600
Batch_Num = 20
## Geometry Setting
geoX = 9.0
geoY = 9.0
## Inp Name
InpName = "../AbaqusInfo/UEL_Test.inp"
## Elastic Properties Setting
seedX = 20
seedY = 40
E1 = 5.0
E2 = 5.0
mu = 0.3
## File Path
dis_type = "2D Exponential"
CreateDir("..\\" + dis_type)
BatchDir = os.path.abspath("..\\" + dis_type + "\\Batch\\")
BatchInp(seedX,seedY,Batch_Num,BatchDir,InpName,E1,E2,geoX,geoY,mu,EleNum,dis_type)