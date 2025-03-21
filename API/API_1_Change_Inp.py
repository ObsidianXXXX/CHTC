from GFE.abqinp import *
properties = [0.2310, 0.2310, 0.3]
jprop = [1]
UEL = Uelinp(properties, jprop)
filename = '../AbaqusInfo/Test.inp'
UEL.ToUel(filename)