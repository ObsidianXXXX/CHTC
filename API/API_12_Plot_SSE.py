from GFE.predata import *
dis_type = "2D Exponential"
Batch_dict = BinToDict("..\\" + dis_type + "\\Batch_dict")
# Geometry Parameter
geoX = 9
geoY = 9
nelx = 40
nely = 40
# Index
Index = 11
Batch_test = Batch_dict[str(Index + 1)]
Batch_test.PltGptStress()
Batch_test.PltGptStrain()
Batch_test.PltNodeModulus()
