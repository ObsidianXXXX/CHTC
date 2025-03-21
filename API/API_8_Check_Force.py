from scipy import sparse
import os

Cur_Dir = os.getcwd()
Batch_Dir = os.path.join(Cur_Dir,'Batch\\1')
os.chdir(Batch_Dir)
Force = sparse.load_npz('1Force.npz')
print(Force)