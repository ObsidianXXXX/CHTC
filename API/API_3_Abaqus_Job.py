import sys
import os
from inspect import getsourcefile
strPythonPath, _ = os.path.split(os.path.abspath(getsourcefile(lambda:0)))
sys.path.append('%s' %(strPythonPath))
os.chdir(strPythonPath)
os.chdir('..\\')
from GFE.odbinfo import *
# Parameter Settings
## Batch Parameters
EleNum = 1600
Batch_Num = 20
## Geometry Parameters
geoX = 9.0
geoY = 9.0
## Subroutine File
subroutine = os.path.join(os.getcwd(),"AbaqusInfo\\UEL_Exy.for")
## Abaqus Setting
cpus = 1
dis_type = "2D Exponential"
BatchDir = os.path.join(os.getcwd(),dis_type + "\\Batch\\")

# Run Batch Inp
job_list = []
broken_list = []
for i in range(Batch_Num):
    if i != 0 and i % 15 == 0 :
        for job in job_list:
            mdb.jobs[job].waitForCompletion()
            if mdb.jobs[job].status == ABORTED:
                broken_list.append(int(job))
            del mdb.jobs[job]
        job_list = []
    job_list.append(str(i + 1))
    RunJob(i, BatchDir, subroutine, cpus)

if job_list is not None:
    for job in job_list:
        mdb.jobs[job].waitForCompletion()
        if mdb.jobs[job].status == ABORTED:
                broken_list.append(int(job))
        del mdb.jobs[job]

# np.save(os.path.join(BatchDir, 'broken_list'), broken_list)
# execfile('...\\AIP\\API_3_Abaqus_Job.py')
# execfile('D:\\CHTC\\API\\API_3_Abaqus_Job.py')
