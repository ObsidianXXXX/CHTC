#coding=UTF-8
"""
Copyright © 2024 Gengxuan Zhu
    Version June. 2024
                    Gengxuan Zhu   zhugx@zju.edu.cn
        Institute of Applied Mechanics, Zhejiang University

Description:
    Functions:
        GetNodeLabel(allnode):
            This function is used to get the node label of the model.
            
        GetEleLabel(allele):
            This function is used to get the element label of the model.
            
        GetFieldValue(field_data,Type):
            This function is used to get the field value of t

        DispFrameValue(cps_ins,uel_frame,cps_frame,node_label):
            This function is used to create a new field output for the displacement of the model.

        GetDisp(nodsx,nodsy,U_temp):
            This function is used to get the displacement of the model. 

        BatchDisp(nodsx,nodsy,Batch_Dir,Batch_Num):
            This function is used to get the displacement of all batches.

        AddUELData(filepath,cps_frame,cps_ins,ele_label):
            This function is used to add UEL data to ODB.

        AddBinData(stress_bin,strain_bin,cps_frame,cps_ins,ele_label,NGauss,type):
            This function is used to add UEL data from binary file to ODB.  

        DisplayShape(cps):
            This function is used to display the deformed shape of the CPS4 model.   

        RunInp(abaqLauncher,abaWorkDir,InpNum,cpus=1,subroutine=None):
            This function is used to run abaqus inp in command.

        RunJob(Batch_Num,BatchDir,subroutine=None,cpus=1):
            This function is used to submit abaqus job in abaqus.
        

"""
from abaqus import *
from odbAccess import *
from textRepr import *
from abaqusConstants import *
import visualization
import numpy as np
import os

def GetNodeLabel(allnode):
    """This function is used to get the node label of the model.

    Args:
        allnode (list): A list of all nodes in the model.   
    
    return:
        node_label (tuple): A tuple of node labels.
    """
    node_label = ()
    for node in allnode:
        node_label = node_label + (node.label,)
    return node_label

def GetEleLabel(allele):
    """This function is used to get the element label of the model.
    
    Args:
        allele (list): A list of all elements in the model.
    
    return:
        ele_label (tuple): A tuple of element labels.
    """
    ele_label = ()
    for ele in allele:
        ele_label = ele_label + (ele.label,)
    return ele_label

def GetFieldValue(field_data,Type):
    """This function is used to get the field value of the model.

    Args:
        field_data (list): A list of field data.
        Type (Constante): The type of the field data.
    
    return:
        data (tuple): A tuple of field data.
    """
    data = ()
    if Type == VECTOR:
        for value in field_data:
            data = data + ((value.data[0],value.data[1]),)
    elif Type == TENSOR_2D_PLANAR:
        for value in field_data:
            data = data + ((value.data[0],value.data[1],value.data[2],value.data[3]),)
    elif Type == SCALAR:
        for value in field_data:
            data = data + ((value.data[0]),)
    return data

def DispFrameValue(cps_ins,uel_frame,cps_frame,node_label):
    """This function is used to create a new field output for the displacement of the model.
    
    Args:
        cps_ins (Instance): The instance of the CPS4 model.
        uel_frame (list): A list of UEL frame data.
        cps_frame (list): A list of CPS4 frame data.
        node_label (tuple): A tuple of node labels.
    
    """
    i = 0
    for frame in cps_frame:
        # Get UEL FieldOutput
        U_temp = uel_frame[i].fieldOutputs['U'].values
        U_field = GetFieldValue(field_data=U_temp,Type=VECTOR)
        New_field = frame.FieldOutput(name='UEL_U',description='UEL Disp Result',type=VECTOR)
        New_field.addData(position=NODAL,instance=cps_ins,labels=node_label,data=U_field)
        i = i + 1

def GetDisp(nodsx,nodsy,U_temp):
    """This function is used to get the displacement of the model.

    Args:
        nodsx (int): The number of nodes in the x direction.
        nodsy (int): The number of nodes in the y direction.
        U_temp (list): A list of displacement data.
        
    return:
        U1 (np.2Darray): The displacement in the x direction.
        U2 (np.2Darray): The displacement in the y direction.

    """
    nod_num = nodsx * nodsy
    ndof = 2 * nod_num
    U1 = np.zeros(nod_num)
    U2 = np.zeros(nod_num)
    disp = np.zeros(ndof)
    U1[:], U2[:] = zip(*(value.data for value in U_temp))
    disp[::2], disp[1::2] = zip(*(value.data for value in U_temp))
    return U1.reshape((nodsy,nodsx)), U2.reshape((nodsy,nodsx)), disp

def BatchDisp(nodsx,nodsy,Batch_Dir,Batch_Num):
    """This function is used to get the displacement of all batches.
    
    Args:
        nodsx (int): The number of nodes in the x direction.
        nodsy (int): The number of nodes in the y direction.
        Batch_Dir (str): The directory of the batch.
        Batch_Num (int): The number of batches.
        
    return:
        U (np.3Darray): The displacement field of all batches.
        disp (np.2Darray) : The displacement vector of all batches
        
    """
    U = np.zeros((Batch_Num,nodsy,nodsx,2))
    disp = np.zeros((2 * nodsx * nodsy, Batch_Num))
    for i in range(Batch_Num):
        abqWorkDir = os.path.join(Batch_Dir,str(i + 1))
        os.chdir(abqWorkDir)
        odb_name = os.path.join(abqWorkDir,str(i + 1) + '.odb')
        odb = openOdb(odb_name)
        last_frame = odb.steps['Step-1'].frames[-1]
        U_temp = last_frame.fieldOutputs['U'].values
        U[i,:,:,0], U[i,:,:,1], disp[:, i] = GetDisp(nodsx,nodsy,U_temp)
        odb.close()
        print('Batch ' + str(i + 1) + ' is done!')
    return U, disp

def AddUELData(filepath,cps_frame,cps_ins,ele_label):
    """This function is used to add UEL data to ODB.
    
    Args:
        filepath (str): The path of the UEL data file.
        cps_frame (list): A list of frames.
        cps_ins (Instance): The instance of the CPS4 model.
        ele_label (tuple): A tuple of element labels.

    """
    # Read UEL Data
    ueldat = open(filepath)
    lines = ueldat.readlines()
    End = " ==================================\n"
    EndLine = [index for (index, item) in enumerate(lines) if item == End]
    Start_index = [index + 1 for index in EndLine[:-1]]
    End_index = [index - 1 for index in EndLine[1:]]
    for n in range(0,len(cps_frame)):
        # Create Stress & Strain Field
        frame = cps_frame[n]
        Strain_Field = frame.FieldOutput(name='UEL_E',description='UEL Strain Result',type=TENSOR_2D_PLANAR)
        Stress_Field = frame.FieldOutput(name='UEL_S',description='UEL Stress Result',type=TENSOR_2D_PLANAR)
        # Read Data from Each Frame
        data_line = lines[Start_index[n]:End_index[n]]
        i = 2
        Strain = ()
        Stress = ()
        while i < len(data_line):
            mod_i = i % 8
            if mod_i == 2 or mod_i == 3 or mod_i == 4 or mod_i == 5:
                line = np.array(data_line[i].split()).astype(np.float32)
                Strain = Strain + ((line[1], line[2], 0.0, line[3]),)
                Stress = Stress + ((line[4], line[5], 0.0, line[6]),)
                i += 1
            elif mod_i == 6:
                i += 4
            else:
                continue
        Strain_Field.addData(position=INTEGRATION_POINT,instance=cps_ins,labels=ele_label,data=Strain)
        Stress_Field.addData(position=INTEGRATION_POINT,instance=cps_ins,labels=ele_label,data=Stress)

def AddBinData(stress_bin,strain_bin,cps_frame,cps_ins,ele_label,NGauss,type):
    """This function is used to add UEL data from binary file to ODB.

    Args:
        stress_bin (str): The path of the binary file of stress data.
        strain_bin (str): The path of the binary file of strain data.
        cps_frame (list): A list of frames.
        cps_ins (Instance): The instance of the CPS4 model.
        ele_label (tuple): A tuple of element labels.
        NGauss (int): The number of Gauss points.
        type (str): The type of the stress and strain data.
    
    """
    # Read Binary File Split By Frame
    frame_num = len(cps_frame)
    stress_data = np.split(np.fromfile(stress_bin,dtype=np.float64),frame_num)
    strain_data = np.split(np.fromfile(strain_bin,dtype=np.float64),frame_num)
    num = len(ele_label) * NGauss
    # Loop to Add Stress & Strain for Each Frame
    for i in range(0,frame_num):
        # Create Stress Field & Strain Field
        frame = cps_frame[i]
        Stress_Field = frame.FieldOutput(name='UEL_S', \
            description='UEL Stress Result',type=type)
        Strain_Field = frame.FieldOutput(name='UEL_E', \
            description='UEL Strain Result',type=type)
        # Split By Gauss Point
        stress_frame = tuple(map(tuple,np.split(stress_data[i],num)))
        strain_frame = tuple(map(tuple,np.split(strain_data[i],num)))
        # Add Stress & Strain Data
        Stress_Field.addData(position=INTEGRATION_POINT, \
            instance=cps_ins,labels=ele_label,data=stress_frame)
        Strain_Field.addData(position=INTEGRATION_POINT, \
            instance=cps_ins,labels=ele_label,data=strain_frame)

def DisplayShape(cps):
    """This function is used to display the deformed shape of the CPS4 model.
    
    Args:
        cps (CPS4): The CPS4 model. 
    
    """
    # create a new viewport
    myViewport = session.Viewport(name='UEL Result', origin=(10, 10), width=150, height=100)
    # open target odb
    myOdb = visualization.openOdb(cps)
    # set default plot object in the new viewport
    myViewport.setValues(displayedObject=myOdb)
    frames = myOdb.steps['Step-1'].frames
    # set the new deformed variable to display for each frame
    for frame in frames:
        field = frame.fieldOutputs['UEL_U']
        myViewport.odbDisplay.setDeformedVariable(field)

# Run Inp File in CMD
def RunInp(abaqLauncher,abaWorkDir,InpNum,cpus=1,subroutine=None):
    """This function is used to run abaqus inp in command.

    Args:
        abaqusLauncher (str): The path of abaqus Launcher.
        abaqusWorkDir (str): The path of abaqus work directory.
        InpNum (int): The number of inp file.
        cpus (int): The number of cpus.
        subroutine (str): The subroutine of abaqus.

    """
    os.chdir(abaWorkDir)
    InpName = os.path.join(abaWorkDir,"%d.inp" %(InpNum))
    abaqus_command = [
        r'"%s" '% abaqLauncher,\
        "job=job-%d " % InpNum,\
        'input="%s" ' % InpName,\
        "cpus=" + str(cpus),\
        " interactive "]
    if subroutine is not None:
        abaqus_command.append('user="%s"' % subroutine)
    command = ''.join(abaqus_command)
    command = command.replace('\\','\\\\')
    print(command)
    t0 = os.path.getatime(__file__)
    process = os.popen(command, 'r')
    elapsed_time = os.path.getatime(__file__) - t0
    stdout = process.read()
    stdout = stdout.decode('utf-8', errors='replace')
    print("Output Message:", stdout)
    if process.close():
        print("Time:", elapsed_time, "s")
    else:
        print("Abaqus Error!!!")

# Run Inp in Abaqus
def RunJob(Batch_Num,BatchDir,subroutine=None,cpus=1):
    """This function is used to submit abaqus job in abaqus.

    Args:
        Batch_Num (int): The number of batch.
        BatchDir (str): The path of batch directory.
        subroutine (str): The subroutine of abaqus.
        cpus (int): The number of cpus.
    
    """
    # Run Inp in Each File
    ## Change to Target File
    abaWorkDir = os.path.join(BatchDir,'%d' % (Batch_Num + 1))
    os.chdir(abaWorkDir)
    ## Get Inp Name
    InpName = os.path.join(abaWorkDir,"%d.inp" %(Batch_Num + 1))
    ## Create Job From Inp
    mdb.JobFromInputFile(name=str(Batch_Num + 1), 
        inputFileName=InpName, 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, 
        userSubroutine=subroutine, 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=cpus, 
        numGPUs=0)
    ## Submit Current Job
    mdb.jobs[str(Batch_Num + 1)].submit(consistencyChecking=OFF)
    ##  Job Running
    # time_start = time.time()
    # print("Job-" + str(Batch_Num + 1) + " Running...")
    # mdb.jobs[str(Batch_Num + 1)].waitForCompletion()
    # del mdb.jobs[str(Batch_Num + 1)]
    ## Calculate Total Time
    # time_end = time.time()
    # time_use = time_end - time_start
    # print("Job-%s Completed, Total Time:%s s" % (str(Batch_Num + 1), str(time_use)))
    ## Delete Job
