"""
Copyright © 2024 Gengxuan Zhu
    Version June. 2024
                    Gengxuan Zhu   zhugx@zju.edu.cn
        Institute of Applied Mechanics, Zhejiang University

Description:
    Class Uelinp:
        This class is used to create a new input file for UEL.

    Functions:
        ChangeInp(properties,jproperties,lines,BatchDir):
            This function is used to change the properites of the uel input file.

        Getlines(InpName):
            This function is used to get the information of the input file
            and create a new directory for the batch.

        CreateDir(BatchDir):
            This function is used to create a new directory.

        BatchInp(seed, Batch_Num, Batch_Dir, InpName, E1, E2, geoY, Poisson, Thickness, nEl):
            This function is used to create a new input file for each batch.
    
    Usage:
        uelinp = Uelinp(properties,jproperties)
        uelinp.ToUel(InpName)

        """
import os
import numpy as np
class Uelinp:
    """This class is used to create a new input file for UEL.
    
    Args:
        properties (np.array): [E_1, E_2,Poisson's ratio, thickness]
        jproperties (np.array): [Number of Elements, batch number]
        
    Attributes:
        properties (np.array): [E_1, E_2,Poisson's ratio, thickness]
        jproperties (np.array): [Number of Elements, batch number]  
        nodes (int): number of nodes in the element.
        coordinates (int): number of coordinates in the element.
        type (str): type of the element.
        variables (int): number of field variable outputs.
        nEl (int): number of elements in the input file.
        uelpath (str): absolute path of the new input file.
        
    Methods:
        ToUel(InpName): This function is used to create a new input file for UEL.
        
        GetInpInfo(InpName): This function is used to get the number of elements
        and coordinates from the input file.
                                
    """    
    def __init__(self,properties,jproperties):
        self._properties = properties
        self._jproperties = jproperties
        self._nodes = 4
        self._coordinates = 2
        self._type = "U24"
        self._variables = 32
    
    @property
    def properties(self):
        return self._properties
    
    @property
    def jproperties(self):
        return self._jproperties
    
    @property
    def nodes(self):
        return self._nodes
    
    @property
    def type(self):
        return self._type

    @property
    def coordinates(self):
        return self._coordinates
    
    @property
    def variables(self):
        return self._variables

    @property
    def nEl(self):
        return self._nEl
    
    @property
    def uelpath(self):
        return self._uelpath

    def ToUel(self,InpName):
        """This function is used to create a new input file for UEL.

        Args:
            InpName (str): absolute path of the input file.
        """        
        # Read Initial Inp
        lines = self.GetInpInfo(InpName)
        print(self.jproperties)
        # Create New Inp File for UEL
        uelfile = Uelinp.ChInpName(InpName)
        self._uelpath = uelfile
        # Target Line
        EleType = "*Element, type=CPS4\n"
        newEle = "*Element, type=%s\n" % (self.type)
        # userEle Line
        userEle = "*User element, nodes=%d, type=%s, Iproperties=%d, properties=%d, coordinates=%d, variables=%d\n" \
                        % (self.nodes,self.type,len(self.jproperties),len(self.properties),self.coordinates,self.variables)
        PartEnd = "*End Part\n"
        strindex = lines.index(EleType)
        lines.insert(strindex, userEle)
        # Add Free Degree Line
        if self.coordinates == 3:
            lines.insert(strindex + 1,"1,2,3\n")
        else:
            lines.insert(strindex + 1,"1,2\n")
        # Add Element Line
        lines.insert(strindex + 2, newEle)
        # Add Element Property Line
        strindex = lines.index(PartEnd)
        uelProp = "*Uel property, elset=Set-1\n"
        lines.insert(strindex, uelProp)
        # Add Element Property Value Line
        prop_str = ""
        for i in self.properties:
            prop_str += "%5.4f, " % (i)
        for i in self.jproperties:
            prop_str += "%5d," % (i)
        prop_str += "\n"
        lines.insert(strindex + 1, prop_str)
        newinp = open(uelfile,"w")
        # Delete CPS4 Lines
        Section = "*Solid section, elset=Set-1, material=Material-1\n"
        Material = "*Material, name=Material-1\n"
        Elastic = "*Elastic\n"
        s1 = lines.index(Section)
        s2 = lines.index(Material)
        s3 = lines.index(Elastic)
        s4 = lines.index(EleType)
        for newline in lines:
            id = lines.index(newline)
            if id == s1 or id == s2 or id == s3 or id == s3 + 1 or id == s4:
                continue
            else:
                newinp.write(newline)
        newinp.close()
    
    def GetInpInfo(self,InpName):
        """This function is used to get the number of elements 
            and coordinates from the input file.

        Args:
            InpName (str): absolute path of the input file.
            
        Returns:
            lines (list): lines of the input file.
        
        """
        # Open CPS4 Inp File
        with open(InpName,'r',encoding='utf-8') as f:
            lines = f.readlines()
        # Get Number of Elements and Coordinates
        EndEle = "*Nset, nset=Set-1, generate\n"
        strindex = lines.index(EndEle)
        EleLine = lines[strindex-1]
        self._nEl = int(EleLine.split(',')[0].strip())
        self._jproperties = [self.nEl] + self.jproperties
        return lines

    @staticmethod
    def ChInpName(InpName):
        """This function is used to change the name of the input file
            by adding "UEL_" before the original file name.

        Args:
            InpName (str): absolute path of the input file.
        
        """
        path_split = os.path.split(InpName)
        filepath, filename = path_split[0], path_split[1]
        filename_split = os.path.splitext(filename)
        filename_ini, extension = 'UEL_' + filename_split[0], filename_split[1]
        newfile = os.path.join(filepath, filename_ini + extension)
        return newfile

def ChangeInp(properties,jproperties,lines,BatchDir):
    """This function is used to change the properites of the uel input file.

    Args:
        properties (np.array): [E_1, E_2,Poisson's ratio, thickness]
        jproperties (np.array): [Number of Elements, batch number]
        lines (list): lines of the input file.
        BatchDir (str): absolute path of the batch directory.    Batch\
    """    
    # Get New Property String
    prop_str = ""
    for i in properties:
        prop_str += "%5.4f, " % (i)
    for i in jproperties:
        prop_str += "%5d," % (i)
    prop_str += "\n"
    # Replace Target String
    PartEnd = "*End Part\n"
    strindex = lines.index(PartEnd)
    tarindex = strindex - 1
    lines[tarindex] = prop_str
    # Move to New Simulation File  Batch\jproperties[1]
    abaWorkDir = os.path.join(BatchDir,'%d' % (jproperties[1]))
    CreateDir(abaWorkDir)
    os.chdir(abaWorkDir)
    # Create New Inp & Write    Batch\jproperties[1]\jproperties[1].inp
    NewInp = "%d.inp" % (jproperties[1])
    NewPath = os.path.join(abaWorkDir,NewInp)
    NewFile = open(NewPath,'w')
    NewFile.write(''.join(lines))
    NewFile.close()

def GetLines(OldInp):
    """_summary_

    Args:
        OldInp (str): initial inp file name    UEL_test.inp

    Returns:
        lines (list): lines of the input file.
    """    
    # Read Initial Inp
    Cur_Dir = os.getcwd()
    OldPath = os.path.join(Cur_Dir,OldInp)
    OldFile = open(OldPath,'r')
    lines = OldFile.readlines()
    # Create Target Filepath 
    Batch_Dir = os.path.join(Cur_Dir,'Batch')
    CreateDir(Batch_Dir)
    OldFile.close()
    return lines

def CreateDir(file):
    """This function is used to create a new directory.

    Args:
        file (str): absolute path of the directory.
    """    
    if not os.path.exists(file):
        os.makedirs(file)
    else:
        print('The directory already exists.')

def BatchExpX(seed, Batch_Num, Batch_Dir, InpName, E1, E2, geoY, Poisson, Thickness, nEl):
    """This function is used to create a batch of input files.
    
    Args:
        seed (int): random seed.
        Batch_Num (int): number of the batch.
        Batch_Dir (str): absolute path of the batch directory.
        InpName (str): initial inp file name.
        E1 (float): minimum Young's modulus of all batch.
        E2 (float): maximum Young's modulus of all batch.
        geoY (float): length of the y direction
        Poisson (float): Poisson's ratio of all batch.
        Thickness (float): thickness of all batch.
        nEl (int): number of elements of all batch.
        
    """
    # Create Young's Modulus Array
    E = np.linspace(E_min, E_max, Batch_Num)
    np.random.seed(seed)
    np.random.shuffle(E)
    np.save('E.npy',E)
    beta = np.log(E / E_min) / geoY
    np.save('beta.npy',beta)
    # Create Batch Directory
    CreateDir(Batch_Dir)
    # Read Initial UEL Inp
    lines = GetLines(InpName)
    # Create Batch Inp
    for i in range(Batch_Num):
        prop = [E1, E[i], Poisson, Thickness]
        jprop = [nEl, i]
        ChangeInp(prop, jprop, lines, Batch_Dir)

def BatchInp(seedX, seedY, Batch_Num, Batch_Dir, InpName, E1, E2, geoX, geoY, Poisson, nEl, dis_type):
    """This function is used to create a batch of input files.
    
    Args:
        seedX (int): random seed.
        seedY (int): random seed.
        Batch_Num (int): number of the batch.
        Batch_Dir (str): absolute path of the batch directory.
        InpName (str): initial inp file name.
        E1 (float): maximum Young's modulus of all batch along X direction.
        E2 (float): maximum Young's modulus of all batch along Y direction.
        geoX (float): length of the x direction
        geoY (float): length of the y direction
        Poisson (float): Poisson's ratio of all batch.
        nEl (int): number of elements of all batch.
        dis_type (str): Type of distribution, '2D Linear' or '2D Exponential'
        
    """
    # Create Young's Modulus Array
    quart = Batch_Num // 4
    if dis_type == '2D Linear':
        ## X - Y -
        Ex1 = np.linspace(0.51, 1.0, quart)
        Ey1 = np.linspace(0.51, 1.0, quart)
        ## X - Y +
        Ex2 = np.linspace(0.51, 1.0, quart)
        Ey2 = np.linspace(1.0, E2, quart)
        ## X + Y -
        Ex3 = np.linspace(1.0, E1, quart)
        Ey3 = np.linspace(0.51, 1.0, quart)
        ## X + Y +
        Ex4 = np.linspace(1.0, E1, quart)
        Ey4 = np.linspace(1.0, E2, quart)
        ## Integrate
        Ex = np.concatenate((Ex1, Ex2, Ex3, Ex4))
        Ey = np.concatenate((Ey1, Ey2, Ey3, Ey4))
        # Shuffle Young's Modulus Array
        np.random.seed(seedX)
        np.random.shuffle(Ex)
        np.random.seed(seedY)
        np.random.shuffle(Ey)
        # Calculate Bilinear Parameter
        alpha = (Ex - 1.0) / geoX
        beta = (Ey - 1.0) / geoY
    elif dis_type == '2D Exponential':
        ## X - Y -
        Ex1 = np.linspace(0.1, 1.0, quart)
        Ey1 = np.linspace(0.1, 1.0, quart)
        ## X - Y +
        Ex2 = np.linspace(0.1, 1.0, quart)
        Ey2 = np.linspace(1.0, E2, quart)
        ## X + Y -
        Ex3 = np.linspace(1.0, E1, quart)
        Ey3 = np.linspace(0.1, 1.0, quart)
        ## X + Y +
        Ex4 = np.linspace(1.0, E1, quart)
        Ey4 = np.linspace(1.0, E2, quart)
        # Integrate
        Ex = np.concatenate((Ex1, Ex2, Ex3, Ex4))
        Ey = np.concatenate((Ey1, Ey2, Ey3, Ey4))
        # Shuffle Young's Modulus Array
        np.random.seed(seedX)
        np.random.shuffle(Ex)
        np.random.seed(seedY)
        np.random.shuffle(Ey)
        # Calculate Expotential Parameter
        alpha = np.log(Ex / 1.0) / geoX
        beta = np.log(Ey / 1.0) / geoY

    else:
        raise ValueError('Invalid Distribution Type')
    # Save Young's Modulus Array
    os.chdir("..\\" + dis_type)
    np.save('Ex.npy', Ex)
    np.save('Ey.npy', Ey)
    # Save Exponential Parameter
    np.save('alpha.npy', alpha)
    np.save('beta.npy', beta)
    # Create Batch Directory
    CreateDir(Batch_Dir)
    # Read Initial UEL Inp
    lines = GetLines(InpName)
    # Create Batch Inp
    for i in range(Batch_Num):
        prop = [alpha[i], beta[i], Poisson]
        jprop = [nEl, i + 1]
        ChangeInp(prop, jprop, lines, Batch_Dir)