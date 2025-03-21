#coding=UTF-8
"""
Copyright © 2024 Gengxuan Zhu
    Version June. 2024
                    Gengxuan Zhu   zhugx@zju.edu.cn
        Institute of Applied Mechanics, Zhejiang University

Description:
    Mechanical Model(Uniaxial Stretch of FGM):
              x  +----------------------------------+  ------> f
              x  |                                  |  ------> f
              x  |      |      |      |      |      |  ------> f
              x  |      |      |      |      |      |  ------> f
        geoY  x  |      |      |      |      |      |  ------> f
              x  |      |      | E(xy)|      |      |  ------> f
              x  |      |      |      |      |      |  ------> f
              x  |      |      |      |      |      |  ------> f
              x  |      V      V      V      V      |  ------> f
              x  |               geoX               |  ------> f
              x  +----------------------------------+  ------> f
    
    Class PLFGMInfo:
        The PlFGMInfo class is designed to 
            1.Store mesh information
                geoX, geoY, nelx, nely, nEl, ex, ey, nodsx, nodsy
                coord, x, y, nGauss, eleNodsID, xEle, yEle
            2.Retrieve stress tensor & strain tensor at Gauss points  
            from the binary file of the graded finite element analysis
            of plane functionally graded material(FGM).
            3.Create modulus data at gauss points according to
            E1(Minimum Modulus) & beta(Grading Ratio).  
            4.Plot each component of stress & strain tensor 
            and modulus distribution.
    
    Functions:
        BatchData(geoX, geoY, nelx, nely, Batch_Dir, Batch_Num, E1, beta):
            This function is used to get data from all batches
            and return a dictionary containing the data from all batches,
            where the key is the batch number and the value is each sample,
            a 4Darray containing S11, S22, S12 of all batches,
            a 4Darray containing E11, E22, E12 of all batches
            and a 2Darray containing modulus of all batches.
            

        DictToBin(tar_dict, filename):
            This function is used to write the data in the dictionary to a binary file.
        
        BinToDict(filename):
            This function is used to read the binary file and return the data as
            a dictionary.

"""
import numpy as np
import pickle
import os
from matplotlib import pyplot as plt
from scipy import sparse

class PlFGMInfo:
    """This class is used to store the mesh information,
        get the stress and strain data at Gauss points from 
        the binary file of the graded finite element analysis
        of plane FGM and create modulus field.
      
    Example:
        >>> Sample = PlFGMInfo(geoX, geoY, nelx, nely)
        >>> Sample.GetStress(stress_bin)
        >>> Sample.GetStrain(strain_bin)
        >>> Sample.GetModulus(E1, beta)
        >>> Sample.PlGptStress()
    
    Attributes:
        geoX (float): length of the x direction
        geoY (float): length of the y direction
        nelx (int): number of elements in the x direction
        nely (int): number of elements in the y direction
        nEl (int): number of elements
        ex (int): length of each element in the x direction
        ey (int): length of each element in the y direction
        nodsx (int): number of nodes in the x direction
        nodsy (int): number of nodes in the y direction
        coord (np.array): node coordinates
        x (np.array): x coordinates of the nodes
        y (np.array): y coordinates of the nodes
        nGauss (int): number of Gauss points
        eleNodsID (np.array): element node ID
        xEle (np.array): x coordinates of the nodes in each element
        yEle (np.array): y coordinates of the nodes in each element
        stress (dict): stress data at Gauss points
        strain (dict): strain data at Gauss points
        modulus (np.array): modulus data at Gauss points
        
    Methods:
        GetStress(): This function reads the plane stress data from the binary file
            and returns the stress data at Gauss points as a dictionary.
            
        GetStrain(): This function reads the plane strain data from the binary file
            and returns the strain data at Gauss points as a dictionary.
            
        GetModulus(): This function create the corresponding modulus data at Gauss points.

        PlGptStress(): This function plots the stress data at Gauss points.
        
        PlGptStrain(): This function plots the strain data at Gauss points.
        
        PlGptModulus(): This function plots the modulus data at Gauss points.

        ShapeFunAtGauss():This function is used to get the shape function value at Gauss Points.

        
    """
    
    # Initialize the mesh information
    def __init__(self, geoX, geoY, nelx, nely):
        """Initialize the mesh information."""
        self._geoX = geoX
        self._geoY = geoY
        self._nelx = nelx
        self._nely = nely
        self._nEl = nelx * nely
        self._NGauss = 4
        self._ex = geoX / nelx
        self._ey = geoY / nely
        self._nodsx = nelx + 1
        self._nodsy = nely + 1
        # Get the node coordinates
        self._X, self._Y = np.meshgrid(np.linspace(0,geoX,self.nodsx),np.linspace(0,geoY,self.nodsy))
        self._coord = np.vstack((self._X.ravel(), self._Y.ravel())).T
        # Get the nodes ID
        meshNods = np.int32(np.reshape(np.arange(1, (self._nodsx * self._nodsy) + 1, 1), (self._nodsy, self._nodsx)))
        self._eleNodsID = np.reshape(meshNods[:-1, :-1] + 1, (self._nEl, 1)) + np.array([-1, 0, self._nely + 1, self._nely])
        # Get the element node ID
        self._xEle = np.take(self._coord[:,0], self._eleNodsID - 1)
        self._yEle = np.take(self._coord[:,1], self._eleNodsID - 1)
        # Get the degree of freedom
        self._nDof = 2 * self.nodsx * self.nodsy
        # Get DOF of each element
        ## Define the first DOF of each element at leftdown
        cVec = np.reshape(2 * meshNods[:-1, :-1] + 1, (self._nEl, 1))
        ## Define other DOF of each element
        cMat = cVec + np.array([-2, -1, 0, 1, 2 * self._nely + 2, 2 * self._nely + 3, 2 * self._nely, 2 * self._nely + 1])
        self._cMat = cMat - 1

    @property
    def geoX(self):
        return self._geoX
    
    @property
    def geoY(self):
        return self._geoY
    
    @property
    def nelx(self):
        return self._nelx
    
    @property
    def nely(self):
        return self._nely
    
    @property
    def nEl(self):
        return self._nEl
    
    @property
    def ex(self):
        return self._ex
    
    @property
    def ey(self):
        return self._ey
    
    @property
    def nodsx(self):
        return self._nodsx
    
    @property
    def nodsy(self):
        return self._nodsy
    
    @property
    def coord(self):
        return self._coord
    
    @property
    def X(self):
        return self._X
    
    @property
    def Y(self):
        return self._Y
    
    @property
    def nGauss(self):
        return self._NGauss
    
    @property
    def eleNodsID(self):
        return self._eleNodsID
    
    @property
    def nDof(self):
        return self._nDof

    @property
    def xEle(self):
        return self._xEle

    @property
    def yEle(self):
        return self._yEle

    @property
    def cMat(self):
        return self._cMat

    @property
    def stress(self):
        return self._stress
    
    @property
    def strain(self):
        return self._strain
    
    @property
    def gaussE(self):
        return self._gaussE

    @property
    def nodeE(self):
        return self._nodeE

    @property
    def eleE(self):
        return self._eleE

    @property
    def odbE(self):
        return self._odbE

    def GetStress(self, stress_bin):
        """This function reads the plane stress data from the binary file
            and returns the stress data at Gauss points as a dictionary.
        Args:
            stress_bin (str): binary file name
        """    
        nEl = self.nEl
        NGauss = self.nGauss
        nelx = self.nelx
        nely = self.nely
        # Read the stress data from the binary file
        Stress_tot = np.fromfile(stress_bin,dtype=np.float64)
        # Get the stress data in the last frame
        split_num = nEl * NGauss
        frame_num = len(Stress_tot) / split_num / 4
        Stress_lf = np.array(np.split(np.split(Stress_tot, frame_num)[-1], split_num))
        # Seperate Components & Store
        S11, S22, S33, S12 = Stress_lf[:, 0], Stress_lf[:, 1], Stress_lf[:, 2], Stress_lf[:, 3]
        # Reshape Stress to Gauss Points Distribution
        ## Gauss Point clockwise from (-1,1)
        start_x = np.array([0, 1, 1, 0])
        start_y = np.array([1, 1, 0, 0])
        SGM11, SGM22, SGM33, SGM12 = [np.zeros((2 * nely, 2 * nelx)) for _ in range(4)]
        # Get Stress Field
        for i in range(NGauss):
            startx = start_x[i]
            starty = start_y[i]
            SGM11[starty::2, startx::2] = S11[i::4].reshape((nely,nelx))
            SGM22[starty::2, startx::2] = S22[i::4].reshape((nely,nelx))
            SGM33[starty::2, startx::2] = S33[i::4].reshape((nely,nelx))
            SGM12[starty::2, startx::2] = S12[i::4].reshape((nely,nelx))
        # Store Stress & Strain Field
        Stress = {key: Comp for key, Comp in zip(['S11', 'S22', 'S33', 'S12'], [SGM11, SGM22, SGM33, SGM12])}
        self._stress = Stress
        return Stress
    
    def GetStrain(self, strain_bin):
        """This function reads the plane strain data from the binary file
            and returns the strain data at Gauss points as a dictionary.
        Args:
            strain_bin (str): binary file name
            nEl (int): Number of elements
            NGauss (int): Number of Gauss points
            nelx (int): Number of elements in x direction
            nely (int): Number of elements in y direction
        """    
        nEl = self.nEl
        NGauss = self.nGauss
        nelx = self.nelx
        nely = self.nely
        # Read the strain data from the binary file
        Strain_tot = np.fromfile(strain_bin,dtype=np.float64)
        # Get the strain data in the last frame
        split_num = nEl * NGauss
        frame_num = len(Strain_tot) / split_num / 4
        Strain_lf = np.array(np.split(np.split(Strain_tot, frame_num)[-1], split_num))
        # Seperate Components & Store
        E11, E22, E33, E12 = Strain_lf[:, 0], Strain_lf[:, 1], Strain_lf[:, 2], Strain_lf[:, 3]
        # Reshape Strain to Gauss Points Distribution
        ## Gauss Point clockwise from (-1,1)
        start_x = np.array([0, 1, 1, 0])
        start_y = np.array([1, 1, 0, 0])
        EGM11, EGM22, EGM33, EGM12 = [np.zeros((2 * nely, 2 * nelx)) for _ in range(4)]
        # Get Strain Field
        for i in range(NGauss):
            startx = start_x[i]
            starty = start_y[i]
            EGM11[starty::2, startx::2] = E11[i::4].reshape((nely,nelx))
            EGM22[starty::2, startx::2] = E22[i::4].reshape((nely,nelx))
            EGM33[starty::2, startx::2] = E33[i::4].reshape((nely,nelx))
            EGM12[starty::2, startx::2] = E12[i::4].reshape((nely,nelx))
        # Store Strain Field
        Strain = {key: Comp for key, Comp in zip(['E11', 'E22', 'E33', 'E12'], [EGM11, EGM22, EGM33, EGM12])}
        self._strain = Strain
        return Strain

    def GetGaussE(self, alpha, beta, dis_type):
        """This function is used to get the modulus data at Gauss points.
        Args:
            alpha (float): Linear Parameter along X direction
            beta (float): Linear Parameter along Y direction
            dis_type (str): Type of distribution, '2D Linear' or '2D Exponential'
        """
        if dis_type == '2D Linear':
            EleModulus = 1.0 + self.xEle * alpha + self.yEle * beta
        elif dis_type == '2D Exponential':
            EleModulus = 1.0 * np.exp(self.xEle * alpha + self.yEle * beta)
        else:
            raise ValueError('Invalid Distribution Type')
        Modulus = np.zeros((2 * self.nely, 2 * self.nelx))
        GaussModulus = np.dot(PlFGMInfo.ShapeValue(), EleModulus.T).reshape((-1), order='F')
        # Gauss Point clockwise from (-1,1)
        start_x = np.array([0, 1, 1, 0])
        start_y = np.array([1, 1, 0, 0])
        for i in range(self.nGauss):
            startx = start_x[i]
            starty = start_y[i]
            Modulus[starty::2, startx::2] = GaussModulus[i::4].reshape((self.nely, self.nelx))
        self._gaussE = Modulus
        return Modulus
    
    def GetOdbE(self, alpha, beta, dis_type):
        """This function is used to get the modulus data at ODB points.
        Args:
            alpha (float): Expotential Parameter along X direction
            beta (float): Expotential Parameter along Y direction
            dis_type (str): Type of distribution, '2D Linear' or '2D Exponential'
        """
        if dis_type == '2D Exponential':
            EleModulus = 1.0 * np.exp(self.X * alpha + self.Y * beta)
        elif dis_type == '2D Linear':
            EleModulus = 1.0 + self.X * alpha + self.Y * beta
        else:
            raise ValueError('Invalid Distribution Type')
        self._odbE = tuple(map(tuple,np.split(EleModulus.flatten(), (self.nodsx * self.nodsy))))
        return EleModulus

    def GetNodeE(self, alpha, beta, dis_type):
        """This function is used to get the modulus data at element nodes.
        Args:
            alpha (float): Linear Parameter along X direction
            beta (float): Linear Parameter along Y direction
            dis_type (str): Type of distribution, '2D Linear' or '2D Exponential'
        """
        if dis_type == '2D Linear':
            EleModulus = 1.0 + self.X * alpha + self.Y * beta
        elif dis_type == '2D Exponential':
            EleModulus = 1.0 * np.exp(self.X * alpha + self.Y * beta)
        else:
            raise ValueError('Invalid Distribution Type')
        self._nodeE = EleModulus
        return EleModulus

    def GetEleE(self, alpha, beta, dis_type):
        """This function is used to get the modulus data at element nodes of each element.
        Args:
            alpha (float): Linear Parameter along X direction
            beta (float): Linear Parameter along Y direction
            dis_type (str): Type of distribution, '2D Linear' or '2D Exponential'
        """
        if dis_type == '2D Linear':
            Eele = 1.0 + self.xEle * alpha + self.yEle * beta
        elif dis_type == '2D Exponential':
            Eele = 1.0 * np.exp(self.xEle * alpha + self.yEle * beta)
        else:
            raise ValueError('Invalid Distribution Type')
        self._eleE = Eele
        return Eele

    def GetGauss(self):
        """This function is used to get the Gauss points coordinates."""
        gptshape, _, _ = PlFGMInfo.ShapeFunAtGauss()       
        Gauss_X = np.dot(gptshape, self.xEle.T).reshape((2 * self.nely, 2 * self.nelx), order='F')
        Gauss_Y = np.dot(gptshape, self.yEle.T).reshape((2 * self.nely, 2 * self.nelx), order='F')
        return Gauss_X, Gauss_Y

    def PltGptStress(self):
        """This function is used to plot the stress data at Gauss points."""
        X = np.linspace(0, self.geoX, 2 * self.nelx)
        Y = np.linspace(0, self.geoY, 2 * self.nely)
        fig, axes = plt.subplots(1, 4, figsize=(16,4))
        i = 0
        for key, comp in self.stress.items():
            pc = axes[i].pcolormesh(X, Y, comp, cmap='RdYlBu')
            fig.colorbar(pc, ax=axes[i])
            axes[i].set_title(key + "  Distribution")
            axes[i].set_xlabel("X")
            axes[i].set_ylabel("Y")
            axes[i].tick_params(axis='both', which='both', length=0)
            i += 1
        plt.tight_layout()
        plt.show()

    def PltGptStrain(self):
        """This function is used to plot the strain data at Gauss points."""
        X = np.linspace(0, self.geoX, 2 * self.nelx)
        Y = np.linspace(0, self.geoY, 2 * self.nely)
        fig, axes = plt.subplots(1, 4, figsize=(16, 4))
        i = 0
        for key, comp in self.strain.items():
            pc = axes[i].pcolormesh(X, Y, comp, cmap='RdYlBu')
            fig.colorbar(pc, ax=axes[i])
            axes[i].set_title(key + "  Distribution")
            axes[i].set_xlabel("X")
            axes[i].set_ylabel("Y")
            axes[i].tick_params(axis='both', which='both', length=0)
            i += 1
        plt.tight_layout()
        plt.show()
    
    def PltGptModulus(self):
        """This function is used to plot the modulus data at Gauss points."""
        X = np.linspace(0, self.geoX, 2 * self.nelx)
        Y = np.linspace(0, self.geoY, 2 * self.nely)
        fig, axes = plt.subplots()
        pc = axes.pcolormesh(X, Y, self.gaussE, cmap='RdYlBu')
        fig.colorbar(pc, ax=axes)
        axes.set_title("Modulus Distribution")
        axes.set_xlabel("X")
        axes.set_ylabel("Y")
        axes.tick_params(axis='both', which='both', length=0)
        plt.show()
    
    def PltNodeModulus(self):
        """This function is used to plot the modulus data at element nodes."""
        X = np.linspace(0, self.geoX, self.nodsx)
        Y = np.linspace(0, self.geoY, self.nodsy)
        fig, axes = plt.subplots()
        pc = axes.pcolormesh(X, Y, self.nodeE, cmap='RdYlBu')
        fig.colorbar(pc, ax=axes)
        axes.set_title("Modulus Distribution")
        axes.set_xlabel("X")
        axes.set_ylabel("Y")
        axes.tick_params(axis='both', which='both', length=0)
        plt.show()

    def GetFGMBD(self, alpha, beta, dis_type):
        """This function is used to get the B matrix & D matrix 
        of FGM throughout all elements.
        Args:
            alpha (float): The first parameter of the FGM.
            beta (float): The second parameter of the FGM.
            dis_type (str): Type of distribution, '2D Linear' or '2D Exponential'
        
        """
        mu = 0.3
        N, dNs, dNt, w = PlFGMInfo.ShapeFunAtGauss()
        Eele = self.GetEleE(alpha, beta, dis_type)
        xEle, yEle = self.xEle, self.yEle
        nIntPoint, nShape = N.shape
        Bgauss = []
        detJgauss = []
        Dgauss = []
        D0 = np.array([[1, mu, 0], [mu, 1, 0], [0, 0, 0.5 * (1 - mu)]]) / (1 - mu ** 2)
        for i in range(nIntPoint):
            # Get the shape function value at each Gauss point
            N_i = N[i, :]
            dNs_i = dNs[i, :]
            dNt_i = dNt[i, :]
            # Get the Jacobian matrix at each Gauss point
            # [1, nEl]
            dxds = np.matmul(dNs_i, xEle.T).reshape((1, -1))
            dxdt = np.matmul(dNt_i, xEle.T).reshape((1, -1))
            dyds = np.matmul(dNs_i, yEle.T).reshape((1, -1))
            dydt = np.matmul(dNt_i, yEle.T).reshape((1, -1))
            detJ = dxds * dydt - dyds * dxdt
            dxds = dxds / detJ
            dxdt = dxdt / detJ
            dyds = dyds / detJ
            dydt = dydt / detJ
            # Get grad(u) & grad(v) at each Gauss point in s-t coord
            duds = np.zeros((1, nShape * 2))
            dudt = np.zeros((1, nShape * 2))
            dvds = np.zeros((1, nShape * 2))
            dvdt = np.zeros((1, nShape * 2))
            duds[0, 0::2] = dNs_i
            dudt[0, 0::2] = dNt_i
            dvds[0, 1::2] = dNs_i
            dvdt[0, 1::2] = dNt_i
            # Get grad(u) & grad(v) at each Gauss point in x-y coord
            # [nEl, 8]
            dudx = np.matmul(dydt.T, duds) - np.matmul(dyds.T, dudt)
            dudy = - np.matmul(dxdt.T, duds) + np.matmul(dxds.T, dudt)
            dvdx = np.matmul(dydt.T, dvds) - np.matmul(dyds.T, dvdt)
            dvdy = - np.matmul(dxdt.T, dvds) + np.matmul(dxds.T, dvdt)
            # Get strain
            ex = dudx
            ey = dvdy
            exy = dudy + dvdx
            # Get B at each Gauss point [nEl, 3, 8]
            Bgauss_i = np.hstack((ex, ey, exy))
            Bgauss.append(Bgauss_i.reshape(-1, 3, 8))
            detJgauss.append(detJ.T)
            # Get D at each Gauss point
            Egauss_i = np.matmul(N_i, Eele.T).reshape((-1, 1))
            D = np.kron(Egauss_i, D0)
            Dgauss.append(np.vsplit(D, self.nEl))
        return Bgauss, detJgauss, Dgauss, w

    def GetFMGKe(self, Bgauss, detJgauss, Dgauss, w):
        """This function is used to get the element stiffness matrix."""
        nEl = self.nEl
        nGauss = len(w)
        # Initiate ke
        ke = [np.zeros((8, 8)) for _ in range(nEl)]
        for i in range(nGauss):
            # Get B, detJ, D at each Gauss point
            B = Bgauss[i]
            detJ = detJgauss[i]
            D = Dgauss[i]
            # Calculate the element stiffness matrix
            ke = [ke[j] + w[i] * np.matmul(np.matmul(B[j].T, D[j]), B[j]) \
                     * detJ[j] for j in range(nEl)]
        keVec = np.concatenate([Ke[np.triu(np.ones(Ke.shape, dtype=bool))] for Ke in ke])
        np.set_printoptions(threshold=np.inf)
        return ke, keVec

    def GetFGMKo(self, alpha, beta, dis_type):
        Bgauss, detJgauss, Dgauss, w = self.GetFGMBD(alpha, beta, dis_type)
        ke, keVec = self.GetFMGKe(Bgauss, detJgauss, Dgauss, w)
        # Generate degree list of global stiffness matrix
        si = []
        sii = []
        for i in range(1,9):
            si.extend(range(i,9))
            sii.extend([i] * (8 - i + 1))
        si = np.array(si)
        sii = np.array(sii)
        ik = self.cMat[:, si - 1].T
        jk = self.cMat[:, sii - 1].T
        ik_vec = ik.reshape((-1,1), order='F')
        jk_vec = jk.reshape((-1,1), order='F')
        Ivar = abs(np.sort(-np.hstack((ik_vec, jk_vec)), axis=1))
        K = sparse.coo_matrix((keVec, (Ivar[:,0], Ivar[:,1])), shape=(self.nDof, self.nDof))
        return K.tocsc()

    def GetForce(self, P):
        """This function is used to get the load vector.
        Args:
            P (float): Pressure Load
        Returns:
            F (sparse.csr): Load vector
        """
        F_Mag = P * self.ey
        loadPoint = np.argwhere(self.coord[:, 0] == self.geoX).reshape(-1)
        loadDof = loadPoint * 2
        F_Vec = np.ones_like(loadDof) * F_Mag
        F_Vec[0] -= F_Mag / 2
        F_Vec[-1] -= F_Mag / 2
        F = sparse.coo_matrix((F_Vec, (loadDof, np.zeros(len(loadDof)))), shape=(self.nDof, 1))
        return F.tocsr()

    @staticmethod
    def ShapeFunAtGauss():
        """This function is used to get the shape function 
            and its derivatives at Gauss points.
        Returns:
            N (np.array): Shape function at Gauss points
            dNs (np.array): Shape function derivatives in x direction at Gauss points
            dNt (np.array): Shape function derivatives in y direction at Gauss points
        """    
        gpt = 1 / np.sqrt(3)
        # Gauss Point clockwise from (-1,1)
        s = gpt * np.array([-1, 1, 1, -1]).reshape((4,1))
        t = gpt * np.array([1, 1, -1, -1]).reshape((4,1))
        w = np.array([1, 1, 1, 1]).reshape((4,1))
        # Node Counterclockwise From LeftDown
        N1 = (1.0-s)*(1.0-t) / 4.0
        N2 = (1.0+s)*(1.0-t) / 4.0
        N3 = (1.0+s)*(1.0+t) / 4.0
        N4 = (1.0-s)*(1.0+t) / 4.0
        N1s = -(1.0 - t) / 4.0
        N1t = -(1.0 - s) / 4.0
        N2s = (1.0 - t) / 4.0
        N2t = -(1.0 + s) / 4.0
        N3s = (1.0 + t) / 4.0
        N3t = (1.0 + s) / 4.0
        N4s = -(1.0 + t) / 4.0
        N4t = (1.0 - s) / 4.0
        N = np.hstack((N1, N2, N3, N4))
        dNs = np.hstack((N1s, N2s, N3s, N4s))
        dNt = np.hstack((N1t, N2t, N3t, N4t))
        return N, dNs, dNt, w

def BatchData(geoX, geoY, nelx, nely, Batch_Dir,Batch_Num, alpha, beta, dis_type):
    """This function is used to get data from all batches.

    Args:
        geoX (float): length of the x direction
        geoY (float): length of the y direction
        nelx (int): _number of elements in the x direction
        nely (int): number of elements in the y direction
        Batch_Dir (str): filepath of batch documents
        Batch_Num (int): number of batches
        alpha (np.1Darray): alpha values of all batches
        beta (np.1Darray): beta values of all batches
        dis_type (str): Type of distribution, '2D Linear' or '2D Exponential'

    Returns:
        Batch_dict: dictionary containing data from all batches
        Stress: np.4Darray containing the stress data from all batches
        Strain: np.4Darray containing the strain data from all batches
        Modulus: np.3Darray containing the modulus data from all batches
    """  
    # Initiate the dictionary, np.4Darray
    # Ignore S33 & E33
    Batch_dict = {} 
    Stress = np.zeros((Batch_Num, 2 * nely, 2 * nelx, 3))
    Strain = np.zeros((Batch_Num, 2 * nely, 2 * nelx, 3))
    Modulus = np.zeros((Batch_Num, nelx + 1, nely + 1))
    for i in range(Batch_Num):
        abqWorkDir = os.path.join(Batch_Dir, str(i + 1))
        os.chdir(abqWorkDir)
        stress_bin = os.path.join(abqWorkDir, str(i + 1) + 'Stress.bin')
        strain_bin = os.path.join(abqWorkDir, str(i + 1) + 'Strain.bin')
        Sample = PlFGMInfo(geoX, geoY, nelx, nely)
        Stress_tmp, Strain_tmp, Modulus[i, :, :] = Sample.GetStress(stress_bin), \
            Sample.GetStrain(strain_bin), Sample.GetNodeE(alpha[i], beta[i], dis_type)   
        Stress[i,:,:,0] = Stress_tmp['S11']
        Stress[i,:,:,1] = Stress_tmp['S22']
        Stress[i,:,:,2] = Stress_tmp['S12']
        Strain[i,:,:,0] = Strain_tmp['E11']
        Strain[i,:,:,1] = Strain_tmp['E22']
        Strain[i,:,:,2] = Strain_tmp['E12']
        Batch_dict[str(i + 1)] = Sample
        print('Batch ' + str(i + 1) + ' is done!')
    return Batch_dict, Stress, Strain, Modulus

def BatchMat(geoX, geoY, nelx, nely, Batch_Dir, Batch_Num, alpha, beta, dis_type):
    """This function is used to get material field & stiffness matrix for each batch.
    
    Args:
        geoX (float): length of the x direction
        geoY (float): length of the y direction
        nelx (int): number of elements in the x direction
        nely (int): number of elements in the y direction
        Batch_Dir (str): filepath of batch documents
        Batch_Num (int): number of batches
        alpha (np.1Darray): alpha values of all batches
        beta (np.1Darray): beta values of all batches
        dis_type (str): Type of distribution, '2D Linear' or '2D Exponential'

    """
    for i in range(Batch_Num):
        Sample = PlFGMInfo(geoX, geoY, nelx, nely)
        Modulus, Stiff = Sample.GetNodeE(alpha[i], beta[i], dis_type), \
            Sample.GetFGMKo(alpha[i], beta[i], dis_type)
        abqWorkDir = os.path.join(Batch_Dir, str(i + 1))
        os.chdir(abqWorkDir)
        # Save in Each Batch File
        np.save(str(i + 1) + 'Modulus', Modulus)
        sparse.save_npz(str(i + 1) + 'Stiff', Stiff, compressed=True)
        print('Batch ' + str(i + 1) + ' is done!')

def BatchForce(geoX, geoY, nelx, nely, Batch_Dir,Batch_Num, P):
    """This function is used to get the force vector for each batch.
    
    Args:
        geoX (float): length of the x direction
        geoY (float): length of the y direction
        nelx (int): number of elements in the x direction
        nely (int): number of elements in the y direction
        Batch_Dir (str): filepath of batch documents
        Batch_Num (int): number of batches
        P (float): pressure load
    """
    for i in range(Batch_Num):
        Sample = PlFGMInfo(geoX, geoY, nelx, nely)
        Force = Sample.GetForce(P)
        abqWorkDir = os.path.join(Batch_Dir, str(i + 1))
        os.chdir(abqWorkDir)
        # Save in Each Batch File
        sparse.save_npz(str(i + 1) + 'Force', Force, compressed=True)
        print('Batch ' + str(i + 1) + ' is done!')

def DictToBin(tar_dict, filename):
    """This function is used to write the data in the dictionary to a binary file.
    
    Args:
        tar_dict (dict): data to be written
        filename (str): binary file name
            
    """

    binary_data = pickle.dumps(tar_dict)
    with open(filename, 'wb') as f:
        f.write(binary_data)

def BinToDict(filename):
    """This function is used to read the binary file and return the data as a dictionary.

    Args:
        filename (str): binary file name

    Returns:
        dict: data in the binary file
    """        
    with open(filename, 'rb') as f:
        binary_data = f.read()
    
    tar_dict = pickle.loads(binary_data)
    return tar_dict


def CreateDir(file):
    """This function is used to create a new directory.

    Args:
        file (str): absolute path of the directory.
    """    
    if not os.path.exists(file):
        os.makedirs(file)
    else:
        print('The directory already exists.')