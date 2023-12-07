# class_iteration.py

import numpy as np

from class_group import group

class iteration:
    '''
    This class represents a single step during a molecular dynamics.

    Attributes:
    - n_atoms: total number of atoms in the simulation
    - cell_dim: dimension of the cell (a.u)
    - ax, ay, az: lattice vectors 
    - L: cell dimension [Lx, Ly, Lz]
    - alat_to_angstrom: conversion between Alat unit and Angstrom
    - t_past: time instant of the last time step
    - t: time instant of the current time step
    - U_pot: potential energy
    - block: array in which the lines of the pwo file about the current file are stored
    - groups: list of the groups in the dynamic
    - RDF_atoms: atoms for which the RDF is calculated 
    - switch_at: boolean that turns on or off the acquisition of the atomic positions (default = False)
    - count: array of the counts during the RDF calculation
    - R: array representing the possible distances considered in the RDF calculation
    - dR: bin length in the RDF histogram
    - norm: normalization values for the RDF
    - e: boolean value representing if the atoms for which the RDF is calculated are the same (True) or not (False)
    '''

    def __init__(self, groups):
        self.n_atoms = 0 #
        self.n_type = 0 #
        self.alat_to_angstrom = 0.0 #
        self.ax = np.zeros(3) #
        self.ay = np.zeros(3) #
        self.az = np.zeros(3) #
        self.groups = groups

        self.L = np.zeros(3)
        self.dt = 0. #
        self.N_iteration = 0 #
        self.t_past = 0
        self.t = 0
        self.U_pot = 0
        self.block = []
        
        self.RDF_atoms = []
        self.switch_at = False 
        self.count = np.zeros(0, dtype=int)
        self.R = np.zeros(0, dtype=float)
        self.dR = 0.0
        self.norm = np.zeros(0, dtype=float)
        self.e = False
    
    def Alat_to_Angstrom(self, celldim):
        '''
        This method defines the conversion parameter between the Alat unit of length (that depends on celldim) and angstrom.
        '''
        self.alat_to_angstrom = celldim * 0.529177249

    def set_mass(self, _type, _mass):#, RDF)
        '''
        This method adds the identification of the group to groups.

        Parameters:
        - _type: name of the atom
        - _mass: mass of the atom
        - RDF: list of atom types for RDF calculation
        '''
        for gr in self.groups:
            if _type in gr.type:
               gr.add_atom(_type, _mass) 
    
        # Check if the atom type is in the RDF list
        #if _type in RDF:
         #   self.RDF_atoms = np.append(self.RDF_atoms, at)
   
    def store(self, line):
        '''
        This method adds the line to block if it is not empty.

        Parameters:
        - line: the line read.
        '''
        if line != '':
            self.block = np.append(self.block, line)
     
    def count_group(self, line):
        '''
        This method counts the number of atoms for each group.

        Parameters:
        - line: the read line of the file
        '''
        atom_type = line[1] 

        for elm in self.groups:
            if atom_type in elm.type:
                for at in elm.atoms:
                    if atom_type == at.name:
                        at.N += 1
                        at.add_id(int(line[0]))
                        at.add_position_past(float(line[6]) * self.alat_to_angstrom, float(line[7]) * self.alat_to_angstrom, float(line[8]) * self.alat_to_angstrom)
                        break
                break

    def forces(self, line):
        '''
        This method extracts the force at a determined time.
    
        Parameters:
        - line: the line in which there are the coordinates of the force
        '''
    
        for elm in self.groups:
            if int(line[1])  in elm.id:
                elm.add_force(float(line[6]), float(line[7]), float(line[8]))
                break
 
    def position(self, line):
        '''
        This method extracts the position at a determined time.
    
        Parameters:
        - line: the line in which there are the coordinates of the position
        '''
        _line = line.split()
        atom_type = _line[0] 
    
        for elm in self.groups:
            if atom_type == elm.type:
                elm.add_position(float(_line[1]), float(_line[2]), float(_line[3]))
                break

    def text(self):
        '''
        This method generates a list in which each element is a line in the output file, representing the current time step.

        Return: 
        - text: the list of the lines
        '''
        text = f'{self.n_atoms}\nLattice(Ang)=\"{self.ax[0]}, {self.ax[1]}, {self.ax[2]}, {self.ay[0]}, {self.ay[1]}, {self.ay[2]}, {self.az[0]}, {self.az[1]}, {self.az[2]}\" t(ps)={self.t} Epot(eV)={self.U_pot / 13.60570398 }'
        body = []
    
        for elm in self.groups:
            _body = elm.generate(self.t - self.t_past)
            elm.DOF = elm.N
            text = text + f" Ek{elm.id_group}(ev)={elm.Ek} DOF{elm.id_group}={elm.DOF} T{elm.id_group}(K)={elm.T} Ftot{elm.id_group}(pN)=\"{elm.Ftot[0]}, {elm.Ftot[1]}, {elm.Ftot[2]}\""
            body.extend(_body)

        text = text.split('\n')
        text.extend(body)
        return text
 
    def set_RDF(self, RDF):
        '''
        This method defines the variable neaded for the RDF calculation ad hoc for the current simulation
        '''
        self.count = np.zeros(RDF[3])
        self.R = np.linspace(0, RDF[0], RDF[3])
        self.dR = self.R[1]
        self.norm = np.multiply([((i + self.dR)**3 - i**3) for i in self.R[:RDF[3]]], np.pi * 4 /3)
        self.e = (RDF[1] == RDF[2])

    def RDF(self, RDF):
        '''
        Calculate distances and store them in the list dist for RDF calculation. 
        Two cases are present, defined by the boolean variable e:
        - If e is True, RDF is calculated within the same group.
        - If e is False, RDF is calculated between two different groups.

        Parameters:
        - RDF: List containing parameters for RDF calculation.
     
        Return:
        - dist: Array of distances.
        '''
        dist = []

        if self.e:
            pos = self.RDF_atoms[0].position
            N = len(pos)

            for k in range(N - 1):
                _pos = pos[k+1:] - pos[k]
                n = np.round(np.divide(_pos, self.L))
                r = np.linalg.norm(_pos - np.multiply(self.L, n))
                dist.extend(r[r < RDF[0]])
        else:
            pos_1 = self.RDF_atoms[0].position
            pos_2 = self.RDF_atoms[1].position

            for i in pos_1:
                _pos = i - pos_2
                n = np.round(np.divide(_pos, self.L))
                r = np.linalg.norm(_pos - np.multiply(self.L, n))
                dist.extend(r[r < RDF[0]])

        dist = np.sort(dist)
        return dist

    def istogram(self, RDF):
        '''
        Add the distances calculated from RDF to the count array for generating the RDF histogram.

        Parameters:
        - RDF: List containing several pieces of information needed to calculate the RDF and define the size of the x-axis of the histogram.
        '''
        self.count += np.histogram(self.RDF(RDF), bins=RDF[3], range=(0, RDF[0]))[0]

    def single_frame(self, RDF):
        '''
        Perform several actions for a single time step:
        - Check the elements of the block to
            - Extract the time of the current step
            - Extract the potential energy
            - Extract the coordinates of the forces against the atoms by calling the method forces
            - Extract the coordinates of the positions of atoms by calling the method position
        - Call the method istogram
        - After all the calculations, set the storing list block as empty

        Return:
        - text: List containing the lines of the current time step for the output file.
        '''
        for line in self.block:
            if 'time      =' in line:
                self.t_past = self.t
                self.t = float(line.split()[2])
            if '!    total energy' in line:
                self.U_pot = float(line.split()[4])
            if '   force =    ' in line:
                self.forces(line)
            if 'ATOMIC_POSITIONS' in line:
                self.switch_at = True
            elif self.switch_at == True:
                _line = line.split()
                if _line == []:
                    self.switch_at = False
                else:
                    self.position(line)
        #self.istogram(RDF)
        text = self.text()
        self.block = []
        return text
    
    def normalization(self):
        '''
        Normalize the count variable.
        The normalization process involves dividing the count variable by a normalization factor.
        The normalization factor (rho) is calculated based on the volume of the simulation cell and the number of atoms involved in the RDF calculation. The specific calculation depends on whether the RDF is calculated between atoms of the same group (e=True) or different groups (e=False).
        - If e=True, rho = (Number of atoms in group 1)^2 / Volume
        - If e=False, rho = (Number of atoms in group 1 * Number of atoms in group 2) / Volume

        The normalized count variable is then obtained by dividing the original count by the normalization factor.
        '''
        V = (self.L[0] * self.L[1] * self.L[2]) ** 2

        if self.e == True:
            rho = self.RDF_atoms[0].N**2 / V
        else:
            rho = (self.RDF_atoms[0].N * self.RDF_atoms[1].N) / V

        self.norm = rho * self.norm
        self.count = np.divide(self.count, self.norm)
        ''' Check if the normalization process is correctly implemented '''

    