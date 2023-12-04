# class_iteration.py

import numpy as np

from class_group import group

class iteration:
    '''
    This class represents a single step dureing a molecular dynamics

    Attributes:
     - n_at: total number of atoms in the simulation
     - celldim: dimension of the the cell (a.u)
     - ax, ay, az: lattice vectors 
     - L: cell dimension [Lx, Ly, Lz]
     - alat_to_angstrom: conversion between Alat unit and Angstrom
     - t_past: time istant of the last time step
     - t: time istant of the current time step
     - U_pot: potential energy
     - blockp: array in which the lines of the pwo file about the current file are stored
     - groups: list of the groups in the dynamic
     - RDF_atoms: atoms for which the RDF is calculated 
     - switch_at: bolean that turn in or off the acquisition of the atomic positions (default = False)
     - countp: array of the counts during the RDF calculation
     - R: array representing the possible distances considered in the RDF calculation
     - dR: bin length in the RDF istogram
     - norm: normalization values for the rDF
     - e: bolean value representing if the atoms for which the RDF is clculated are the same (True) or not (False)
    '''

    def __init__(self):
        self.n_at = 0
        self.celldim = 0
        self.ax = [0, 0, 0]
        self.ay = [0, 0, 0]
        self.az = [0, 0, 0]
        self.L = [0, 0, 0]
        self.alat_to_angstrom = 0.
        self.t_past = 0
        self.t = 0
        self.U_pot = 0
        self.block = np.array([], dtype=str)
        self.groups = []
        self.RDF_atoms = []
        self.switch_at = False 
        self.count = np.array([], dtype = int)
        self.R = np.array([], dtype=float)
        self.dR = 0.
        self.norm = np.array([], dtype=float)
        self.e = False
    
    def Alat_to_Angstrom(self):
        '''
        This method define the conversion parameter between the Alat unit of length (that depends on celldim) and angstrom
        '''
        self.alat_to_angstrom = self.celldim * 0.529177249

    def store(self, line):
        '''
        This method add to block the line in the current iteration, if it is not empty
        
        Parameters:
            line: the line readed
        '''
        if line != '':
            self.block = np.append(self.block, line)
    
    def add_group(self, _type, _mass, RDF):#, idx):
        '''
        This method add to groups the identification of the group
        
        Parameters:
            _type: name of the atom
            _mass: mass of the atom
        '''
        at = group(_type, _mass)
        self.groups = np.append(self.groups, at )
        if _type in RDF:
            self.RDF_atoms = np.append(self.RDF_atoms, at)
    
    def count_group(self, line):
        '''
        This method counts the number of atoms for each group
        
        Parameters:
            line: the readed line of the file
            '''
        for elm in self.groups:
            if elm.type == line[1]:
                elm.N += 1
                elm.add_id(int(line[0]))
                elm.add_position_past(float(line[6]) * self.alat_to_angstrom, float(line[7]) * self.alat_to_angstrom, float(line[8]) * self.alat_to_angstrom)
                break
    
    def forces(self, line):
        '''
        This method extract the force at a determined t
        
        Parameters:
            line: the line in which there are the coordinates of the force        
        '''
        _line = line.split()
        for elm in self.groups:
            if int(_line[1]) in elm.id:
                elm.add_force(_line[6], _line[7], _line[8])
                break
    
    def position(self, line):
        '''
        This method extract the position at a determined t
        
        Parameters:
            line: the line in which there are the coordinates of the position        
        '''
        _line = line.split()
        for elm in self.groups:
            if _line[0] == elm.type:
                elm.add_position(_line[1], _line[2], _line[3])
                break

    def text(self):
        '''
        This method generates a list in which each elements is a line in the output file, representing the current time step.
        
        Return: 
         - text: the list of the lines
        '''
        text = f'{self.n_at}\nLattice(Ang)=\"{self.ax[0]}, {self.ax[1]}, {self.ax[2]}, {self.ay[0]}, {self.ay[1]}, {self.ay[2]}, {self.az[0]}, {self.az[1]}, {self.az[2]}\" t(ps)={self.t} Epot(eV)={self.U_pot / 13.60570398 }'
        body = np.array([], dtype=str)
        for elm in self.groups:
            _body = elm.generate(self.t - self.t_past)
            elm.DOF = elm.N
            text = text + f" Ek{elm.id_group}(ev)={elm.Ek} DOF{elm.id_group}={elm.DOF} T{elm.id_group}(K)={elm.T} Ftot{elm.id_group}(pN)=\"{elm.Ftot[0]}, {elm.Ftot[1]}, {elm.Ftot[2]}\""
            body = np.append(body, _body)
        text = text.split('\n')
        text = np.append(text, body)
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
        This method calculates the distances and store them into the list dist, in order to calculate the RDF. 
        Two cases are present, defined through the bolean variable e: if e = True, we are considering RDF is calculated between the same group, while if e = False, two different group are considered.
        during the calculation.

        Return:
         - dist: array of distances
        '''
        dist = []
        if self.e == True:
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
            #print(pos_1)
            #print(pos_2)
            for i in pos_1:
                #print('i ', i, '   pos2 ', pos_2)
                _pos = i - pos_2
                #print(_pos)
                n = np.round(np.divide(_pos, self.L))
                #print(_pos - np.multiply(self.L, n))
                #print(np.linalg.norm(_pos - np.multiply(self.L, n)))
                r = np.linalg.norm(_pos - np.multiply(self.L, n))
                dist.extend(r[r < RDF[0]])
        dist = np.sort(dist)
        return dist

    def istogram(self, RDF):
        '''
        This method add to count the list given in input, in order to generate the RDF istogram
        
        Parameters:
         - RDF: list with several information neaded to calculate the RDF and define the size of the x-axis of the istogram
        '''
        self.count += np.histogram(self.RDF(RDF), bins=RDF[3], range=(0, RDF[0]))[0]

    def single_frame(self, RDF):
        '''
        In this method several things are performed:
         - the elements of block are checked in order to
           - extract the time of the current step
           - extract the potential energy
           - extract the coordinates of the forces against the atoms, calling the method forces
           - extract the coordinates of the positions of atoms, calling the method position
         - call the method istogram
         - after all the calculation the storing list block is setted again as empty
        
        Return:
         - text: the list in which the lines of the current time step for the output file are storred
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
        self.istogram(RDF)
        text = self.text()
        self.block = []
        return text
    
    def normalization(self):
        '''
        This method normalize count
        .. dire ancora come ..
        '''
        V = (self.L[0] * self.L[1] * self.L[2])**2
        #print(V)
        
        #print(self.RDF_atoms[0].N, self.RDF_atoms[1].N)
        if self.e == True:
            rho = self.RDF_atoms[0].N **2/ V
        else:
            rho = self.RDF_atoms[0].N * self.RDF_atoms[1].N / V
        #print(rho)
        #print(self.norm)
        self.norm = rho * self.norm
        self.count = np.divide(self.count, self.norm)
        '''qui non mi torna qualcosa'''
    