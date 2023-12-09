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
        self.groups = groups

        self.n_atoms = 0 #
        self.n_type = 0 #
        self.alat_to_angstrom = 0.0 #
        self.ax = np.zeros(3) #
        self.ay = np.zeros(3) #
        self.az = np.zeros(3) #
        
        self.U_pot = 0
        self.dt = 0. 
        self.N_iteration = 0
  
    def Alat_to_Angstrom(self, celldim):
        '''
        This method defines the conversion parameter between the Alat unit of length (that depends on celldim) and angstrom.
        '''
        self.alat_to_angstrom = celldim * 0.529177249

    def set_mass(self, _type, _mass):#, RDF)
        '''
        This method creates the atom given the name and the mass.

        Parameters:
        - _type: name of the atom
        - _mass: mass of the atom
        '''
        for gr in self.groups:
            if _type in gr.type:
               gr.Add_atom(_type, _mass) 

    def count_group(self, line):
        '''
        This method counts the number of atoms for each group and extract the initial position.

        Parameters:
        - line: the read line of the file
        '''
        atom_type = line[1] 
        for gr in self.groups:
            if atom_type in gr.type:
                gr.id_tot = np.append(gr.id_tot, int(line[0]))
                for at in gr.atoms:
                    if atom_type == at.name:
                        at.N += 1
                        at.id = np.append(at.id, int(line[0]))
                        at.add_position_past(float(line[6]) * self.alat_to_angstrom, float(line[7]) * self.alat_to_angstrom, float(line[8]) * self.alat_to_angstrom)
                        break
                break
    
    def set_DOF(self):
        for gr in self.groups:
            gr.DOF = 3 * len(gr.id_tot)
            if not gr.id_group == 0 :
                gr.DOF = gr.DOF - 3

    def forces(self, line):
        '''
        This method extracts the force at a determined time.
    
        Parameters:
        - line: the line in which there are the coordinates of the force
        '''
        
        atom_type = int(line[1])
                
        for gr in self.groups:
            if atom_type in gr.id_tot:
                for at in gr.atoms:
                    if atom_type in at.id:
                        at.add_force(float(line[6]), float(line[7]), float(line[8]))
                        gr.Add_force(float(line[6]), float(line[7]), float(line[8]))
                        break
                break

    def positions(self, line):
        atom_type = line[0] 
        
        for gr in self.groups:
            if atom_type in gr.type:
                for at in gr.atoms:
                    if atom_type == at.name:
                        at.add_position(float(line[1]), float(line[2]), float(line[3]))
                        break
                break
        
    
    def single_frame(self, RDF):
        '''
        This method generates a list in which each element is a line in the output file, representing the current time step.

        Return: 
        - text: the list of the lines
        '''        
        text = [f'{self.n_atoms}', f'Lattice(Ang)=\"{self.ax[0]:.3f}, {self.ax[1]:.3f}, {self.ax[2]:.3f}, {self.ay[0]:.3f}, {self.ay[1]:.3f}, {self.ay[2]:.3f}, {self.az[0]:.3f}, {self.az[1]:.3f}, {self.az[2]:.3f}\" dt(ps)={self.dt:.6f} N={self.N_iteration} Epot(eV)={(self.U_pot / 13.60570398):.3f}']

        for gr in self.groups:
            '''for at in gr.atoms:
                print(at.name, at.position_past)'''
            
            body = gr.Generate(self.dt)
            text[1] = text[1] + f' Ek{gr.id_group}(ev)={gr.Ek:.3f} DOF{gr.id_group}={gr.DOF} T{gr.id_group}(K)={gr.T:.3f} Ftot{gr.id_group}(pN)=\"{gr.Ftot[0]}, {gr.Ftot[1]}, {gr.Ftot[2]}\"'
            text.extend(body)

        return text

'''quando fai rdf e devi calcolare la matrice, dopo formatta i vettori "{:.4f}".format(n)'''