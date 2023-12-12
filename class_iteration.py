# class_iteration.py

import numpy as np

#from class_RDF import RDF

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
        """
        Initialize a System object.

        Parameters
        ----------
        groups : list
                List of group instances in the system.

        Attributes
        ----------
        groups : list
                List of group instances in the system.
        n_atoms : int
                Number of atoms in the system.
        n_type : int
                Number of different atom types in the system.
        ax : numpy.ndarray
                Vector representing the 'ax' lattice vector.
        ay : numpy.ndarray
                Vector representing the 'ay' lattice vector.
        az : numpy.ndarray
                Vector representing the 'az' lattice vector.
        U_pot : float
                Potential energy of the system.
        dt : float
                Time interval.
        N_iteration : int
                Number of iterations.
        alat_to_angstrom : float
                Conversion factor from lattice constant to angstrom.
        Ryau_to_pN : float
                Conversion factor from Rydberg atomic units to piconewtons.

        Notes
        -----
        This constructor sets up a 'MolecularDynamicsStep' object with the specified attributes, including the list of group instances in the system, the number of atoms, the number of different atom types, lattice vectors ('ax', 'ay', 'az'), time interval, and conversion factors.
        The elements of the list 'group' will be updated at each cycle, together with the potential energy and the number of iterations.

        Examples
        --------
        >>> group_instance_1 = group("example_type_1", 1)
        >>> group_instance_2 = group("example_type_2", 2)
        >>> iteration_instance = iteration([group_instance_1, group_instance_2])
        >>> print(iteration_instance.groups)
        [group_instance_1, group_instance_2]
        >>> print(iteration_instance.n_atoms)
        0
        >>> print(iteration_instance.n_type)
        0
        >>> print(iteration_instance.ax)
        [0. 0. 0.]
        >>> print(iteration_instance.ay)
        [0. 0. 0.]
        >>> print(iteration_instance.az)
        [0. 0. 0.]
        >>> print(iteration_instance.U_pot)
        0
        >>> print(iteration_instance.dt)
        0.0
        >>> print(iteration_instance.N_iteration)
        0
        >>> print(iteration_instance.alat_to_angstrom)
        0.0
        >>> print(iteration_instance.Ryau_to_pN)
        41778.2489644
        """
        self.groups = groups
        self.n_atoms = 0
        self.n_type = 0
        self.ax = np.zeros(3)
        self.ay = np.zeros(3)
        self.az = np.zeros(3)
        self.U_pot = 0
        self.dt = 0.0
        self.N_iteration = 0
        self.alat_to_angstrom = 0.0
        self.Ryau_to_pN = 4.17782489644e4

  
    def convert_alat_to_angstrom(self, celldim):
        """
        Convert lattice constant from Alat unit to Angstrom.

        Parameters
        ----------
        celldim : float
            Lattice constant in atomic units (a.u.).

        Notes
        -----
        This method sets the conversion factor from lattice constant to angstrom.

        Examples
        --------
        >>> iteration_istance = iteration(groups)
        >>> iteration_istance.convert_alat_to_angstrom(5.0)
        >>> print(iteration_istance.alat_to_angstrom)
        2.645886245
        """
        self.alat_to_angstrom = celldim * 0.529177249


    def set_mass(self, _type, _mass):
        """
        Set the mass of atoms in the specified group.

        Parameters
        ----------
        _type : str
            Name of the atom.
        _mass : float
            Mass of the atom.

        Notes
        -----
        This method iterates over the groups in the system and adds atoms with the specified name and mass.

        Examples
        --------
        >>> iteration_instance = iteration(groups)
        >>> iteration_instance.set_mass("C", 12.0)
        >>> iteration_instance.set_mass("O", 16.0)
        """
        for gr in self.groups:
            if _type in gr.type:
                gr.Add_atom(_type, _mass)


    def count_group(self, line):
        """
        Count the number of atoms for each group and extract the initial position.

        Parameters
        ----------
        line : list
            The read line of the file.

        Notes
        -----
        This method iterates over the groups in the system, counts the number of atoms for each group, and extracts the initial position of the atoms.

        Examples
        --------
        >>> iteration_instance = iteration(groups)
        >>> line_data = ["1", "C", ...]  # Example line data from the file
        >>> iteration_instance.count_group(line_data)
        """
        atom_type = line[1] 
        # Check on the groups
        for gr in self.groups:
            if atom_type in gr.type:
                gr.id_tot = np.append(gr.id_tot, int(line[0]))
                # Check on the group's elements
                for at in gr.atoms:
                    if atom_type == at.name:
                        at.N += 1
                        at.id = np.append(at.id, int(line[0]))
                        at.add_position_past(float(line[6]) * self.alat_to_angstrom, float(line[7]) * self.alat_to_angstrom, float(line[8]) * self.alat_to_angstrom)
                        break
                break
    

    def set_DOF(self):
        """
        Set the degrees of freedom for each group in the system.

        Notes
        -----
        This method iterates over the groups in the system, calculates the degrees of freedom based on the total number of atoms in the group, and adjusts it if the group's `id_group` is not equal to 0.

        Examples
        --------
        >>> iteration_instance = iteration(groups)
        >>> iteration_instance.set_DOF()
        """
        for gr in self.groups:
            gr.DOF = 3 * len(gr.id_tot)
            if not gr.id_group == 0:
                gr.DOF = gr.DOF - 3


    def forces(self, line):
        """
        Extract the force acting on atoms at a specific time.

        Parameters
        ----------
        line : list
            The line containing the force acting on the atom.

        Notes
        -----
        This method iterates over the groups in the system, finds the corresponding atom, and adds the force to the atom and group.

        Examples
        --------
        >>> iteration_instance = iteration(groups)
        >>> line = [0, 1, ..., 2.0, 3.0, 4.0, ...]
        >>> iteration_instance.forces(line)
        """
        atom_type = int(line[1])
        
        # Check on the groups
        for gr in self.groups:
            if atom_type in gr.id_tot:
                gr.Add_force(float(line[6]) * self.Ryau_to_pN, float(line[7]) * self.Ryau_to_pN, float(line[8]) * self.Ryau_to_pN)
                # Check on the group's elements
                for at in gr.atoms:
                    if atom_type in at.id:
                        at.add_force(float(line[6]) * self.Ryau_to_pN, float(line[7]) * self.Ryau_to_pN, float(line[8]) * self.Ryau_to_pN)
                        break
                break

    def positions(self, line, istogram):
        """
        Extract and store the positions of atoms at a determined time.

        Parameters
        ----------
        line : list
            The line containing the atom positions.
        istogram : RDF
            The RDF instance with which the radial distribution function is calculated.

        Notes
        -----
        This method iterates over the groups in the system, finds the corresponding atom, and adds the position to the atom.
        Additionally, it checks the atom type and adds positions to the RDF istance 'istogram' (if applicable).

        Examples
        --------
        >>> iteration_instance = iteration(['C', 'H', 'P', 'K'])
        >>> line = [0, 1.0, 2.0, 3.0, ...]
        >>> istogram = RDF("filename", 10, ['C', 'H'], 500)
        >>> iteration_instance.positions(line, istogram)
        """
        atom_type = line[0] 

        # Check if for the atom in line RDF will be calculated
        if atom_type == istogram.type[0]:
            istogram.add_position1(float(line[1]), float(line[2]), float(line[3]))
        if atom_type == istogram.type[1]:
            istogram.add_position2(float(line[1]), float(line[2]), float(line[3]))
        
        # Check on the groups
        for gr in self.groups:
            if atom_type in gr.type:
                # Check on the group's elements
                for at in gr.atoms:
                    if atom_type == at.name:
                        at.add_position(float(line[1]), float(line[2]), float(line[3]))
                        break
                break
        
    def single_frame(self):
        """
        Generate a list where each element is a line in the output file, representing the current time step.

        Returns
        -------
        text : list
            The list of lines.

        Notes
        -----
        This method prints information about the lattice, time step, number of iterations, and potential energy in the second line of the xyz file.
        It iterates over the groups in the system, adding information about kinetic energy, degrees of freedom, temperature, and total force.
        At the end, it adds the coordinates, velocity, and forces of each atom using the 'generate' method of the group.

        Examples
        --------
        >>> iteration_instance = iteration(groups)
        >>> lines = iteration_instance.single_frame()
        >>> print(lines)
        ['n_atoms Lattice(Ang)="...", ' dt(ps)=..., N=..., Epot(eV)=...']
        ['...']
        ['...']
        ['...']
        """    
        # Generate the general information of the time step
        text = [f'{self.n_atoms}', f'Lattice(Ang)=\"{self.ax[0]:.3f}, {self.ax[1]:.3f}, {self.ax[2]:.3f}, {self.ay[0]:.3f}, {self.ay[1]:.3f}, {self.ay[2]:.3f}, {self.az[0]:.3f}, {self.az[1]:.3f}, {self.az[2]:.3f}\" dt(ps)={self.dt:.6f} N={self.N_iteration} Epot(eV)={(self.U_pot / 13.60570398):.3f}']

        # Generate the information for each group  of the time step
        for gr in self.groups:            
            body = gr.Generate(self.dt)
            text[1] = text[1] + f' Ek{gr.id_group}(ev)={gr.Ek:.3f} DOF{gr.id_group}={gr.DOF} T{gr.id_group}(K)={gr.T:.3f} Ftot{gr.id_group}(pN)=\"{gr.Ftot[0]:.3f}, {gr.Ftot[1]:.3f}, {gr.Ftot[2]:.3f}\"'
            text.extend(body)

        return text