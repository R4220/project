# class_group
import numpy as np

class group:
    '''
    This class represents a group of atoms used to divide the whole atoms in a molecular dynamics simulation.
    In particular, the implemented methods calculate several properties of the group.

    Attributes:
     - type: atomic type in the group (str)
     - mass: mass of the element in the group (float)
     - id_group: identification number of the group (int, default = 0)
     - id: list of the identification number of the atoms in the group ([int])
     - N: number of atoms in a group
     - force: bidimensional array representing the forces against the atoms in the group 
     - position_past: bidimensional array representing the position of the atoms in the group during the last time step
     - position: bidimensional array representing the position of the atoms in the group during the current time step
     - velocity: bidimensional array representing the velocity of the atoms in the group
     - Ek: kinetic energy of the group
     - Ftot: total force acting on the group
     - T: temperature of the group
     - dof: degree of freedom of the group
    '''

    def __init__(self, type, id_group, move = True):
        '''
        Initialize a Group instance.

        Parameters:
         - type (str): Atomic type in the group.
         - mass (float): Mass of the element in the group.
        '''
        self.atoms = [] #
        self.type = type #
        self.move = move #
        self.id_group = id_group #
        self.id_tot = []#
        #self.atom_dict = {} 
        self.force = np.array([], dtype=float).reshape(0, 3)

    def add_atom(self, name, mass):
        self.atoms = np.append(self.atoms, self.atom(name, mass))

    def get_atom_by_id(self, atom_id):
        return self.atom_dict.get(atom_id, None) 
    
    def add_position(self, x, y, z):
        """
        Add the coordinates of the considered atom from the current time step to the list position_past.

        Parameters:
        - x (float): X-coordinate of the atom.
        - y (float): Y-coordinate of the atom.
        - z (float): Z-coordinate of the atom.
        """
        self.position = np.vstack([self.position, np.array([x, y, z], dtype=float)])
    
    def add_force(self, x, y, z):
        """
        Add the coordinates of the force acting on the considered atom from the current time step to the force list.

        Parameters:
        - x (float): X-coordinate of the force.
        - y (float): Y-coordinate of the force.
        - z (float): Z-coordinate of the force.
        """
        self.force = np.vstack([self.force, np.array([x, y, z], dtype=float)])

    class atom:
        def __init__(self, name, mass):
            self.name = name
            self.mass = mass
            self.N = 0
            self.id = []

            self.position_past = np.array([], dtype=float).reshape(0, 3)
            self.position = np.array([], dtype=float).reshape(0, 3)
            self.velocity = np.array([], dtype=float).reshape(0, 3)
        
        def add_id(self, i): 
            """
            Add the number 'i' to the list 'id', representing the identification number of the atom from the PWO file.

            Parameters:
            - i (int): Input integer number.
            """
            self.id = np.append(self.id, i)
            #self.id_tot = np.append(self.id_tot, i)

        def add_position_past(self, x, y, z):
            """
            Add the coordinates of the considered atom from the starting configuration to the list position_past.

            Parameters:
            - x (float): X-coordinate of the atom.
            - y (float): Y-coordinate of the atom.
            - z (float): Z-coordinate of the atom.
            """
            self.position_past = np.vstack([self.position_past, np.array([x, y, z], dtype=float)])
        
        def add_position(self, x, y, z):
            """
            Add the coordinates of the considered atom from the starting configuration to the list position_past.

            Parameters:
            - x (float): X-coordinate of the atom.
            - y (float): Y-coordinate of the atom.
            - z (float): Z-coordinate of the atom.
            """
            self.position = np.vstack([self.position, np.array([x, y, z], dtype=float)])

