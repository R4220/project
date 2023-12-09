# class_atom
import numpy as np

class atom:
    def __init__(self, name, mass, move, id_group):
        self.name = name
        self.mass = mass
        self.move = move
        self.id_group = id_group
        self.DOF = 0
        self.N = 0
        self.id = []

        self.position_past = np.array([], dtype=float).reshape(0, 3)
        self.position = np.array([], dtype=float).reshape(0, 3)
        self.velocity = np.array([], dtype=float).reshape(0, 3)
        self.force = np.array([], dtype=float).reshape(0, 3)
        
    def add_id(self, i): 
        """
        Add the number 'i' to the list 'id', representing the identification number of the atom from the PWO file.

        Parameters:
        - i (int): Input integer number.
        """
        self.id = np.append(self.id, i)
        self.N += 1

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
        
    def add_force(self, x, y, z):
        """
        Add the coordinates of the considered atom from the starting configuration to the list position_past.

        Parameters:
        - x (float): X-coordinate of the atom.
        - y (float): Y-coordinate of the atom.
        - z (float): Z-coordinate of the atom.
        """
        self.force = np.vstack([self.force, np.array([x, y, z], dtype=float)])

    def generarte_velocity(self, dt):
        """
        Generate the bidimensional array representing the velocity of each atom. Velocities are calculated as the difference between the current position and the position from the previous time step over the time interval.

        Parameters:
        - Dt (float): Time interval.
        """
        self.velocity = (np.array(self.position, dtype=float) - np.array(self.position_past, dtype=float)) / dt
