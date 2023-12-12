# class_group
import numpy as np

from class_atom import atom

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

    def __init__(self, type, id_group, move):
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
        self.id_tot = np.array([], dtype=int)
        self.DOF = 0.
        self.Ek = 0.

        self.Ftot = np.array([], dtype=float).reshape(0, 3)
        self.force = np.array([], dtype=float).reshape(0, 3)

        self.Vtot = np.array([], dtype=float).reshape(0, 3)
        self.velocity = np.array([], dtype=float).reshape(0, 3)

    def Add_atom(self, name, mass):
        self.atoms = np.append(self.atoms, atom(name, mass, self.move, self.id_group))
    print
    def Add_force(self, x, y, z):
        """
        Add the coordinates of the force acting on the considered atom from the current time step to the force list.

        Parameters:
        - x (float): X-coordinate of the force.
        - y (float): Y-coordinate of the force.
        - z (float): Z-coordinate of the force.
        """
        self.force = np.vstack([self.force, np.array([x, y, z], dtype=float)])

    def Kinetic_energy(self, dt):
        # Generate the velocity array for each atom type
        for at in self.atoms:
            at.generarte_velocity(dt)
            self.velocity = np.append(self.velocity, at.velocity)
        
        self.Vtot = np.sum(self.velocity, axis=0)/len(self.id_tot)
        #print(self.Vtot)
        for at in self.atoms:
            self.Ek += 0.5 * float(at.mass) * np.sum(np.linalg.norm(at.velocity - self.Vtot, axis=1)**2) * 0.0001036426948415943
        
    def Generate(self, dt):
        # Calculate the kinetic energy
        '''for at in self.atoms:
                print(at.name, at.position_past)'''
        self.Kinetic_energy(dt)

        # Calculate temperature
        self.T = (2 * self.Ek) / (self.DOF * 8.617333262145e-5)

        # Calculate total force
        self.Ftot = np.sum(self.force, axis=0) # le unità delle forze sono sbaghliarte

        # Generate output line
        body = np.array([], dtype=str)
        for at in self.atoms:
            for i in range(at.N):
                #print(at.name)
                #print('pos', at.position)
                #print('vel', at.velocity)
                #print('for', at.force)
                body = np.append(body, f'{at.name}\t  {at.position[i][0]}\t  {at.position[i][1]}\t  {at.position[i][2]}\t  {at.velocity[i][0]}\t  {at.velocity[i][1]}\t  {at.velocity[i][2]}\t  {at.force[i][0]}\t  {at.force[i][1]}\t  {at.force[i][2]}\t  {at.id_group}')
            at.position_past = at.position
            at.position = np.array([], dtype=float).reshape(0, 3)
            at.force = np.array([], dtype=float).reshape(0, 3)
        self.force = np.array([], dtype=float).reshape(0, 3)

        # Prepare for the next time step
        return body

    