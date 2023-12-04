# class_group
import numpy as np

class group:
    '''
    This class represents a group of atoms used to divide the whole atoms in a molecular dinamics. In particular the methods impelemented calculate several properties of the group
    
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
     - DOF: degree of freedom of the group
    '''

    def __init__(self, type, mass):
        self.type = type
        self.mass = mass
        self.id_group = 0
        for i in self.type:
            if i.isdigit():
                self.id_group = int(i)
        self.id = np.array([], dtype=int)
        self.N = 0
        self.force = np.array([], dtype=float).reshape(0, 3)
        self.position_past = np.array([], dtype=float).reshape(0, 3)
        self.position = np.array([], dtype=float).reshape(0, 3)
        self.velocity = np.array([], dtype=float).reshape(0, 3)
        self.Ek = 0
        self.Ftot = np.array([0, 0, 0], dtype=float)
        self.T = 0
        self.DOF = 0
        
    def add_id(self, i):
        '''
        This method add to the list id the number i, that represents the identification number of the atom from the pwo file

        Parameters:
         - i: input integer number
        
        '''
        self.id = np.append(self.id, i)
    
    def add_position_past(self, x, y, z):
        '''
        This method add to the list position_past the coordinates of the considered atom from the starting configuration

        Parameters:
         - x, y, z: coordinates of the considered atom
        
        '''
        self.position_past = np.vstack([self.position_past, np.array([x, y, z], dtype=float)])
    
    def add_position(self, x, y, z):
        '''
        This method add to the list position_past the coordinates of the considered atom from the current time step

        Parameters:
         - x, y, z: coordinates of the considered atom
        
        '''
        self.position = np.vstack([self.position, np.array([x, y, z], dtype=float)])
    
    def add_force(self, x, y, z):
        '''
        This method add to the list force the coordinates of force acting on the considered atom from the current time step

        Parameters:
         - x, y, z: force coordinates of the considered atom
        
        '''
        self.force = np.vstack([self.force, np.array([x, y, z], dtype=float)])

    def Velocity(self, Dt):
        '''
        This method generates the bidimensional array representing the velocity of each atoms. The velocities are calculated as the difference between the current position and the one of the time step before over the time intervall

        Parameters:
         - Dt: time intervall
        
        '''
        self.velocity = (np.array(self.position, dtype=float) - np.array(self.position_past, dtype=float)) / Dt

    def generate(self, Dt):
        '''
        In this method several things are performed:
         - the array of the velocity is generated
         - the kinetic energy is caluìculated
         - the temperature is calculated
         - the total force is calculated
         - the line printed in the ouput file is generated
        Then the group is prepared for the next time step:
         - position_past becams equal to position
         - force and positioon are recasted as empty array
          
        Parameters:
         - Dt: time intervall
        '''
        self.Velocity(Dt)
        self.Ek = 0.5 * float(self.mass) * np.sum(np.linalg.norm(self.velocity, axis=1)**2) * 0.0001036426948415943
        self.DOF = 3* self.N
        self.T = (2 * self.Ek) / (self.DOF * 8.617333262145e-5)
        self.Ftot = np.sum(self.force, axis=0)
        body = np.array([], dtype=str)
        for i in range(self.N):
            body = np.append(body, f'{self.type}\t  {self.position[i][0]}\t  {self.position[i][1]}\t  {self.position[i][2]}\t  {self.velocity[i][0]}\t  {self.velocity[i][1]}\t  {self.velocity[i][2]}\t  {self.force[i][0]}\t  {self.force[i][1]}\t  {self.force[i][2]}\t  {self.id_group}')
        self.position_past = self.position
        self.force = np.array([], dtype=float).reshape(0, 3)
        self.position = np.array([], dtype=float).reshape(0, 3)
        return body