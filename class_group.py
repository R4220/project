# class_group

import numpy as np

class group:

    id_group = 0
    id = []
    N = 0
    force = []
    position_past = []
    position = []
    velocity = []

    Ek = 0
    Ftot = [0, 0, 0]
    T = 0
    DOF = 0

    def __init__(self, type, mass):
        #estrai il numero dalla stringa 
        self.type = type
        self.mass = mass
        for i in self.type:
            if i.isdigit():
                self.id_group = int(i)

    def add_id(self, i):
        self.id = np.append(self.id, i)
        return
    
    def add_position_past(self, x, y, z):
        self.position_past = np.append(self.position_past, [x, y, z])
        return
    
    def add_position(self, x, y, z):
        self.position = np.append(self.position, [x, y, z])
        return
    
    def add_force(self, x, y, z):
        self.force = np.append(self.force, [x, y, z])
        return
    
    def Velocity_(self, Dt):
        self.velocity = (self.position - self.position_past)/Dt
        return
    
    def generate(self, Dt):
        self.Velocity_(Dt)

        self.Ek = 0.5 * self.mass * np.sum(np.linalg.norm(self.velocity, axis=1)**2)

        self.T = (2 * self.Ek) / (self.DOF * 1.380649e-23) # ocio alle unità di misura, sono sbagliate ora

        self.Force_tot = np.sum(self.force, axis=0)

        body = []

        for i in self.N:
            body = np.append(body, f'{self.type}    {self.position[i][0]}   {self.position[i][1]}   {self.position[i][2]}   {self.velocity[i][0]}   {self.velocity[i][1]}   {self.velocity[i][2]}   {self.force[i][0]}   {self.force[i][1]}   {self.force[i][2]}')

        self.position_past = self.position
        self.position = []
        self.force = []

        return body