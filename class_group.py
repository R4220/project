# class_group

import numpy as np

class group:

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
        self.id = np.append(self.id, i)
    
    '''def check_id(self, id):
        return self.id == id'''
    
    def add_position_past(self, x, y, z):
        self.position_past = np.vstack([self.position_past, np.array([x, y, z], dtype=float)])
    
    def add_position(self, x, y, z):
        self.position = np.vstack([self.position, np.array([x, y, z], dtype=float)])
    
    def add_force(self, x, y, z):
        self.force = np.vstack([self.force, np.array([x, y, z], dtype=float)])

    def Velocity(self, Dt):
        self.velocity = (np.array(self.position, dtype=float) - np.array(self.position_past, dtype=float)) / Dt

    def generate(self, Dt):
        self.Velocity(Dt)
        self.Ek = 0.5 * float(self.mass) * np.sum(np.linalg.norm(self.velocity, axis=1)**2) * 0.0001036426948415943
        self.DOF = self.N
        self.T = (2 * self.Ek) / (self.DOF * 8.617333262145e-5)
        self.Ftot = np.sum(self.force, axis=0)
        body = np.array([], dtype=str)
        for i in range(self.N):
            body = np.append(body, f'{self.type}\t  {self.position[i][0]}\t  {self.position[i][1]}\t  {self.position[i][2]}\t  {self.velocity[i][0]}\t  {self.velocity[i][1]}\t  {self.velocity[i][2]}\t  {self.force[i][0]}\t  {self.force[i][1]}\t  {self.force[i][2]}\t  {self.id_group}')
        self.position_past = self.position
        self.force = np.array([], dtype=float).reshape(0, 3)
        self.position = np.array([], dtype=float).reshape(0, 3)
        return body