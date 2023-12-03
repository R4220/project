# class_group

import numpy as np

class group:

    id_group = 0
    id = []
    N = 0
    force =  np.array([]).reshape(0, 3)
    position_past = np.array([]).reshape(0, 3)
    position = np.array([]).reshape(0, 3)
    velocity = np.array([]).reshape(0, 3)

    Ek = 0
    Ftot = [0, 0, 0]
    T = 0
    DOF = 0

    Ry_to_eV = 13.60570398  # 1 Ry = 13.60570398 eV
    eV_j = 1.60e-19         # 1 eV = 1.60e-19 j

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
    
    def check_id(self, id):
        return self.id == id
    
    def add_position_past(self, x, y, z):
        self.position_past = np.vstack([self.position_past, np.array([x, y, z], dtype=float)])
        return
    
    def add_position(self, x, y, z):
        self.position = np.vstack([self.position, np.array([x, y, z], dtype=float)])
        return
    
    def add_force(self, x, y, z):
        #print(self.force)
        self.force = np.vstack([self.force, np.array([x, y, z], dtype=float)])
        return
    
    def Velocity_(self, Dt):
        #print(np.array(self.position, dtype=float) - np.array(self.position_past, dtype=float))
        #print(self.position_past, '\n', Dt)
        self.velocity = (np.array(self.position, dtype=float) - np.array(self.position_past, dtype=float)) / Dt

        return
    
    def generate(self, Dt):
        self.Velocity_(Dt)
        #print(self.velocity)
        #print(np.linalg.norm(self.velocity, axis=1), '\n')
        #print(np.linalg.norm(self.velocity, axis=1)**2, '\n')
        self.Ek = 0.5 * float(self.mass) * np.sum(np.linalg.norm(self.velocity, axis=1)**2) * 0.0001036426948415943
        # (A/ps -> m/s)^2 = 1e4
        # g/mol -> Kg/mol = 1e-3
        # eV -> j = 1.602176634e-19
        # Na = 6.02214086e23
        # Kb = 1.3806503e-23 jK-1 = 8.617333262145e-5 eVK-1
        #print(10 /(1.602176634e-19 * 6.02214086e23))
        self.DOF = self.N
        #print(self.DOF)
        self.T = (2 * self.Ek) / (self.DOF * 8.617333262145e-5)
        #print('bef')
        #print(self.velocity)
        self.Force_tot = np.sum(self.force, axis=0)
        #print('aft')
        body = []
        #print(self.force)
        for i in range(self.N):
            #print(f'{self.type}    {self.position[i][0]}   {self.position[i][1]}   {self.position[i][2]}   {self.velocity[i][0]}   {self.velocity[i][1]}   {self.velocity[i][2]}  {self.force[i][0]}   {self.force[i][1]}   {self.force[i][2]}')
            body = np.append(body, f'{self.type}\t  {self.position[i][0]}\t  {self.position[i][1]}\t  {self.position[i][2]}\t  {self.velocity[i][0]}\t  {self.velocity[i][1]}\t  {self.velocity[i][2]}\t  {self.force[i][0]}\t  {self.force[i][1]}\t  {self.force[i][2]}\t  {self.id_group}')

        self.position_past = self.position
        self.force =  np.array([]).reshape(0, 3)
        self.position = np.array([]).reshape(0, 3)

        return body