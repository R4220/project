# class RDF
import numpy as np
import matplotlib.pyplot as plt

class RDF:
    def __init__(self, filename, Rmax, atoms, N):
        self.filename = filename
        self.Rmax = Rmax
        self.type = atoms
        self.N = N

        self.count = np.zeros(N)
        self.R = np.linspace(0, Rmax, N)
        self.dR = self.R[1]
        self.norm = np.multiply([((i + self.dR)**3 - i**3) for i in self.R[:N]], np.pi * 4 /3)
        self.condition = (atoms[0] == atoms[1])
        self.at1 = np.array([], dtype=float).reshape(0, 3)
        self.N1 = 0
        self.at2 = np.array([], dtype=float).reshape(0, 3)
        self.N2 = 0

        self.matrix = None      

    def add_position1(self, x, y, z):
        self.at1 = np.vstack([self.at1, np.array([x, y, z], dtype=float)])
    
    def add_position2(self, x, y, z):
        self.at2 = np.vstack([self.at2, np.array([x, y, z], dtype=float)])

    def RDF(self):
        '''
        Calculate distances and store them in the list dist for RDF calculation. 
        Two cases are present, defined by the boolean variable e:
        - If e is True, RDF is calculated within the same group.
        - If e is False, RDF is calculated between two different groups.

        Parameters:
        - RDF: List containing parameters for RDF calculation.
     
        Return:
        - dist: Array of distances.
        '''
        dist = []

        if self.condition:
            n = len(self.at1)
            #print(self.at1)
            #print(self.at2)
            rpos = self.at1
            #pos = rpos

            for i, _pos in enumerate(self.at1):
                rpos[i] = np.dot(np.linalg.inv(self.matrix), _pos)
                #print(i,'\npos\n', _pos, '\nrpos\n', rpos[i])
            
            for k in range(n - 1):
                #print(rpos)
                rdiff = rpos[k+1:] - rpos[k]
                #print(rdiff)
                int_pos = np.round(rdiff)
                #print(int_pos)
                rdiff = rdiff - int_pos
                #print(rdiff)
                diff = rdiff
                for i, _rpos in enumerate(rdiff):
                    diff[i] = np.dot(self.matrix, _rpos)
                    #print(_rpos)
                #print(diff)
                r = np.linalg.norm(diff)
                #print(r)
                  #print('iter')
                dist.extend(r[r < self.Rmax])
            #print(dist)
        else:
            n1 = len(self.at1)
            n2 = len(self.at2)
            #print(self.at1)
            #print(self.at2)

            rpos1 = self.at1
            rpos2 = self.at2

            for i, _pos in enumerate(self.at1):
                rpos1[i] = np.dot(np.linalg.inv(self.matrix), _pos)
            for i, _pos in enumerate(self.at2):
                rpos2[i] = np.dot(np.linalg.inv(self.matrix), _pos)

            for i in rpos1:
                rdiff = i - rpos2
                rdiff = rdiff - np.round(rdiff)
                diff = rdiff
                for i, _rpos in enumerate(rdiff):
                    diff[i] = np.dot(self.matrix, _rpos)
                    #print(_rpos)
                #print(diff)
                r = np.linalg.norm(diff)
                #print(r)
                #print('iter')
                dist.extend(r[r < self.Rmax])
            #print(dist)
        
        self.N1 = len(self.at1)
        self.N2 = len(self.at2)
        
        self.at1 = np.array([], dtype=float).reshape(0, 3)
        self.at2 = np.array([], dtype=float).reshape(0, 3)

        #dist = np.sort(dist)
        self.count += np.histogram(dist, bins=self.N, range=(0, self.Rmax))[0]
        #return dist
    

    def normalization(self, iteration_obj):
        '''
        Normalize the count variable.
        The normalization process involves dividing the count variable by a normalization factor.
        The normalization factor (rho) is calculated based on the volume of the simulation cell and the number of atoms involved in the RDF calculation. The specific calculation depends on whether the RDF is calculated between atoms of the same group (e=True) or different groups (e=False).
        - If e=True, rho = (Number of atoms in group 1)^2 / Volume
        - If e=False, rho = (Number of atoms in group 1 * Number of atoms in group 2) / Volume

        The normalized count variable is then obtained by dividing the original count by the normalization factor.
        '''
        V = np.dot(iteration_obj.az, np.cross(iteration_obj.ax, iteration_obj.ay)) ** 2

        if self.condition == True:
            print(V, '\n' , self.N1,'\n', iteration_obj.N_iteration)
            rho_1 = V / (self.N1 * (self.N1 -1 ) * iteration_obj.N_iteration)
        else:
            rho_1 = V / (self.N1 * self.N2 * iteration_obj.N_iteration)

        self.norm = rho_1 * self.norm
        self.count = np.divide(self.count, self.norm)
        ''' Check if the normalization process is correctly implemented '''

    def plot(self):
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.R, self.count, label='Power spectrum')
        ax.set_xlabel('r ($A$)')
        ax.set_ylabel('g(r)')
        ax.grid()
        plt.savefig(f'{self.filename}.png')
        plt.show()