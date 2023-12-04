# class_iteration.py

import numpy as np

from class_group import group

class iteration:

    def __init__(self):
        self.n_at = 0
        self.celldim = 0
        self.ax = [0, 0, 0]
        self.ay = [0, 0, 0]
        self.az = [0, 0, 0]
        self.L = [0, 0, 0]
        self.t_past = 0
        self.t = 0
        self.U_pot = 0
        self.input = ''
        self.block = np.array([], dtype=str)
        self.groups = []#np.array([], dtype=str)
        self.RDF_atoms = []#np.array([], dtype=str)
        self.switch = False
        self.Ry_to_eV = 13.60570398
        self.bohr_angstrom = 0.529177249
        self.alat_to_angstrom = 0
        self.count = np.array([], dtype = int)
        self.R = np.array([], dtype=float)
        self.dR = 0.
        self.norm = np.array([], dtype=float)
        self.e = False

    def Alat_to_Angstrom(self):
        self.alat_to_angstrom = self.celldim * self.bohr_angstrom

    '''def first_two_lines(self):
        
        This method generate the first part of the two starting line of an interation in the xyz file, accoring to the setted parameters
        
        Args:
            None
            
        Return:
            None
        

        input = f'{self.n_at}\nLattice(Ang)=\"{self.ax[0]}, {self.ax[1]}, {self.ax[2]}, {self.ay[0]}, {self.ay[1]}, {self.ay[2]}, {self.az[0]}, {self.az[1]}, {self.az[2]}\"'
        return '''
    
    def store(self, line):
        '''
        This method add to block the line in the current iteration, if it is not empty
        
        Arg:
            line: the line readed
        
        Return:
            None
        '''
        if line != '':
            self.block = np.append(self.block, line)
    
    def add_group(self, _type, _mass, RDF):#, idx):
        '''
        This method add to groups the identification of the group
        
        Arg:
            _type: name of the atom
            _mass: mass of the atom

        Return:
            None
        '''
        
        at = group(_type, _mass)
        self.groups = np.append(self.groups, at )
        if _type in RDF:
            self.RDF_atoms = np.append(self.RDF_atoms, at)

    
    def count_group(self, line):
        '''
        This method counts the number of atoms for each group
        
        Args:
            line: the readed line of the file
            
        Return:
            None
            '''
        for elm in self.groups:
            if elm.type == line[1]:
                elm.N += 1
                elm.add_id(int(line[0]))
                elm.add_position_past(float(line[6]) * self.alat_to_angstrom, float(line[7]) * self.alat_to_angstrom, float(line[8]) * self.alat_to_angstrom)
                break
    
    def forces(self, line):
        '''
        This method extract the force at a determined t
        
        Args:
            line: the line in which there are the coordinates of the force
        
        Return:
            None
        
        '''
        _line = line.split()
        for elm in self.groups:
            if int(_line[1]) in elm.id:
                elm.add_force(_line[6], _line[7], _line[8])
                break
    
    def position(self, line):
        '''
        This method extract the position at a determined t
        
        Args:
            line: the line in which there are the coordinates of the position
        
        Return:
            None
        
        '''
        _line = line.split()
        for elm in self.groups:
            if _line[0] == elm.type:
                elm.add_position(_line[1], _line[2], _line[3])
                break

    def text(self):
        text = np.array([f'{self.n_at}', f"Lattice(Ang)=\"{self.ax[0]}, {self.ax[1]}, {self.ax[2]}, {self.ay[0]}, {self.ay[1]}, {self.ay[2]}, {self.az[0]}, {self.az[1]}, {self.az[2]}\" t(ps)={self.t} Epot(eV)={self.U_pot / self.Ry_to_eV }"])
        body = np.array([], dtype=str)
        for elm in self.groups:
            _body = elm.generate(self.t - self.t_past)
            elm.DOF = elm.N
            text[1] += f" Ek{elm.id_group}(ev)={elm.Ek} DOF{elm.id_group}={elm.DOF} T{elm.id_group}(K)={elm.T} Ftot{elm.id_group}(pN)=\"{elm.Ftot[0]}, {elm.Ftot[1]}, {elm.Ftot[2]}\""
            body = np.append(body, _body)

        text = np.append(text, body)
        return text
    
        # estrazione delle coordinate degli atomi
    '''def coordinates(Block, ind):
        ''fa cose''
        positions = [Block[i] for i in ind]
        for idx, i in enumerate(positions):
            matches = re.findall(r'(-?\d+(\.\d+)?)', i)
            positions[idx] = [float(match[0]) for match in matches[:3]]
        return positions'''
    
    
    # calcolo le varie distanze di un blocco
    def RDF(self, RDF):
        dist = []
        

        if self.e == True:  # same species
            pos = self.RDF_atoms[0].position
            N = len(pos)
            for k in range(N - 1):
                _pos = pos[k+1:] - pos[k]
                n = np.divide(_pos, self.L)
                print(n)#####quiii
                for i in n:
                    i = round(i)
                r = np.linalg.norm(_pos - np.multiply(self.L, n))
                dist.extend(r[r < RDF[0]])

        else:
            pos_1 = self.RDF_atoms[0].position
            pos_2 = self.RDF_atoms[1].position
            for i in pos_1:
                r = np.linalg.norm(i - pos_2, axis=1)
                dist.extend(r[r < RDF[0]])

        dist = np.sort(dist)
        return dist

    def set_RDF(self, RDF):
        self.count = np.zeros(RDF[3])
        self.R = np.linspace(0, RDF[0], RDF[3])
        self.dR = self.R[1]
        self.norm = np.multiply([((i + self.dR)**3 - i**3) for i in self.R[:RDF[3]]], np.pi * 4 /3)
        self.e = (RDF[1] == RDF[2])

    def istogram(self, RDF):
        d = self.RDF(RDF)
        self.count += np.histogram(d, bins=RDF[3], range=(0, RDF[0]))[0]
        
        
    
        '''for i in blocks:
            d = RDF(i, ind, e, r_max)
            count += np.histogram(d, bins=N, range=(0, r_max))[0]

        print('distances generated')
    
        V = 4 * np.pi * R[-1] #Lx, Ly, Lz
        rho_1 = V/n_at
        count = np.divide(count, norm) * rho_1
'''


    def single_frame(self, RDF):
        for line in self.block:
            if 'time      =' in line:
                self.t_past = self.t
                self.t = float(line.split()[2])
            if '!    total energy' in line:
                self.U_pot = float(line.split()[4])
            if '   force =    ' in line:
                self.forces(line)
            if 'ATOMIC_POSITIONS' in line:
                self.switch = True
            elif self.switch == True:
                _line = line.split()
                if _line == []:
                    self.switch = False
                else:
                    self.position(line)

        self.istogram(RDF)

        text = self.text()
        self.block = []
        return text
    
    def normalization(self):
        V = (self.L[0] * self.L[1] * self.L[2])**2
        print(self.L)
        print(self.RDF_atoms[0].N, self.RDF_atoms[1].N)
        if self.e == True:
            rho = self.RDF_atoms[0].N **2/(V)
        else:
            rho = self.RDF_atoms[0].N * self.RDF_atoms[1].N /(V)
        
        self.norm = rho * self.norm
        self.count = np.divide(self.count, self.norm)
    