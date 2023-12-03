# class_iteration.py

import numpy as np
from class_group import group

class iteration:
    n_at = 0
    celldim = 0
    ax = [0, 0, 0] #A
    ay = [0, 0, 0] #A
    az = [0, 0, 0] #A
    t_past = 0 #ps
    t = 0 #ps
    U_pot = 0
    input = ''
    block = []
    groups = []
    switch = False

    # Conversion
    Ry_to_eV = 13.60570398  # 1 Ry = 13.60570398 eV
    #eV_j = 1.60e-19         # 1 eV = 1.60e-19 j
    bohr_angstrom =	0.529177249 # 1 A = 0.529177249 Bohr
    alat_to_angstrom = 0


    pass

    def Alat_to_Angstrom(self):
        self.alat_to_angstrom = self.celldim * self.bohr_angstrom
        return 

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
        return
    
    def add_group(self, _type, _mass):#, idx):
        '''
        This method add to groups the identification of the group
        
        Arg:
            _type: name of the atom
            _mass: mass of the atom

        Return:
            None
        '''
        
        self.groups = np.append(self.groups, group(_type, _mass))#, idx))
        #print(_mass, _type)#, idx)
        return
    
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
                #print(line[0])
                elm.add_id(int(line[0]))
                #elm.add_position_past(line[6] * self.alat_to_angstrom, line[7] * self.alat_to_angstrom, line[8] * self.alat_to_angstrom)
                elm.add_position_past(float(line[6]) * self.alat_to_angstrom, float(line[7]) * self.alat_to_angstrom, float(line[8]) * self.alat_to_angstrom)
                break
        return
    
    def forces(self, line):
        '''
        This method extract the force at a determined t
        
        Args:
            line: the line in which there are the coordinates of the force
        
        Return:
            None
        
        '''
        _line = line.split()
        #print(_line)
        for elm in self.groups:
            #print(_line[2])
            if int(_line[1]) in elm.id:
                #print(_line[6], _line[7], _line[8])
                elm.add_force(_line[6], _line[7], _line[8])
                break
        return
    
    def position(self, line):
        '''
        This method extract the position at a determined t
        
        Args:
            line: the line in which there are the coordinates of the position
        
        Return:
            None
        
        '''
        _line = line.split()
        #print(_line)
        for elm in self.groups:
            if _line[0] == elm.type:
                elm.add_position( _line[1], _line[2], _line[3])
                break
        return

    def text(self):
        text = [f'{self.n_at}', 
                f"Lattice(Ang)=\"{self.ax[0]}, {self.ax[1]}, {self.ax[2]}, {self.ay[0]}, {self.ay[1]}, {self.ay[2]}, {self.az[0]}, {self.az[1]}, {self.az[2]}\" t(ps)={self.t} Epot(eV)={self.U_pot / self.Ry_to_eV }"]
        body = []
        for elm in self.groups:
            _body = elm.generate(self.t - self.t_past)
            elm.DOF = elm.N
            text[1] += f" Ek{elm.id_group}(ev)={elm.Ek} DOF{elm.id_group}={elm.DOF} T{elm.id_group}(K)={elm.T} Ftot{elm.id_group}(pN)=\"{elm.Ftot[0]}, {elm.Ftot[1]}, {elm.Ftot[2]}\""

            body = np.append(body, _body)
        
        text = np.append(text, body)

        return text

    def single_frame(self):
        for line in self.block:

            if 'time      =' in line:
                self.t_past = self.t
                self.t = float(line.split()[2]) #ps 

            if '!    total energy' in line:
                self.U_pot = float(line.split()[4]) #Ry
            
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
        
        text =self.text() 

        self.block = []

        return text