# class_iteration.py

import numpy as np
from class_group import group

class iteration:
    n_at = 0
    ax = [0, 0, 0] #A
    ay = [0, 0, 0] #A
    az = [0, 0, 0] #A
    t_past = 0
    t = 0 #ps
    U_pot = 0
    input = ''
    block = []
    groups = []
    switch = False

    pass

    def first_two_lines(self):
        '''
        This method generate the first part of the two starting line of an interation in the xyz file, accoring to the setted parameters
        Args:
            None
            
        Return:
            None
        '''
        input = f'{self.n_at}\nLattice = \"{self.ax[0]}, {self.ax[1]}, {self.ax[2]}, {self.ay[0]}, {self.ay[1]}, {self.ay[2]}, {self.az[0]}, {self.az[1]}, {self.az[2]}\"'
        return 
    
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
                #print(line[8])
                elm.add_id(line[0])
                elm.add_position_past( line[7], line[8], line[9])
                break
        return
    
    def forces(self, line):
        _line = line.split()
        #print(_line)
        for elm in self.groups:
            if _line[2] in elm.id:
                elm.add_force( _line[7], _line[8], _line[9])
                break
        return
    
    def position(self, line):
        _line = line.split()
        #print(_line)
        for elm in self.groups:
            if _line[0] == elm.type:
                elm.add_position( _line[1], _line[2], _line[3])
                break
        return

    def single_frame(self):
        for id, line in enumerate(self.block):

            if 'time      =' in line:
                self.t_past = self.t
                self.t = float(line.split()[2]) * 1e3 #ps 

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

         
#2000
        
        self.block = []

        return