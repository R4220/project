# main
import numpy as np

from class_iteration import iteration

def setup():
    '''
    This function extracts the setup information from Setup.txt

    Args:
        None

    Returns:
        filename: the name of the pwo file that we'll convert into a xyz file

    '''

    with open('Setup.txt', 'r') as fset:
        for line in fset:
            if 'Filename:' in line:
                filename = line.split()[1]

    return filename

def xyz_gen(fout, fin):
    '''
    Thin function check the fin file (pwo form) in order to estract all the values and generate the correspective fout file (xyz form)

    Args:
        fout:   xyz file
        fin:    pwo file

    '''
    iter = iteration()
    e_storage = False
    e_gen = False
    e_group = False
    e_count = False
    #idx = 1
    
    for line in fin:

        '''define the setup information'''

        if 'number of atoms/cell' in line:
            iter.n_at = int(line.split()[4])
            #print(iter.n_at)
        elif 'celldm(1)= ' in line:
            iter.celldim = float(line.split()[1])
            iter.Alat_to_Angstrom()
            #print(iter.celldim)
        elif 'a(1)' in line:
            x, y, z = map(float, line.split()[3:6])
            iter.ax = np.multiply([x, y, z], iter.alat_to_angstrom)
            #print(iter.ax)
        elif 'a(2)' in line:
            x, y, z = map(float, line.split()[3:6])
            iter.ay = np.multiply([x, y, z], iter.alat_to_angstrom)
            #print(iter.ay)
        elif 'a(3)' in line:
            x, y, z = map(float, line.split()[3:6])
            iter.az = np.multiply([x, y, z], iter.alat_to_angstrom)
            #print(iter.az)

        '''defining the groups'''
        if 'atomic species   valence' in line:
            e_group = True
        elif e_group == True:
            #print(line)
            _line = line.split()
            if _line == []:
                e_group = False
            else:
                iter.add_group(_line[0], _line[2])#, idx)
                #idx += 1
        elif 'site n.' in line:
            e_count = True
        elif e_count == True:
            _line = line.split()
            if _line == []:
                e_count = False
            else:
                #print(_line)
                iter.count_group(_line)#, idx)
            

        '''defining the configuration at each time'''

        if 'Self-consistent Calculation' in line:
            if e_gen == True:
                fout.writelines(["%s\n" % i for i in iter.single_frame()])
            e_gen = True
            e_storage = True
        elif e_storage == True:
            iter.store(line)

        #Self constant calculation

        #print(n_at, bravais, lattice_param, V)
        #fout.writelines()
    return

if __name__ == "__main__":
    filename = setup()
    with open(filename + '.xyz', "w+") as fout:
        with open(filename + '.pwo', 'r') as fin:
            xyz_gen(fout, fin)