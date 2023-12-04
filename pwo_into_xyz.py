# main
import numpy as np
import matplotlib.pyplot as plt
from class_iteration import iteration

def setup():
    '''
    This function extracts the setup information from Setup.txt

    Args:
        None

    Returns:
        filename: the name of the pwo file that we'll convert into a xyz file
        Rmax: maximum distance to consider while RDF is calculated
        at1, at2: atoms between which the RDF is calculated

    '''

    with open('Setup.txt', 'r') as fset:
        for line in fset:
            _line = line.split()
            if 'Filename:' in line:
                filename = _line[1]
            if 'Rmax:' in line:
                Rmax = float(_line[1])
            if 'Particles:' in line:
                at1 = _line[1]
                at2 = _line[2]
            if 'N:' in line:
                N = int(_line[1])

    return filename, Rmax, at1, at2, N

def xyz_gen(fout, fin, RDF, filename):
    '''
    Thin function check the fin file (pwo form) in order to estract all the values and generate the correspective fout file (xyz form)

    Args:
        fout:   xyz file
        fin:    pwo file

    '''
    iter = iteration()
    iter.set_RDF(RDF)
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
            a = np.multiply([x, y, z], iter.alat_to_angstrom)
            iter.ax = a
            iter.L[0] = np.linalg.norm(a)
            #print(iter.ax)
        elif 'a(2)' in line:
            x, y, z = map(float, line.split()[3:6])
            a = np.multiply([x, y, z], iter.alat_to_angstrom)
            iter.ay = a
            iter.L[1] = np.linalg.norm(a)
            #print(iter.ay)
        elif 'a(3)' in line:
            x, y, z = map(float, line.split()[3:6])
            a = np.multiply([x, y, z], iter.alat_to_angstrom)
            iter.az = a
            iter.L[2] = np.linalg.norm(a)
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
                iter.add_group(_line[0], _line[2], RDF)#, idx)
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
                fout.writelines(["%s\n" % i for i in iter.single_frame(RDF)])
            e_gen = True
            e_storage = True
        elif e_storage == True:
            iter.store(line)

        #Self constant calculation

        #print(n_at, bravais, lattice_param, V)
        #fout.writelines()
    
    iter.normalization()
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(iter.R, iter.count, label='Power spectrum')
    ax.set_xlabel('r ($A$)')
    ax.set_ylabel('g(r)')
    ax.grid()
    plt.savefig(f'{filename}.png')
    plt.show()
    
    return

if __name__ == "__main__":
    filename, Rmax, at1, at2, N = setup()
    RDF = [Rmax, at1, at2, N]
    with open(filename + '.xyz', "w+") as fout:
        with open(filename + '.pwo', 'r') as fin:
            xyz_gen(fout, fin, RDF, filename)