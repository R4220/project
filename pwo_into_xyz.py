# main
import numpy as np
import matplotlib.pyplot as plt
from class_iteration import iteration
from class_group import group
from class_RDF import RDF

def setup():
    '''
    Extract setup information from Setup.txt.

    Returns:
    - filename (str): Name of the pwo file to convert into an xyz file.
    - Rmax (float): Maximum distance to consider while calculating RDF.
    - at1, at2 (str): Atoms between which the RDF is calculated.
    - N (int): Number of atoms in the simulation.

    '''
    filename = ''

    Rmax = 0.0
    at1, at2 = '', ''
    N = 0
    groups = []

    with open('Setup.txt', 'r') as fset:
        id = 0
        move = True
        for line in fset:
            _line = line.split()

            if 'Filename:' in line:
                filename = _line[1]
            if 'Rmax:' in line:
                Rmax = float(_line[1])
            if 'Particles:' in line:
                at1, at2 = _line[1], _line[2]
            if 'N:' in line:
                N = int(_line[1])
            if 'Group' in line:
                id = _line[1]
                if 'ext' in line:
                    move = False
                else:
                    move = True
                _line = fset.readline().split()
                groups = np.append(groups, group(_line, id, move))

    return filename, Rmax, at1, at2, N, groups

def preamble(fin, iteration_obj, preamble_switch):
    '''
    This method extract the information form the preamble of the input file
    '''
    line = fin.readline()

    if 'number of atoms/cell' in line:
        iteration_obj.n_atoms = int(line.split()[4])
        preamble_switch[0] = False
    elif 'number of atomic types' in line:
        iteration_obj.n_type = int(line.split()[5])
        preamble_switch[1] = False
    elif 'celldm(1)= ' in line:
        iteration_obj.Alat_to_Angstrom(float(line.split()[1]))
        preamble_switch[2] = False
    elif 'a(1)' in line:
        x, y, z = map(float, line.split()[3:6])
        a = np.multiply([x, y, z], iteration_obj.alat_to_angstrom)
        iteration_obj.ax = a
        preamble_switch[3] = False
    elif 'a(2)' in line:
        x, y, z = map(float, line.split()[3:6])
        a = np.multiply([x, y, z], iteration_obj.alat_to_angstrom)
        iteration_obj.ay = a
        preamble_switch[4] = False
    elif 'a(3)' in line:
        x, y, z = map(float, line.split()[3:6])
        a = np.multiply([x, y, z], iteration_obj.alat_to_angstrom)
        iteration_obj.az = a
        preamble_switch[5] = False
    
    # defining the groups
    elif 'atomic species   valence' in line:
        for _ in range(iteration_obj.n_type):
            line = fin.readline().split()
            iteration_obj.set_mass(line[0], line[2])
        preamble_switch[6] = False

    # defining the index of atoms and the initial positions
    elif 'site n.' in line:
        for _ in range(iteration_obj.n_atoms):
            line = fin.readline()
            iteration_obj.count_group(line.split())
        preamble_switch[7] = False
        iteration_obj.set_DOF()

    return preamble_switch

def extract_forces(line, iteration_obj):
    while not 'negative rho' in line:
        line = fin.readline()
    
    for _ in range(iteration_obj.n_atoms):
        line = fin.readline().split()
        iteration_obj.forces(line)

def extract_positions(line, iteration_obj, istogram):
    for _ in range(iteration_obj.n_atoms):
        line = fin.readline().split()
        iteration_obj.positions(line, istogram)

def body(iteration_obj, istogram):
    generation_switch = [True] * 4
    for line in fin:
        
        # extracting the potential energy
        if '!    total energy' in line:
            iteration_obj.U_pot = float(line.split()[4])
            generation_switch[0] = False

        # extracting the forces on atoms
        if 'Forces a' in line:
            extract_forces(line, iteration_obj)
            generation_switch[1] = False

        # extracting the time step
        elif 'Time step' in line:
            iteration_obj.dt = float(line.split()[3])  * 4.8378e-5 #s

        # extracting the iteration number
        elif 'Entering' in line:
             iteration_obj.N_iteration = int(line.split()[4])
             generation_switch[2] = False
             
        # extracting the atomic position
        elif 'ATOMIC_POSITIONS' in line:
            extract_positions(line, iteration_obj, istogram)
            generation_switch[3] = False
            
        if not any(generation_switch):
            fout.writelines(["%s\n" % i for i in iteration_obj.single_frame(istogram)])
            #print('ok')
            istogram.RDF()
            
            generation_switch = [True] * 4

def xyz_gen(fout, fin, RDF_, groups ):
    '''
    This function checks the fin file (pwo form) to extract all the values and generate the corresponding fout file (xyz form)

    Args:
    - fout: xyz file
    - fin: pwo file
    - RDF: list with several information needed to calculate the RDF and define the size of the x-axis of the histogram
    - filename: the name of the pwo file that we'll convert into a xyz file
    '''
    iteration_obj = iteration(groups)
    preamble_switch = [True] * 8

    # define the setup information
    while any(preamble_switch):
        preamble_switch = preamble(fin, iteration_obj, preamble_switch)
    
    istogram = RDF(RDF_[0], RDF_[1], [RDF_[2], RDF_[3]], RDF_[4])
    #print(iteration_obj.ax, iteration_obj.ay, iteration_obj.az)
    istogram.matrix = np.column_stack((iteration_obj.ax, iteration_obj.ay, iteration_obj.az))
    #print(istogram.matrix)
    # defining the configuration at each time
    body(iteration_obj, istogram)

            
        
    istogram.normalization(iteration_obj)
    istogram.plot()
    
    return

if __name__ == "__main__":
    filename, Rmax, at1, at2, N, groups = setup()
    RDF_ = [filename, Rmax, at1, at2, N]
    with open(filename + '.xyz', "w+") as fout:
        with open(filename + '.pwo', 'r') as fin:
            xyz_gen(fout, fin, RDF_, groups)