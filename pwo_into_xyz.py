# main
import numpy as np
import matplotlib.pyplot as plt
from class_iteration import iteration

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

    with open('Setup.txt', 'r') as fset:
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

    return filename, Rmax, at1, at2, N

def xyz_gen(fout, fin, RDF, filename):
    '''
    This function checks the fin file (pwo form) to extract all the values and generate the corresponding fout file (xyz form)

    Args:
    - fout: xyz file
    - fin: pwo file
    - RDF: list with several information needed to calculate the RDF and define the size of the x-axis of the histogram
    - filename: the name of the pwo file that we'll convert into a xyz file
    '''
    iteration_obj = iteration()
    iteration_obj.set_RDF(RDF)
    reading_atomic_species = False
    counting_atoms = False
    generation_switch = False
    storage_switch = False
    
    
    for line in fin:

        # define the setup information
        
        if 'number of atoms/cell' in line:
            iteration_obj.n_at = int(line.split()[4])
        elif 'celldm(1)= ' in line:
            iteration_obj.celldim = float(line.split()[1])
            iteration_obj.Alat_to_Angstrom()
        elif 'a(1)' in line:
            x, y, z = map(float, line.split()[3:6])
            a = np.multiply([x, y, z], iteration_obj.alat_to_angstrom)
            iteration_obj.ax = a
            iteration_obj.L[0] = np.linalg.norm(a)
        elif 'a(2)' in line:
            x, y, z = map(float, line.split()[3:6])
            a = np.multiply([x, y, z], iteration_obj.alat_to_angstrom)
            iteration_obj.ay = a
            iteration_obj.L[1] = np.linalg.norm(a)
        elif 'a(3)' in line:
            x, y, z = map(float, line.split()[3:6])
            a = np.multiply([x, y, z], iteration_obj.alat_to_angstrom)
            iteration_obj.az = a
            iteration_obj.L[2] = np.linalg.norm(a)

        # defining the groups

        elif 'atomic species   valence' in line:
            reading_atomic_species = True
        elif reading_atomic_species == True:
            _line = line.split()
            if _line == []:
                reading_atomic_species = False
            else:
                iteration_obj.add_group(_line[0], _line[2], RDF)

        # defining the index of atoms and the initial positions

        elif 'site n.' in line:
            counting_atoms = True
        elif counting_atoms:
            _line = line.split()
            if _line == []:
                counting_atoms = False
            else:
                iteration_obj.count_group(_line)

        # defining the configuration at each time

        elif 'Self-consistent Calculation' in line:
            if generation_switch == True:
                fout.writelines(["%s\n" % i for i in iteration_obj.single_frame(RDF)])
            generation_switch = True
            storage_switch = True
        elif storage_switch == True:
            iteration_obj.store(line)
    
    iteration_obj.normalization()
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(iteration_obj.R, iteration_obj.count, label='Power spectrum')
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