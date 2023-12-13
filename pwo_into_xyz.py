# main
import numpy as np

from class_iteration import iteration
from class_group import group
from class_graph import graph


def setup():
    """
    Read configuration parameters from 'Setup.txt' and initialize simulation setup.

    Returns
    -------
    tuple
        A tuple containing the following information:
        - filename (str): The name of the file.
        - Rmax (float): The maximum distance for radial distribution function calculation.
        - at1 (str): The type of the first particle for RDF.
        - at2 (str): The type of the second particle for RDF.
        - N (int): The total number of iterations.
        - groups (list): A list of group instances initialized based on configuration.

    Notes
    -----
    This function reads configuration parameters from 'Setup.txt' and initializes the simulation setup.
    It extracts information such as the filename, maximum distance for RDF, particle types for RDF calculation, total iterations, and group instances.

    Examples
    --------
    >>> filename, Rmax, at1, at2, N, groups = setup()
    >>> print(filename)
    'example_file.txt'
    >>> print(Rmax)
    10.0
    >>> print(at1)
    'C'
    >>> print(at2)
    'O'
    >>> print(N)
    1000
    """
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
                _line = fset.readline().split()
                groups = np.append(groups, group(_line, id))

    return filename, Rmax, at1, at2, N, groups


def preamble(fin, iteration_obj, preamble_switch):
    """
    Process the preamble of the input file and set up the iteration object.

    Parameters
    ----------
    fin : file
        The input file.
    iteration_obj : iteration
        The iteration object to be initialized.
    preamble_switch : list
        A list of boolean switches indicating which part of the preamble has been processed.

    Returns
    -------
    list
        Updated preamble switch indicating the processed parts.

    Notes
    -----
    This function processes the preamble of the input file, updating the iteration object attributes based on the information extracted.
    The preamble switch list is used to keep track of which parts have been processed.

    Examples
    --------
    >>> fin = open('input_file.txt', 'r')
    >>> iteration_obj = iteration(groups)
    >>> preamble_switch = [True, True, True, True, True, True, True, True]
    >>> updated_switch = preamble(fin, iteration_obj, preamble_switch)
    >>> print(updated_switch)
    [False, False, False, False, False, False, False, False]
    >>> fin.close()
    """
    line = fin.readline()

    # Extraction of the number of atom in the cell
    if 'number of atoms/cell' in line:
        iteration_obj.n_atoms = int(line.split()[4])
        preamble_switch[0] = False

    # Extraction of the number of atomic types 
    elif 'number of atomic types' in line:
        iteration_obj.n_type = int(line.split()[5])
        preamble_switch[1] = False

    # Extraction of the celldim parameter
    elif 'celldm(1)= ' in line:
        iteration_obj.convert_alat_to_angstrom(float(line.split()[1]))
        preamble_switch[2] = False

    # Extraction of the first lattice vector
    elif 'a(1)' in line:
        x, y, z = map(float, line.split()[3:6])
        a = np.multiply([x, y, z], iteration_obj.alat_to_angstrom)
        iteration_obj.ax = a
        preamble_switch[3] = False

    # Extraction of the second lattice vector
    elif 'a(2)' in line:
        x, y, z = map(float, line.split()[3:6])
        a = np.multiply([x, y, z], iteration_obj.alat_to_angstrom)
        iteration_obj.ay = a
        preamble_switch[4] = False

    # Extraction of the third lattice vector
    elif 'a(3)' in line:
        x, y, z = map(float, line.split()[3:6])
        a = np.multiply([x, y, z], iteration_obj.alat_to_angstrom)
        iteration_obj.az = a
        preamble_switch[5] = False
    
    # Setting of the mass of the atoms inside the groups
    elif 'atomic species   valence' in line:
        for _ in range(iteration_obj.n_type):
            line = fin.readline().split()
            iteration_obj.set_mass(line[0], line[2])
        preamble_switch[6] = False

    # Extraction of the index of atoms and the initial positions
    elif 'site n.' in line:
        for _ in range(iteration_obj.n_atoms):
            line = fin.readline()
            iteration_obj.count_group(line.split())
        preamble_switch[7] = False
        iteration_obj.set_DOF()

    return preamble_switch


def extract_forces(line, iteration_obj):
    """
    Extract forces from the input file and update the forces in the iteration object.

    Parameters
    ----------
    line : str
        The current line from the input file.
    iteration_obj : iteration
        The iteration object to be updated with forces.

    Notes
    -----
    This function reads the input file from the current line until the end of the forces section ('negative rho').
    For each line containing force information, it updates the forces in the iteration object.

    Examples
    --------
    >>> fin = open('input_file.txt', 'r')
    >>> iteration_obj = iteration(groups)
    >>> line = fin.readline()
    >>> extract_forces(line, iteration_obj)
    >>> fin.close()
    """
    while not 'negative rho' in line:
        line = fin.readline()
    
    for _ in range(iteration_obj.n_atoms):
        line = fin.readline().split()
        iteration_obj.forces(line)


def extract_positions(line, iteration_obj, istogram):
    """
    Extract positions from the input file and update the positions in the iteration object and RDF histogram.

    Parameters
    ----------
    line : str
        The current line from the input file.
    iteration_obj : iteration
        The iteration object to be updated with positions.
    istogram : RDF
        The RDF histogram object to be updated with positions.

    Notes
    -----
    This function reads the input file for the specified number of atoms, Extraction of their positions.
    For each line containing position information, it updates the positions in the iteration object and, if applicable, in the RDF histogram.

    Examples
    --------
    >>> fin = open('input_file.txt', 'r')
    >>> iteration_obj = iteration(groups)
    >>> istogram = RDF("filename", 10, ['C', 'H'], 500)
    >>> line = fin.readline()
    >>> extract_positions(line, iteration_obj, istogram)
    >>> fin.close()
    """
    for _ in range(iteration_obj.n_atoms):
        line = fin.readline().split()
        iteration_obj.positions(line, istogram)


def body(fout, iteration_obj, graphs):
    """
    Extract relevant information from the input file and write output to the output file.

    Parameters
    ----------
    fout : file
        The output file object.
    iteration_obj : iteration
        The iteration object containing system information.
    graphs : graph
        The graph object in which the plotting of the graphs happens.

    Notes
    -----
    This function reads the input file line by line, extracting relevant information such as potential energy, forces on atoms, time step, iteration number, and atomic positions.
    It updates the information in the iteration object and writes the corresponding output to the output file.
    Additionally, it calculate the radial distribution function and add the parameter of the timestep in order to draw the graph of energies, temperature and forces in function of time.

    Examples
    --------
    >>> fin = open('input_file.txt', 'r')
    >>> fout = open('output_file.txt', 'w')
    >>> iteration_obj = iteration(groups)
    >>> istogram = RDF("filename", 10, ['C', 'H'], 500)
    >>> body(fout, iteration_obj, istogram)
    >>> fin.close()
    >>> fout.close()
    """
    generation_switch = [True] * 5
    for line in fin:
        
        # Extraction of the potential energy
        if '!    total energy' in line:
            iteration_obj.U_pot = float(line.split()[4]) * 13.605703976
            generation_switch[0] = False

        # Extraction of the forces on atoms
        elif 'Forces a' in line:
            extract_forces(line, iteration_obj)
            generation_switch[1] = False

        # Extraction of the time step
        elif 'Time step' in line:
            iteration_obj.dt = float(line.split()[3])  * 4.8378e-5 #s

        # Extraction of the iteration number
        elif 'Entering' in line:
             iteration_obj.N_iteration = int(line.split()[4])
             generation_switch[2] = False
             
        # Extraction of the atomic position
        elif 'ATOMIC_POSITIONS' in line:
            extract_positions(line, iteration_obj, graphs)
            generation_switch[3] = False
        
        elif 'temperature           =' in line:
            graphs.T = np.append(graphs.T, float(line.split()[2]))
            generation_switch[4] = False
        
        # Generation of the the text printed on the output file and performing the RDF calculation
        if not any(generation_switch):
            fout.writelines(["%s\n" % i for i in iteration_obj.single_frame()])
            graphs.extracting_values(iteration_obj)
            graphs.RDF()
            generation_switch = [True] * 5

def xyz_gen(fout, fin, RDF_, groups):
    """
    Generate an XYZ file from a PWO file, extracting relevant information and updating the output file.

    Parameters
    ----------
    fout : file
        The output XYZ file object.
    fin : file
        The input PWO file object.
    RDF_ : list
        List containing information needed for RDF calculation.
    groups : list
        List of group instances in the system.

    Returns
    -------
    None

    Notes
    -----
    This function checks the input PWO file to extract all the values and generates the corresponding output XYZ file.
    It initializes the iteration object and RDF histogram object, extracts setup information, and then iterates through the file to extract and update information.
    The RDF is calculated, normalized, and plotted at the end within the other graphs.

    Examples
    --------
    >>> fin = open('input_pwo_file.txt', 'r')
    >>> fout = open('output_xyz_file.xyz', 'w')
    >>> RDF_info = ["filename", 10, ['C', 'H'], 500]
    >>> groups = [group_instance_1, group_instance_2]
    >>> xyz_gen(fout, fin, RDF_info, groups)
    >>> fin.close()
    >>> fout.close()
    """
    iteration_obj = iteration(groups)

    # Extraction the setup informations
    preamble_switch = [True] * 8
    while any(preamble_switch):
        preamble_switch = preamble(fin, iteration_obj, preamble_switch)
    
    # Setting of the RDF object
    graphs = graph(RDF_[0], RDF_[1], [RDF_[2], RDF_[3]], RDF_[4])
    graphs.matrix = np.column_stack((iteration_obj.ax, iteration_obj.ay, iteration_obj.az))
    
    # Generation of the configuration at each time
    body(fout, iteration_obj, graphs)

    graphs.normalization(iteration_obj)
    graphs.plot_RDF()
    graphs.plot_energy()
    graphs.plot_forces(iteration_obj)
    graphs.plot_temperature(iteration_obj)

if __name__ == "__main__":
    # Extract setup information
    filename, Rmax, at1, at2, N, groups = setup()
    RDF_ = [filename, Rmax, at1, at2, N]

    # Open output and input files
    with open(filename + '.xyz', "w+") as fout:
        with open(filename + '.pwo', 'r') as fin:
             # Generate XYZ file from PWO file
            xyz_gen(fout, fin, RDF_, groups)