# main

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
    iter = iteration()
    e_storage = False
    e_gen = False
    
    for idx,line in enumerate(fin):
        
        if 'number of atoms/cell' in line:
            iter.n_at = int(line.split()[4])
            print(iter.n_at)
        elif 'a(1)' in line:
            x, y, z = map(float, line.split()[3:6])
            iter.ax = [x, y, z]
            #print(iter.ax)
        elif 'a(2)' in line:
            x, y, z = map(float, line.split()[3:6])
            iter.ay = [x, y, z]
            #print(iter.ay)
        elif 'a(3)' in line:
            x, y, z = map(float, line.split()[3:6])
            iter.az = [x, y, z]
            #print(iter.az)

        if 'Self-consistent Calculation' in line:
            fout.writelines(iter.first_two_lines())
            if e_gen == True:
                #calculation
            e_gen = True
            e_storage = True
        
        if e_storage == True:
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