import re
#da inserire eventuali altri lat par, o volume e finire gli alòtri casi
#fai dizionario
def ibrav(i, a):
    bravais_mapping = {
        0: '  Free',
        1: f'  Cubic (sc): V1 = {a}(1, 0, 0), V2 = {a}(0, 1, 0), V3 = {a}(0, 0, 1)',
        2: '  Cubic (fcc)',
        3: f'  Cubic (bcc): V1 = {a/2}(1, 1, 1), V2 = {a/2}(-1, 1, 1), V3 = {a/2}(-1, -1, 1)',
        -3: '  Cubic (bcc)',
        4: f'  Hexagonal: V1 = {a}(1, 0, 0), V2 = {a}(-1/2, sqrt(3)/2, 0), V3 = (0, 0, c)',
        5: '  Trigonal',
        -5: '  Trigonal',
        6: '  Tetragonal (st)',
        7: '  Tetragonal (bct)',
        8: '  Orthorhombic',
        9: '  Orthorhombic (bco)',
        -9: '  Orthorhombic (bco)',
        91: '  Orthorhombic one-face base-centered',
        10: '  Orthorhombic face-centered',
        11: '  Orthorhombic body-centered',
        12: '  Monoclinic',
        -12: '  Monoclinic',
        13: '  Monoclinic (bc)',
        -13: '  Monoclinic (bc)',
        14: '  Triclinic',
        -14: '  Triclinic',
    }
    return bravais_mapping.get(i, '  stai sbagliando') 
   
def read_info(fin):
    n_at, bravais, lattice_param = 0, 0, 0.0
    for line in fin:
        if 'number of atoms/cell' in line:
            n_at = int(re.search(r'-?\d+(\.\d+)?', line).group())
        elif 'bravais-lattice index' in line:
            bravais = int(re.search(r'-?\d+(\.\d+)?', line).group())
        elif 'lattice parameter' in line:
            lattice_param = float(re.search(r'-?\d+(\.\d+)?', line).group())
    return n_at, bravais, lattice_param

def first_two_lines(fin):
    n_at, brav, lattice_param = read_info(fin)
    return [n_at, ibrav(brav, lattice_param)]

def gen_list_for_xyz(filein):
    with open(filein, 'r') as fin:
        store = fin.read().split('\n')

    fstored = first_two_lines(store)
    t, n_at, ctrl, counter = 0, fstored[0], 0, 0

    for line in store:
        if 'time      =   ' in line:
            t = float(re.search(r'-?\d+(\.\d+)?', line).group())
        elif 'Entering Dynamics:    iteration =' in line:
            N = int(re.search(r'-?\d+(\.\d+)?', line).group())
        elif ctrl > n_at:
            ctrl = 0
        elif ctrl > 0:
            fstored.append(line)
            ctrl += 1
        elif 'ATOMIC_POSITIONS' in line:
            if N == 1:
                fstored[1] = f'- t = {t} ps (iter. {N}) - {fstored[1]}'
            else:
                fstored.append(str(n_at))
                fstored.append(f'- t = {t} ps (iter. {N}) -')
            ctrl = 1
            counter += 1

    return fstored

if __name__ == "__main__":
    file_to_open = str(input('Insert filename: ____.pwo \n'))
    with open(file_to_open + '.xyz', "w+") as fout:
        fout.writelines(["%s\n" % i for i in gen_list_for_xyz(file_to_open + '.pwo')])
