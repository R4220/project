import re
#da inserire eventuali altri lat par, o volume e finire gli alòtri casi

conversion = 0.529177249

def ibrav_diz(i, a):
    a = a * conversion
    bravais_mapping = {
        0: '  Free',
        1: f'  Cubic (sc): V1 = {a}(1, 0, 0) $\AA$, V2 = {a}(0, 1, 0) $\AA$, V3 = {a}(0, 0, 1) $\AA$',
        2: '  Cubic (fcc)',
        3: f'  Cubic (bcc): V1 = {a/2}(1, 1, 1) $\AA$, V2 = {a/2}(-1, 1, 1) $\AA$, V3 = {a/2}(-1, -1, 1) $\AA$',
        -3: '  Cubic (bcc)',
        4: f'  Hexagonal: V1 = {a}(1, 0, 0) $\AA$, V2 = {a}(-1/2, sqrt(3)/2, 0) $\AA$, V3 = (0, 0, c) $\AA$',
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

def ibrav(i, a, V):
    a = a * conversion
    V = V * conversion**3
    if i == 0 : 
        return '  Free'  
    elif i == 1 :
        return f'  Cubic (sc): V1 = {a}(1, 0, 0) Ang, V2 = {a}(0, 1, 0) Ang, V3 = {a}(0, 0, 1) Ang'  
    elif i == 2 :
        return '  Cubic (fcc)'  
    elif i == 3 :
        return f'  Cubic (bcc): V1 = {a/2}(1, 1, 1) Ang, V2 = {a/2}(-1, 1, 1) Ang, V3 = {a/2}(-1, -1, 1) Ang'  
    elif i == -3:
        return '  Cubic (bcc)'  
    elif i == 4 :
        c = V *2 / (a**2 * (sqrt(3)))
        return f'  Hexagonal: V1 = {a}(1, 0, 0) Ang, V2 = {a}(-1/2, sqrt(3)/2, 0) Ang, V3 = {c}(0, 0, 1) Ang'  
    elif i == 5 :
        return '  Trigonal'  
    elif i == -5 :
        return '  Trigonal'  
    elif i == 6 :
        return '  Tetragonal (st)'  
    elif i == 7 :
        return '  Tetragonal (bct)'  
    elif i == 8 :
        return '  Orthorhombic'  
    elif i == 9 :
        return '  Orthorhombic (bco)'  
    elif i == -9 :
        return '  Orthorhombic (bco)'  
    elif i == 91 :
        return '  Orthorhombic one-face base-centered'  
    elif i == 10 :
        return '  Orthorhombic face-centered'  
    elif i == 11 :
        return '  Orthorhombic body-centered'   
    elif i == 12 : 
        return '  Monoclinic'  
    elif i == -12 :
        return '  Monoclinic'  
    elif i == 13 :
        return '  Monoclinic (bc)'  
    elif i == -13 :
        return '  Monoclinic (bc)'  
    elif i == 14 :
        return '  Triclinic'  
    else :
        return '  stai sbagliando'  
 
def read_info(fin):
    n_at, bravais, lattice_param, V = 0, 0, 0.0
    for line in fin:
        if 'number of atoms/cell' in line:
            n_at = int(re.search(r'-?\d+(\.\d+)?', line).group())
        elif 'bravais-lattice index' in line:
            bravais = int(re.search(r'-?\d+(\.\d+)?', line).group())
        elif 'lattice parameter' in line:
            lattice_param = float(re.search(r'-?\d+(\.\d+)?', line).group())
        elif 'unit-cell volume' in line:
            V = float(re.search(r'-?\d+(\.\d+)?', line).group())
    return n_at, bravais, lattice_param, V

def first_two_lines(fin):
    n_at, brav, lattice_param = read_info(fin)
    return [n_at, ibrav(brav, lattice_param, V)]

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
