import re
import numpy as np
import matplotlib.pyplot as plt

#generazione lisa con screenshot ai vari t
def gen_list_of_atoms(filein):
    with open(filein, 'r') as fin:
        store = fin.read().split('\n')
    if store[-1] == '':
        store = store[:-1]

    n_at = int(store[0])
    blocks = [store[i + 2: i + 2 + n_at] for i in range(0, len(store), n_at + 2)]    

    return blocks, n_at

# identifica le specie chimiche nel file di input
def atomic_species(block):
    species = []
    indices = []
    for i, line in enumerate(block):
        atom = re.sub(r'[^a-zA-Z]', '', line)
        if atom not in species:
            species.append(atom)
            indices.append([i])
        else:
            j = species.index(atom)
            indices[j].append(i)
    return [species, indices]

# selezione delle specie chimiche da considerare
def setup(block):
    at, ind = atomic_species(block)

    print('The available atomic species are:')
    for i in at:
        print(i)

    print('Choose two of them (including also duplicates if you want)')
    at_1 = str(input())
    at_1 = re.sub(r'[^a-zA-Z]', '', at_1)
    at_2 = str(input())
    at_2 = re.sub(r'[^a-zA-Z]', '', at_2)

    chosen_atoms = [at.index(at_1)] if at_1 == at_2 else [at.index(at_1), at.index(at_2)]
    chosen_indices = [ind[i] for i in chosen_atoms]

    return [chosen_indices, at_1 == at_2]

# estrazione delle coordinate degli atomi
def coordinates(Block, ind):
    positions = [Block[i] for i in ind]
    for idx, i in enumerate(positions):
        matches = re.findall(r'(-?\d+(\.\d+)?)', i)
        positions[idx] = [float(match[0]) for match in matches[:3]]
    return positions

# calcolo le varie distanze di un blocco
def RDF(Block, ind, e, r_max):
    dist = []

    if e == True:  # same species
        pos = np.array(coordinates(Block, ind[0]))
        N = len(pos)
        for k in range(N - 1):
            r = np.linalg.norm(pos[k+1:] - pos[k])
            dist.extend(r[r < r_max])

    else:
        pos_1 = np.array(coordinates(Block, ind[0]))
        pos_2 = np.array(coordinates(Block, ind[1]))

        for i in pos_1:
            r = np.linalg.norm(i - pos_2, axis=1)
            dist.extend(r[r < r_max])

    dist = np.sort(dist)
    return dist

def istogram(blocks, ind, e, r_max, N, n_at):
    count = np.zeros(500)
    R = np.linspace(0, r_max, N)
    dR = R[1]
    norm = np.multiply([((i + dR)**3 - i**3) for i in R], np.pi * 4 /3)
    
    for i in blocks:
        d = RDF(i, ind, e, r_max)
        count += np.histogram(d, bins=N, range=(0, r_max))[0]
    print('distances generated')
    n_at = sum(count)
    V = 4 * np.pi * R[-1]
    rho_1 = V/n_at
    count = np.divide(count, norm) * rho_1

    return [R, count]


def main():
    file_to_open = str(input('Insert filename: ____.xyz \n'))
    blocks, n_at = gen_list_of_atoms(file_to_open + '.xyz')
    ind, e = setup(blocks[0])
    r_max = float(input('Insert Rmax:\n'))
    data = istogram(blocks, ind, e, r_max, 500, n_at)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(data[0], data[1], label='Power spectrum')
    ax.set_xlabel('r ($A$)')
    ax.set_ylabel('g(r)')
    ax.grid()
    plt.savefig(f'{file_to_open}.png')
    plt.show()

if __name__ == "__main__":
    main()
