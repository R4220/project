#calcolo gdr
import re
import numpy as np
from numpy import pi, sqrt, rint, zeros, int_
import matplotlib.pyplot as plt

# normalizzazione - gdr 
def write_gdr(self, N, T, rho, gdr_out='gdr.out'):
      """ here L and rho from time averages ? """
      from numpy import zeros, pi, savetxt, column_stack
      V = zeros(self.kg) 
      r = zeros(self.kg)
      g = zeros(self.kg) 
      for lm in range(self.kg) :
          V[lm] = 4./3.*pi*(self.ldel**3)*(3*lm*lm +3*lm + 1); 
          g[lm] = self.gcount[lm]/(V[lm]*(N -1)*T*rho)
          r[lm] = (lm+0.5)*self.ldel
      gout = column_stack( (r, g) )
      savetxt(gdr_out, gout , fmt=('%12.7g ','%12.7g'), header="    'r'     'g(r)'" ) 


#generazione lisa con screenshot ai vari t
def gen_list_of_atoms(filein):
   fin = open(filein, 'r')
   store = []
   for line in fin.read().split('\n'):
      store.append(line)
   fin.close()
   n_at = int(store[0])
   if store[-1] == '':
      store = store[:-1]
   N_iter = int(np.size(store)/(n_at + 2))
   # aggiungi assert che N_iter deve essere intero
   blocks = []
   for i in range(N_iter):
      j = int(i * (n_at + 2))
      blocks.append(store[j + 2: j + 2 + n_at])
   return blocks

# identifica le specie chimiche nel file di input
def atomic_species(block):
   species = []
   indices = []
   for i in range(np.size(block)):#enumerate
      atom = re.sub(r'[^a-zA-Z]', '', block[i])
      if atom not in species:
         species.append(atom)
         indices.append([i])
      else:
         j = species.index(atom)
         indices[j].append(i)
   return [species, indices]

# selezione delle specie chimiche da considerare
def selection_of_couples(atoms):
   print('The aviable atomic species are:')
   for i in atoms:
      print(i)
   print('Choose two of them (including also duplicates if you want)')
   at_1 = str(input())
   at_1 = re.sub(r'[^a-zA-Z]', '', at_1)
   at_2 = str(input())  
   at_2 = re.sub(r'[^a-zA-Z]', '', at_2) 
   if at_1 == at_2:
      return[[atoms.index(at_1)], 0]
   else:
      return [[atoms.index(at_1), atoms.index(at_2)], 1]

# definizione del setup
def setup(Block):
   at, ind = atomic_species(Block)
   chosen_atoms, e = selection_of_couples(at)
   #print(ind)
   #print(chosen_atoms)
   chosen_indices = [ind[i] for i in chosen_atoms]
   #print(chosen_indices)
   return [chosen_indices, e]

# estrazione delle coordinate degli atomi
def coordinates(Block, ind):
    positions = [Block[i] for i in ind]
    for idx, i in enumerate(positions):
        matches = re.findall(r'(-?\d+(\.\d+)?)', i)
        #print(matches)
        positions[idx] = [float(match[0]) for match in matches[:3]]
    #print(positions)
    return positions

# calcolo le varie distanze di un blocco
def RDF(Block, ind, e, r_max):
   dist = []
   if e == 0 : # same species
      ind = ind[0]
      #print(Block)
      pos = np.array(coordinates(Block, ind))
      N = len(pos)
      #print(pos)
      for k in range(N-1) :
         for j in range(k + 1, N):
            #print(k, '  ', j)
            r = np.linalg.norm(pos[k]-pos[j])
            if r < r_max:
               dist.append(r)  
               #print(pos[k], '\n', pos[j], '\n', pos[k]-pos[j], '\n', np.linalg.norm(pos[k]-pos[j]), '\n')
   else:
      ind_1 = ind[0]
      ind_2 = ind[1]
      pos_1 = np.array(coordinates(Block, ind_1))
      pos_2 = np.array(coordinates(Block, ind_2))
      for i in pos_1:
         for j in pos_2:
            r = np.linalg.norm(i - j)
            if r < r_max:
               dist.append(r)
   dist = np.sort(dist)
   return dist

def istogram(Blocks, ind, e, r_max, N):
   count = np.zeros(500)
   for i in Blocks:
      d = RDF(i, ind, e, r_max)
      count += np.histogram(d, bins=N, range=(0, r_max))[0]
   norm, R = normalize(N, r_max)
   count = np.divide(count, norm )
   return [R, count]

def normalize(N_bin, r_max):
   R = np.linspace(0, r_max, N_bin)
  # print(R)
   dR = R[1]
#   print(dR)
 #  print([(i + dR)**3 - i**3 for i in R])
   norm = np.multiply([(i + dR)**3 - i**3 for i in R], pi * 4 /3)
   return [norm, R]


file_to_open=str(input('Insert filename: ____.xyz \n'))
#file_to_open = 'test'

Blocks = gen_list_of_atoms(file_to_open+'.xyz')
ind, e = setup(Blocks[0])
r_max = float(input('Instert Rmax:\n')) # variabile globale?
data = istogram(Blocks, ind, e, r_max, 500)

fig = plt.figure(figsize=(6, 8))
ax = fig.add_subplot(1, 1, 1)
ax.plot(data[0], data[1], label = 'Power spectrum')
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Power density (V$^2$/Hz)')
ax.grid()
plt.show()