#calcolo gdr
from numpy import sqrt, rint, zeros, int_
import re
import numpy as np

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

def selection_of_couples(atoms):
   print('The aviable atomic species are:', atoms)
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

def setup(Block):
   at, ind = atomic_species(Block)
   chosen_atoms, e = selection_of_couples(at)
   #print(ind)
   #print(chosen_atoms)
   chosen_indices = [ind[i] for i in chosen_atoms]
   #print(chosen_indices)
   return [chosen_indices, e]

def coordinates(Block, ind):
    positions = [Block[i] for i in ind]
    for idx, i in enumerate(positions):
        matches = re.findall(r'-?\d+(\.\d+)?', i)
        positions[idx] = [float(match) for match in matches[:3]]
    #print(positions)
    return positions

def RDF(Block, ind, e, r_max):
   dist = []
   if e == 0 : # same species
      ind = ind[0]
      pos = np.array(coordinates(Block, ind))
      N = len(pos)
      print(N)
      for k in range(N-1) :
         for j in range(k + 1, N):
            print(k, '  ', j)
            r = np.linalg.norm(pos[k]-pos[j])
            if r < r_max:
               dist.append(r)  
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
   return dist

file_to_open=str(input('Insert filename: ____.xyz \n'))
#file_to_open = 'test'

Blocks = gen_list_of_atoms(file_to_open+'.xyz')
ind, e = setup(Blocks[0])
r_max = float(input('Instert Rmax:\n'))
print(RDF(Blocks[0], ind, e, r_max))