#calcolo gdr
from numpy import sqrt, rint, zeros, int_
import re
import numpy as np



def calcgdr(self, N ):
    for k in range(N-1) :
        j=k+1
        #for j in range(k+1,N) :
        #  rx,ry,rz  reduced coordinates  x/L,y/L,z/L  , L = length of cubic box
        dx = self.rx[k]-self.rx[j:N]
        dy = self.ry[k]-self.ry[j:N]
        dz = self.rz[k]-self.rz[j:N]
        # Periodic boundary conditions    
        dx[...]-= rint(dx)
        dy[...]-= rint(dy)
        dz[...]-= rint(dz)
        dx[...] = dx*self.L
        dy[...] = dy*self.L
        dz[...] = dz*self.L
        r2 = dx*dx + dy*dy + dz*dz
        # using the mask array "b" for speedup
        b = r2 < self.r2max
        lm  = sqrt(r2[b])
        #if lm<self.kg :
        for elm in lm :
            self.gcount[int(elm/self.ldel)]+=2.  # factor of 2 for gdr normalization
    #return




# normalizzazione - gdr 
def write_gdr(self, N, T, rho, gdr_out='gdr.out'):
      """ here L and rho from time averages ? """
      from numpy import zeros, pi, savetxt, column_stack
      V = zeros(self.kg) 
      r = zeros(self.kg)
      g = zeros(self.kg) 
      for lm in range(self.kg) :
          V[lm] = 4./3.*pi*(self.ldel**3)*(3*lm*lm +3*lm + 1); 
          g[lm] = self.gcount[lm]/(V[lm]*(N -1)*T*rho);
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
   for i in range(np.size(block)):
      atom = re.sub(r'[^a-zA-Z]', '', block[i])
      if atom not in species:
         species.append(atom)
         indices.append([i+1])
      else:
         j = species.index(atom)
         indices[j].append(i+1)
   return [species, indices]
         
    

#file_to_open=str(input('Insert filename: ____.xyz \n'))
file_to_open = 'CH4'

Blocks = gen_list_of_atoms(file_to_open+'.xyz')
print(atomic_species(Blocks[0]))