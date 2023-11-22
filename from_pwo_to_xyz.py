import re

def ibrav(i, a):
   if i == 0 : 
      return 'Free'
   elif i == 1 :
      return 'Cubic (sc)'
   elif i == 2 :
      return 'Cubic (fcc)'
   elif i == 3 :
      return 'Cubic (bcc)'
   elif i == -3:
      return 'Cubic (bcc)'
   elif i == 4 :
      return f'hexagonal: V1 = {a}(1, 0, 0), V2 = {a}(-1/2, sqrt(3)/2, 0), V3 = (0, 0, c)'
   elif i == 5 :
      return 'Trigonal'
   elif i == -5 :
      return 'Trigonal'
   elif i == 6 :
      return 'Tetragonal (st)'
   elif i == 7 :
      return 'Tetragonal (bct)'
   elif i == 8 :
      return 'Orthorhombic'
   elif i == 9 :
      return 'Orthorhombic (bco)'
   elif i == -9 :
      return 'Orthorhombic (bco)'
   elif i == 91 :
      return 'Orthorhombic one-face base-centered'
   elif i == 10 :
      return 'Orthorhombic face-centered'
   elif i == 11 :
      return 'Orthorhombic body-centered' 
   elif i == 12 : 
      return 'Monoclinic'
   elif i == -12 :
      return 'Monoclinic'
   elif i == 13 :
      return 'Monoclinic (bc)'
   elif i == -13 :
      return 'Monoclinic (bc)'
   elif i == 14 :
      return 'Triclinic'
   else :
      return 'stai sbagliando'
   
def read_info(fin):
   for i in fin:
      if 'number of atoms/cell' in i:
         match_n_at = re.search(r'-?\d+(\.\d+)?', i)
         n_at = int(match_n_at.group())
      if 'bravais-lattice index' in i:
         match_bli = re.search(r'-?\d+(\.\d+)?', i)
         bravais = int(match_bli.group())
      if 'lattice parameter' in i:
         match_lp = re.search(r'-?\d+(\.\d+)?', i)
         lattice_param = float(match_lp.group())
   return [n_at, bravais, lattice_param]

def first_two_lines(fin):
   n_at, brav, latt = read_info(fin)
   return [n_at, ibrav(brav, latt)]

def gen_list_for_xyz(filein):
   fin = open(filein, 'r')
   store = []
   for line in fin.read().split('\n'):
      store.append(line)
   fin.close()

   fstored = first_two_lines(store)
   t = 0
   n_at = fstored[0]
   ctrl = 0
   counter = 0
   for i in store:
      if 'time      =   ' in i:
         match_t = re.search(r'-?\d+(\.\d+)?', i)
         t = float(match_t.group())
      if 'Entering Dynamics:    iteration =' in i:
         match_N = re.search(r'-?\d+(\.\d+)?', i)
         N = int(match_N.group())
      if ctrl > n_at:
         ctrl = 0
         fstored.append('\n')
      if ctrl > 0:
         fstored.append(i)
         ctrl += 1
      if 'ATOMIC_POSITIONS' in i:
         fstored.append(n_at)
         if N == 1:
            fstored[1] = f'- t = {t} ps (iter. {N}) -' + fstored[1]
         else:
            fstored.append(f'- t = {t} ps (iter. {N}) -')
         ctrl = 1
         counter += 1
   return(fstored)


#file_to_open=str(input('Insert filename: ____.pwo '))
file_to_open = 'test'

fout=open(file_to_open + '.xyz', "w+")

fout.writelines(["%s\n" % i  for i in gen_list_for_xyz(file_to_open + '.pwo')])
fout.close()