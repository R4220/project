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
   for line in fin.read().split('\n'):
      if 'number of atoms/cell' in line:
         match_n_at = re.search(r'-?\d+(\.\d+)?', line)
         n_at = int(match_n_at.group())
      if 'bravais-lattice index' in line:
         match_bli = re.search(r'-?\d+(\.\d+)?', line)
         bravais = int(match_bli.group())
      if 'lattice parameter' in line:
         match_lp = re.search(r'-?\d+(\.\d+)?', line)
         lattice_param = float(match_lp.group())
   return [n_at, bravais, lattice_param]

def first_two_lines(fin):
   n_at, brav, latt = read_info(fin)
   return [n_at, ibrav(brav, latt)]

def gen_list_for_xyz(filein):
   fin = open(filein, 'r')
   fstored = first_two_lines(fin)
   print(fstored)
   fin.close()

   '''natomi = n_at
   altratab = []
   ctrl = 0
   counter = 0
   for line in fstored:
      if ctrl == natomi+1:
         ctrl = 0
      if ctrl>0:
         altratab.append(line)#.split('\n')
         ctrl += 1
      if 'ATOMIC_POSITION' in line:
         #altratab.append(line)
         ctrl = 1
         counter += 1
         altratab.append(natomi)
         altratab.append(float(counter))

   altratab.append(' ')
   return altratab'''
   return


#file_to_open=input('Insert filename: ')
file_to_open = 'test.pwo'
gen_list_for_xyz(file_to_open)

'''fout=open(str(file_to_open), "w+")

fout.writelines(["%s\n" % i  for i in gen_list_for_xyz(str(file_to_open)+"/"+str(file_to_open)+".pwo", 212)])
fout.close()
fout.close()'''