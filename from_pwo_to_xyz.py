def gen_list_for_xyz(filein, n_at):

   fin = open(filein, 'r')
   fstored = []
   for line in fin.read().split('\n'):
      fstored.append(line)
   fin.close()

   natomi = n_at
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
   return altratab

file_to_open="Ni_PAO_3Gpa_C1C1_300K_md"


print(file_to_open)

fout=open(str(file_to_open)+"/"+str(file_to_open)+".xyz", "w+")

fout.writelines(["%s\n" % i  for i in gen_list_for_xyz(str(file_to_open)+"/"+str(file_to_open)+".pwo", 212)])
fout.close()
fout.close()
