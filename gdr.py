#calcolo gdr
  def calcgdr(self, N ):
    from numpy import sqrt, rint, zeros, int_
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
