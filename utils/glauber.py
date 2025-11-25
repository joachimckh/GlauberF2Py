import numpy as np
from utils import nucphys 

#params au:
# w = 0 # true for gold, I dont really use this, but keep for completeness
def rho(r, R, a, w=0):
  rho0 = 1 # nuclear density in the center (just use 1) 
  num = 1 + w #w = w(r/R)=0 for gold (in general spherical nuclei) -> then rho decays to woods_saxon
  denum = 1 + np.exp((r - R) / a)
  return rho0 * num/denum

def sample_f(r, R, a):
  return 4*np.pi*r**2 * rho(r, R, a)

def sample_woods_saxon(A, R, a, dmin):
  """
    Rejection sampling of nucleon coordinates from Woods-Saxon distribution
  """
  rmax = 15.0 # same as TGlauber

  # r_grid = np.linspace(0.0, rmax, 1000)
  # fmax = np.max(sample_f(r_grid, R, a))
  fmax = 400.0 # precomputed maximum of sample_f for R=6.38, a=0.535, rmax=15.0 - should find a better way to do this

  coords = []
  while len(coords) < A:
    r = rmax * np.random.rand()
    costh = 2*np.random.rand() - 1 # [-1,1] for cos
    phi = 2*np.pi * np.random.rand()
    u = np.random.rand() * fmax
    if u < sample_f(r, R, a):
      s = np.sqrt(1 - costh*costh) # sin(theta)
      x = r * s * np.cos(phi)
      y = r * s * np.sin(phi)
      z = r * costh
      ok = True
      for (xx, yy, zz) in coords:
        if (x-xx)**2 + (y-yy)**2 + (z-zz)**2 < dmin**2:
          ok = False
          break
      if ok:
        coords.append((x, y, z))
  return np.array(coords)

def find_collisions(nuc1,nuc2,sig_nn=10.0, doBC=False):
  """ 
  Find collisions between two sets of nucleons.
  nuc1, nuc2: arrays of shape (N,3) with nucleon coordinates
  sig_nn: nucleon-nucleon cross section in fm^2
  Returns:
    part1: boolean array of shape (N1,) indicating which nucleons in nuc1 participated in collisions
    part2: boolean array of shape (N2,) indicating which nucleons in nuc2 participated in collisions
  """
  rcut = np.sqrt(sig_nn/np.pi)
  part1 = np.zeros(len(nuc1),bool)
  part2 = np.zeros(len(nuc2),bool)
  if doBC:
    bc = []
  for i in range(len(nuc1)):
    for j in range(len(nuc2)):
      dx = nuc1[i,0]-nuc2[j,0]
      dy = nuc1[i,1]-nuc2[j,1]
      dz = nuc1[i,2]-nuc2[j,2]
      if dx*dx+dy*dy+dz*dz<rcut*rcut:
        part1[i] = True
        part2[j] = True
        if doBC:
          bc.append(((nuc1[i,0]+nuc2[j,0])/2,
                   (nuc1[i,1]+nuc2[j,1])/2,
                   (nuc1[i,2]+nuc2[j,2])/2))
  if doBC:
    return part1, part2, np.array(bc)
  return part1, part2



#Nucleus	Model	<r2>1/2 [fm]	c or a [fm]	z or ? [fm]	w	q-range	reference
# 206Pb	FB	5.49	N/A	N/A	N/A	
class Nuclei:
  def __init__(self, par):
    self.coords = None
    self.A = par.A # nucleon number
    self.R = par.R #1.25 * A**(1/3) # R = r0 * A^(1/3), r0 = 1.25 fm
    self.a = par.a #0.5 # diffuseness parameter #found on wiki
    self.dmin = par.dmin # 0.4 # minimum inter-nucleon distance

  def sample_ws(self):
    self.coords = sample_woods_saxon(self.A, self.R, self.a, self.dmin)
    tmp = np.array(self.coords)

  def sample_fws(self):
    self.coords = np.zeros((self.A, 3), 'f', order='F')
    #fwoodsax.nucmath.woods_saxon(self.R, self.a, self.dmin, self.coords)
    nucphys.nucmath.woods_saxon(self.R, self.a, self.dmin, self.coords, self.A)

  def __getitem__(self, index):
    return self.coords[index]

  def __setitem__(self, index, value):
    self.coords[index] = value

  def __len__(self):
    return len(self.coords)

  @property
  def shape(self):
    return self.coords.shape

def calcPsi2Ecc2(nuc1, nuc2, wounded1, wounded2):
  # collect wounded nucleons from both nuclei
  A = nuc1[wounded1][:, :2] 
  B = nuc2[wounded2][:, :2]
  P = np.vstack([A, B])
  if P.shape[0] == 0:
    return -999.0, -999.0

  x = P[:,0] - np.mean(P[:,0])
  y = P[:,1] - np.mean(P[:,1])

  r = np.sqrt(x*x + y*y)
  phi = np.arctan2(y, x)
  w = r**2

  cos2 = np.sum(w * np.cos(2*phi))
  sin2 = np.sum(w * np.sin(2*phi))
  r2 = np.sum(w)

  if r2 == 0:
    return -999, -999 

  eps2 = np.sqrt(cos2**2 + sin2**2) / r2
  psi2 = 1/2.0 * (np.arctan2(sin2,cos2) + np.pi) # shift from -pi,pi to 0,2pi

  return eps2, psi2

if __name__ == "__main__":
  print(nucphys.nucmath.woods_saxon.__doc__)




