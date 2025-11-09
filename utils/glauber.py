import numpy as np


def sample_woods_saxon(A, R, a, dmin, V0 = 1.0): # we assume V0 = 1
  coords=[]
  while len(coords)<A:
    rmax = 15.0 # should be enough for all nuclei
    r = rmax*np.random.rand() 
    costh = 2*np.random.rand()-1 # [-1,1] for cos(theta) 
    phi = 2*np.pi*np.random.rand()
    prob = V0/(1.0+np.exp((r-R)/a))
    if np.random.rand()<prob:
      x = r*np.sqrt(1-costh**2)*np.cos(phi)
      y = r*np.sqrt(1-costh**2)*np.sin(phi)
      z = r*costh
      ok = True
      for (xx,yy,zz) in coords:
        if (x-xx)**2+(y-yy)**2+(z-zz)**2<dmin**2:
          ok = False
          break
      if ok:
        coords.append((x,y,z))
  return np.array(coords)

def find_collisions(nuc1,nuc2,sig_nn=10.0):
  """ 
  Find collisions between two sets of nucleons.
  nuc1, nuc2: arrays of shape (N,3) with nucleon coordinates
  sig_nn: nucleon-nucleon cross section in fm^2
  Returns:
    part1: boolean array of shape (N1,) indicating which nucleons in nuc1 participated in collisions
    part2: boolean array of shape (N2,) indicating which nucleons in nuc2 participated in collisions
    bc: array of shape (M,3) with coordinates of binary collisions

  """
  rcut = np.sqrt(sig_nn/np.pi)
  part1 = np.zeros(len(nuc1),bool)
  part2 = np.zeros(len(nuc2),bool)
  bc = []
  for i in range(len(nuc1)):
    for j in range(len(nuc2)):
      dx = nuc1[i,0]-nuc2[j,0]
      dy = nuc1[i,1]-nuc2[j,1]
      dz = nuc1[i,2]-nuc2[j,2]
      if dx*dx+dy*dy+dz*dz<rcut*rcut:
        part1[i] = True
        part2[j] = True
        bc.append(((nuc1[i,0]+nuc2[j,0])/2,
                   (nuc1[i,1]+nuc2[j,1])/2,
                   (nuc1[i,2]+nuc2[j,2])/2))
  return part1, part2, np.array(bc)


#Nucleus	Model	<r2>1/2 [fm]	c or a [fm]	z or ? [fm]	w	q-range	reference
# 206Pb	FB	5.49	N/A	N/A	N/A	
class Nuclei:
  def __init__(self, A = 50):
    self.coords = None
    self.A = A # nucleon number
    self.R = 1.25 * A**(1/3) # R = r0 * A^(1/3), r0 = 1.25 fm
    self.a = 0.5 # diffuseness parameter #found on wiki
    self.dmin = 0.4 # minimum inter-nucleon distance

  def sample_ws(self):
    self.coords = sample_woods_saxon(self.A, self.R, self.a, self.dmin)

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
  A = nuc1[wounded1][:, :2] # check this
  B = nuc2[wounded2][:, :2]
  P = np.vstack([A, B])
  if P.shape[0] == 0:
    return 0.0, 0.0

  # recentre
  x = P[:,0] - np.mean(P[:,0])
  y = P[:,1] - np.mean(P[:,1])

  r = np.sqrt(x*x + y*y)
  phi = np.arctan2(y, x)
  w = r**2

  cos2 = np.sum(w * np.cos(2*phi))
  sin2 = np.sum(w * np.sin(2*phi))
  r2   = np.sum(w)

  if r2 == 0:
    return 0.0, 0.0

  eps2 = np.sqrt(cos2**2 + sin2**2) / r2
  psi2 = 1/2.0 * (np.arctan2(sin2, cos2) + np.pi)

  return eps2, psi2







