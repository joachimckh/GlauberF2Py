import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse

from utils.glauber import sample_woods_saxon, find_collisions, Nuclei, calcPsi2Ecc2

def main():
  np.random.seed(11)
  A=206 # nucleon number (
  b=9.0 # impact parameter (change to centrality at some point)
  nuc1 = Nuclei(A)
  nuc1.sample_fws()
  nuc2 = Nuclei(A)
  nuc2.sample_fws()

  #nuc1[:,1]+=b/4
  nuc2[:,1]+=b/2
  nuc1[:,0]-=b/2
  nuc2[:,0]+=b/2

  wounded1, wounded2, _ = find_collisions(nuc1, nuc2)
  eps2, psi2 = calcPsi2Ecc2(nuc1, nuc2, wounded1, wounded2)

  f = plt.figure(figsize=(6,6))
  plt.scatter(nuc1[wounded1,0], nuc1[wounded1,1], c='green', s=20, label='A wounded')
  plt.scatter(nuc2[wounded2,0], nuc2[wounded2,1], c='orange', s=20, label='B wounded')
  plt.scatter(nuc1[~wounded1,0], nuc1[~wounded1,1], c='blue', s=10, alpha=0.3)
  plt.scatter(nuc2[~wounded2,0], nuc2[~wounded2,1], c='red', s=10, alpha=0.3)

  c1x = np.mean(nuc1[:,0]); c1y = np.mean(nuc1[:,1])
  c2x = np.mean(nuc2[:,0]); c2y = np.mean(nuc2[:,1])
  bx = c2x - c1x
  by = c2y - c1y
  psi_rp = np.arctan2(by, bx)

  cx = np.mean(np.concatenate([nuc1[wounded1,0], nuc2[wounded2,0]]))
  cy = np.mean(np.concatenate([nuc1[wounded1,1], nuc2[wounded2,1]]))
  L = 10

  plt.plot([cx - L*np.cos(psi_rp), cx + L*np.cos(psi_rp)],
           [cy - L*np.sin(psi_rp), cy + L*np.sin(psi_rp)],
           color='tab:purple', lw=2, label=r'$\Psi^\mathrm{RP}$')

  plt.plot([cx - L*np.cos(psi2), cx + L*np.cos(psi2)],
           [cy - L*np.sin(psi2), cy + L*np.sin(psi2)],
           color='black', lw=2, label=r'$\Psi_2$')



  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  plt.axis('equal')
  plt.legend(frameon=False)
  plt.title(fr'$\varepsilon_2={eps2:.3f}$,  $\Psi_2={(psi2):.3f}$')


  x_w = np.concatenate([nuc1[wounded1,0], nuc2[wounded2,0]])
  y_w = np.concatenate([nuc1[wounded1,1], nuc2[wounded2,1]])

  sigma = 3.
  sx = np.std(x_w) * sigma
  sy = np.std(y_w) * sigma

  ellipse = Ellipse((cx, cy),
                    width=sx, height=sy,
                    angle=np.degrees(psi2),
                    edgecolor='k', facecolor='none', lw=2, ls='--')

  plt.gca().add_patch(ellipse)

  r_nucleus = nuc1.R
  circle1 = Circle((c1x, c1y), r_nucleus, edgecolor='blue',alpha=0.5, facecolor='none', lw=1., ls='-')
  circle2 = Circle((c2x, c2y), r_nucleus, edgecolor='red',alpha=0.5, facecolor='none', lw=1., ls='-')
  plt.gca().add_patch(circle1)
  plt.gca().add_patch(circle2)


  plt.show()
  f.savefig('figures/one_glauber_col_with_RP_Psi2.png', dpi=300, bbox_inches='tight')

if __name__ == "__main__":
  main()

