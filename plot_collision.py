import sys

import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse

from utils.glauber import sample_woods_saxon, find_collisions, Nuclei, calcPsi2Ecc2
from utils import nucphys 
from main import setup

def draw_psi2_plane(ax, psi2, zmin=-25, zmax=25, r=15):
  nx = np.cos(psi2 + np.pi/2)
  ny = np.sin(psi2 + np.pi/2)

  L = r
  corners = np.array([
    [-L*ny,  L*nx, zmin],
    [ L*ny, -L*nx, zmin],
    [ L*ny, -L*nx, zmax],
    [-L*ny,  L*nx, zmax],
  ])
  plane = Poly3DCollection([corners], alpha=0.2, facecolor='tab:purple')
  ax.add_collection3d(plane)
  return plane

def collision_plane():
  np.random.seed(15) # does not matter when using fortran
  #A=206 # nucleon number (
  par = setup()
  b=9.0 
  nuc1 = Nuclei(par)
  nuc1.sample_ws()
  nuc2 = Nuclei(par)
  nuc2.sample_ws()

  #nuc1[:,1]+=b/4
  nuc2[:,1]+=b/2
  nuc1[:,0]-=b/2
  nuc2[:,0]+=b/2

  wounded1 = np.zeros(par.A, 'f', order='F') 
  wounded2 = np.zeros(par.A, 'f', order='F') 
  nucphys.nucutils.find_collisions(nuc1.coords,nuc2.coords, 10.0, wounded1, wounded2, par.A) 
  wounded1 = wounded1.astype(bool)
  wounded2 = wounded2.astype(bool)
  #wounded1, wounded2 = find_collisions(nuc1, nuc2)
  eps2, psi2 = calcPsi2Ecc2(nuc1, nuc2, wounded1, wounded2)

  f = plt.figure(figsize=(6,6))
  plt.scatter(nuc1[wounded1,0], nuc1[wounded1,1], c='green', s=20, label='A wounded')
  plt.scatter(nuc2[wounded2,0], nuc2[wounded2,1], c='orange', s=20, label='B wounded')
  plt.scatter(nuc1[~wounded1,0], nuc1[~wounded1,1], c='blue', s=10, alpha=0.3)
  plt.scatter(nuc2[~wounded2,0], nuc2[~wounded2,1], c='red', s=10, alpha=0.3)

  c1x = np.mean(nuc1[:,0])
  c1y = np.mean(nuc1[:,1])
  c2x = np.mean(nuc2[:,0])
  c2y = np.mean(nuc2[:,1])
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

  sigma = 3. # ellipse size scaling
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
  
  plt.xlim(-15,15)
  plt.ylim(-15,15)

  plt.show()

  f.savefig('figures/collision.png', dpi=300, bbox_inches='tight')

def col_anim():
  #A=206 # nucleon number (
  par = setup()
  b=10.0 # impact parameter (change to centrality at some point)
  z0=20.0 # initial z offset
  nuc1 = Nuclei(par)
  nuc1.sample_fws()
  nuc2 = Nuclei(par)
  nuc2.sample_fws()

  nuc1[:,0]-=b/2
  nuc2[:,0]+=b/2
  wounded1, wounded2 = find_collisions(nuc1, nuc2)
  eps2, psi2 = calcPsi2Ecc2(nuc1, nuc2, wounded1, wounded2)

  def init():
    sc1._offsets3d=([],[],[])
    sc2._offsets3d=([],[],[])
    scb._offsets3d=([],[],[])
    return sc1,sc2,scb

  def update(frame, nuc1, nuc2, collisions, wounded1, wounded2):
    #nonlocal nuc1, nuc2, collisions, wounded1, wounded2
    nuc1[:,2] -= 0.5
    nuc2[:,2] += 0.5
    part1, part2, bc  = find_collisions(nuc1, nuc2, doBC= True)

    p1 = np.asarray(part1)
    p2 = np.asarray(part2)
    if p1.dtype == np.bool_ and p1.shape == (nuc1.shape[0],):
      wounded1 = np.logical_or(wounded1, p1)
    else:
      wounded1[p1] = True
    if p2.dtype == np.bool_ and p2.shape == (nuc2.shape[0],):
      wounded2 = np.logical_or(wounded2, p2)
    else:
      wounded2[p2] = True

    colors1 = np.where(wounded1, 'limegreen', 'darkgreen')
    colors2 = np.where(wounded2, 'red', 'darkred')

    sc1._offsets3d = (nuc1[:,0], nuc1[:,1], nuc1[:,2])
    sc1.set_color(colors1)
    sc2._offsets3d = (nuc2[:,0], nuc2[:,1], nuc2[:,2])
    sc2.set_color(colors2)

    if bc is not None and bc.shape[0] > 0:
      collisions.append(bc)
    if len(collisions) > 0:
      all_bc = np.vstack(collisions)
      scb._offsets3d = (all_bc[:,0], all_bc[:,1], all_bc[:,2])

    return sc1, sc2#, scb
  # draw Psi2 line
  cx = np.mean(np.concatenate([nuc1[wounded1,0], nuc2[wounded2,0]]))
  cy = np.mean(np.concatenate([nuc1[wounded1,1], nuc2[wounded2,1]]))
  L = 10

  nuc1[:,2]+=z0
  nuc2[:,2]-=z0

  fig=plt.figure()
  ax=fig.add_subplot(111,projection='3d')
  ax.set_xlim(-15,15); ax.set_ylim(-15,15); ax.set_zlim(-25,25)
  ax.set_xlabel('x [fm]'); ax.set_ylabel('y [fm]'); ax.set_zlabel('z [fm]')
  #ax.view_init(elev=180, azim=30)
  plane = draw_psi2_plane(ax, psi2)

  sc1=ax.scatter([],[],[],c='darkgreen',s=15)
  sc2=ax.scatter([],[],[],c='darkred',s=15)
  sc_wounded1 = ax.scatter([],[],[],c='lime',s=20)
  sc_wounded2 = ax.scatter([],[],[],c='orangered',s=20)
  scb=ax.scatter([],[],[],c='darkorange',s=5,alpha=0.5)

  
  collisions = []
  wounded1 = np.zeros(len(nuc1), dtype=bool)
  wounded2 = np.zeros(len(nuc2), dtype=bool)

  
  ani=FuncAnimation(fig,update,frames=80,init_func=init,fargs=(nuc1,nuc2,collisions,wounded1,wounded2))
  ani.save("figures/glauber.gif", writer="pillow", fps=20)
  

if __name__ == "__main__":
  arg_parser = argparse.ArgumentParser()
  arg_parser.add_argument('--v', help="Run collision plot or animation", choices=['col', 'anim'], default='col')
  args = arg_parser.parse_args()
  if args.v == 'col':
    collision_plane()
  elif args.v == 'anim':
    col_anim() 
  else:
    print("Unknown option")

