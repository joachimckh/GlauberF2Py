import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from utils.glauber import sample_woods_saxon, find_collisions, Nuclei, calcPsi2Ecc2


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

A=206 # nucleon number (
b=10.0 # impact parameter (change to centrality at some point)
z0=20.0 # initial z offset
nuc1 = Nuclei(A)
nuc1.sample_ws()
nuc2 = Nuclei(A)
nuc2.sample_ws()

nuc1[:,0]-=b/2
nuc2[:,0]+=b/2
wounded1, wounded2, _ = find_collisions(nuc1, nuc2)
eps2, psi2 = calcPsi2Ecc2(nuc1, nuc2, wounded1, wounded2)


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

def init():
  sc1._offsets3d=([],[],[])
  sc2._offsets3d=([],[],[])
  scb._offsets3d=([],[],[])
  return sc1,sc2,scb

collisions = []
wounded1 = np.zeros(len(nuc1), dtype=bool)
wounded2 = np.zeros(len(nuc2), dtype=bool)

def update(frame):
  global nuc1, nuc2, collisions, wounded1, wounded2
  nuc1[:,2] -= 0.5
  nuc2[:,2] += 0.5
  part1, part2, bc = find_collisions(nuc1, nuc2)

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

  return sc1, sc2, scb

ani=FuncAnimation(fig,update,frames=80,init_func=init,interval=100,blit=False)
ani.save("figures/glauber.gif", writer="pillow", fps=20)









