import sys
import numpy as np
import time
import argparse
from dataclasses import dataclass, field

from utils.plotting import plot_res
from utils.glauber import sample_woods_saxon, find_collisions, Nuclei, calcPsi2Ecc2
from utils import nucphys 


@dataclass
class setup:
  A: int = 206  # nucleon number
  b_array: np.ndarray = field(default_factory=lambda: np.arange(0, 10, 1))

  n_events: int = 10
  # save_path: str = "glauber_results.npz"

def main_purefort():
  par = setup()
  #print(nucphys.simulation.__doc__)

  cols = np.zeros(par.b_array.shape, 'f', order='F')
  ecc2 = np.zeros(par.b_array.shape, 'f', order='F')
  psi2 = np.zeros(par.b_array.shape, 'f', order='F')
  nucphys.simulation.proc(par.n_events, par.b_array, par.A, cols, psi2, ecc2)

  print(cols)
  # plot_res(par, cols, ecc2, psi2, np.zeros_like(ecc2), np.zeros_like(psi2))




def main_fortran():
  # print(nucphys.nucutils.__doc__)
  par = setup()


  cols = []
  total_psi2, total_ecc2 = [], []
  total_psi2_std = []
  total_ecc2_std = []

  time0 = time.time()
  for b in par.b_array:
    bcols = []
    psi2, ecc2 = [], []
    for _ in range(par.n_events):
      sys.stdout.write(f"\r  Event {_+1}/{par.n_events} for b={b:.1f} fm")
      sys.stdout.flush()
      nuc1 = Nuclei(par.A)
      nuc1.sample_ws()
      nuc2 = Nuclei(par.A)
      nuc2.sample_ws()
      
      nuc1[:,0]-=b/2
      nuc2[:,0]+=b/2

      wounded1 = np.zeros(par.A, 'f', order='F') 
      wounded2 = np.zeros(par.A, 'f', order='F') 
      nucphys.nucutils.find_collisions(nuc1.coords,nuc2.coords, 10.0, wounded1, wounded2, par.A) 
      wounded1 = wounded1.astype(bool)
      wounded2 = wounded2.astype(bool)


      #wounded1, wounded2, bc = find_collisions(nuc1, nuc2)
      bcols.append(np.count_nonzero(wounded1) + np.count_nonzero(wounded2))
      psi, ecc = calcPsi2Ecc2(nuc1, nuc2, wounded1, wounded2)
      psi2.append(psi)
      ecc2.append(ecc)
    cols.append(bcols)
    total_psi2.append(np.mean(psi2))
    total_ecc2.append(np.mean(ecc2))
    total_psi2_std.append(np.std(psi2))
    total_ecc2_std.append(np.std(ecc2))

  time1 = time.time()
  print(f"\nFortran version took {time1 - time0:.2f} seconds")
  
  cols = np.array(cols)
  # print(cols.mean(axis=1))
  plot_res(par, cols, total_ecc2, total_psi2, total_ecc2_std, total_psi2_std)
 


def main_py():
  par = setup()
    
  cols = []
  total_psi2, total_ecc2 = [], []
  total_psi2_std = []
  total_ecc2_std = []

  time0 = time.time()
  for b in par.b_array:
    bcols = []
    psi2, ecc2 = [], []
    for _ in range(par.n_events):
      sys.stdout.write(f"\r  Event {_+1}/{par.n_events} for b={b:.1f} fm")
      sys.stdout.flush()
      nuc1 = Nuclei(par.A)
      nuc1.sample_ws()
      nuc2 = Nuclei(par.A)
      nuc2.sample_ws()
      
      nuc1[:,0]-=b/2
      nuc2[:,0]+=b/2
      wounded1, wounded2, bc = find_collisions(nuc1, nuc2)
      bcols.append(np.count_nonzero(wounded1) + np.count_nonzero(wounded2))
      psi, ecc = calcPsi2Ecc2(nuc1, nuc2, wounded1, wounded2)
      psi2.append(psi)
      ecc2.append(ecc)
    cols.append(bcols)
    total_psi2.append(np.mean(psi2))
    total_ecc2.append(np.mean(ecc2))
    total_psi2_std.append(np.std(psi2))
    total_ecc2_std.append(np.std(ecc2))

  time1 = time.time()
  print(f"\nPython version took {time1 - time0:.2f} seconds")
  
  cols = np.array(cols)
  plot_res(par, cols, total_ecc2, total_psi2, total_ecc2_std, total_psi2_std)
  
if __name__ == "__main__":
  arg_parser = argparse.ArgumentParser()
  arg_parser.add_argument('--v', help='run either python or fortran version', choices=['py', 'f', 'pf'], default='py')

  args = arg_parser.parse_args()
  if args.v == 'py':
    main_py()
  elif args.v == 'f':
    main_fortran()
  elif args.v == 'pf':
    main_purefort()
