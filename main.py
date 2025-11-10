import sys

import numpy as np
import matplotlib.pyplot as plt

from utils.glauber import sample_woods_saxon, find_collisions, Nuclei, calcPsi2Ecc2

import argparse
from dataclasses import dataclass, field


def plot_res(par, cols, total_ecc2, total_psi2, total_ecc2_std, total_psi2_std):

  f,axs = plt.subplots(1,2, figsize=(14,6))
  
  ax = axs[0]
  ax.errorbar(par.b_array, np.mean(cols, axis=1), yerr=np.std(cols, axis=1), fmt='o-', capsize=5)
  
  ax.set_xlabel("Impact parameter b (fm)", fontsize=14)
  ax.set_ylabel("Number of collisions", fontsize=14)
  
  
  ax = axs[1]
  
  ax.errorbar(par.b_array, total_ecc2, yerr=total_ecc2_std, fmt='o-', label=r'Eccentricity $\epsilon_2$')
  ax.errorbar(par.b_array, total_psi2, yerr=total_psi2_std, fmt='s-', label=r'Participant plane angle $\Psi_2$')
  ax.set_xlabel("Impact parameter b (fm)", fontsize=14)
  ax.set_ylabel("Eccentricity / Participant plane angle", fontsize=14)
  ax.legend(fontsize=12)


  plt.show()


@dataclass
class setup:
  A: int = 206  # nucleon number
  b_array: np.ndarray = np.arange(0,10,1)
  b_array: np.ndarray = field(default_factory=lambda: np.arange(0, 10, 1))

  n_events: int = 1
  save_path: str = "glauber_results.npz"


def main_fortran():
  par = setup()


  cols = []
  total_psi2, total_ecc2 = [], []
  total_psi2_std = []
  total_ecc2_std = []
  for b in par.b_array:
    bcols = []
    psi2, ecc2 = [], []
    for _ in range(par.n_events):
      sys.stdout.write(f"\r  Event {_+1}/{par.n_events} for b={b:.1f} fm")
      sys.stdout.flush()
      nuc1 = Nuclei(par.A)
      nuc1.sample_fws()
      nuc2 = Nuclei(par.A)
      nuc2.sample_fws()
      
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
  
  cols = np.array(cols)
  plot_res(par, cols, total_ecc2, total_psi2, total_ecc2_std, total_psi2_std)
 


def main_py():
  par = setup()
    
  cols = []
  total_psi2, total_ecc2 = [], []
  total_psi2_std = []
  total_ecc2_std = []
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
  
  cols = np.array(cols)
  plot_res(par, cols, total_ecc2, total_psi2, total_ecc2_std, total_psi2_std)
  
if __name__ == "__main__":
  arg_parser = argparse.ArgumentParser()
  arg_parser.add_argument('--v', help='run either python or fortran version', choices=['py', 'fortran'], default='py')

  args = arg_parser.parse_args()
  if args.v == 'py':
    main_py()
  elif args.v == 'fortran':
    main_fortran()

