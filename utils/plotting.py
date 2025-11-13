import matplotlib.pyplot as plt
import numpy as np

def plot_res(par, cols, total_ecc2 = None, total_psi2 = None, ftype = "python"):

  f,axs = plt.subplots(1,2, figsize=(14,6))
  
  ax = axs[0]
  #ax.errorbar(par.b_array, np.mean(cols, axis=1), yerr=np.std(cols, axis=1), fmt='o-', capsize=5)
  for i, b in enumerate(par.b_array):
    x = np.full(cols.shape[1], b)
    y = cols[i]
    ax.scatter(x, y, s=1, color="gray")
  
  ax.set_xlabel("Impact parameter b (fm)", fontsize=14)
  ax.set_ylabel("Number of collisions", fontsize=14)
  
  
  ax = axs[1]
  
  if total_ecc2 is not None and total_psi2 is not None:
    ax.errorbar(par.b_array, np.mean(total_ecc2,axis=1), yerr=np.std(total_ecc2,axis=1), fmt='o-', label=r'Eccentricity $\epsilon_2$')
    ax.errorbar(par.b_array, np.mean(total_psi2, axis=1), yerr=np.std(total_psi2,axis=1), fmt='s-', label=r'Participant plane angle $\Psi_2$')
    ax.set_xlabel("Impact parameter b (fm)", fontsize=14)
    ax.set_ylabel("Eccentricity / Participant plane angle", fontsize=14)
    ax.legend(fontsize=12)

    ax.set_xlim(0, par.b_array.max())

  plt.show()

  f.savefig(f"figures/nc_impb_alg{ftype}.png", bbox_inches='tight')
