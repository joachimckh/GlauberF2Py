import matplotlib.pyplot as plt
import numpy as np

def plot_res(par, cols, total_ecc2 = None, total_psi2 = None, ftype = "python", total_time = 0):

  f,axs = plt.subplots(1,3, figsize=(26,6))
  
  ax = axs[0]
  for i, b in enumerate(par.b_array):
    x = np.full(cols.shape[1], b)
    y = cols[i]
    ax.scatter(x, y, s=1, color="gray")
  
  ax.set_ylabel("Number of collisions", fontsize=14)
  
  total_events = par.n_events * len(par.b_array)
  ax.annotate("Algorithm: "+ftype, xy=(0.5,0.9), xycoords="axes fraction")
  ax.annotate("Total events: "+str(total_events), xy=(0.5,0.85), xycoords="axes fraction")
  ax.annotate(f"Total time: {total_time:.2f} s", xy=(0.5,0.8), xycoords="axes fraction")
  ax.annotate(f"Time/event: {total_time/total_events*1000:.2f} ms", xy=(0.5,0.75), xycoords="axes fraction")
  
  ax = axs[1]
  ax2 = axs[2]
  
  if total_ecc2 is not None and total_psi2 is not None:
    total_ecc2 = np.array(total_ecc2)
    total_psi2 = np.array(total_psi2)
    for i, b in enumerate(par.b_array):
      x = np.full(total_ecc2.shape[1], b)
      y = total_ecc2[i]
      ax.scatter(x, y, s=1, color="limegreen", alpha=1)
      y2 = total_psi2[i]
      ax2.scatter(x, y2, s=1, color="darkred", alpha=1)

    ax.set_ylabel(r"$\varepsilon_2$", fontsize=14)
    ax2.set_ylabel(r"$\Psi_2$", fontsize=14)
    

    ax.set_ylim(-0.5,3.5)
    ax2.set_ylim(-0.5,3.5)

  for ax in axs:
    ax.set_xlabel("Impact parameter b (fm)", fontsize=14)
    ax.set_xlim(0, par.b_array.max())

  plt.show()

  f.savefig(f"figures/nc_impb_alg{ftype}.png", bbox_inches='tight')
