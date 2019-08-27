import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def plot_histogram_overlap(bins, hists, args, kT):
  ## check for overlap so that iteration doesn't diverge
  fig, ax = plt.subplots()
  lines = ax.plot(bins, hists)
  ax.set_xlabel(r'$N$')
  ax.set_ylabel(r'$P(N)$')
  fig.tight_layout()
  fig.savefig(args.path + '{}_overlap.svg'.format(args.output))

def plot_wham_results(bins, full_histogram, iter_errs, n, args, kT):
  ## plot WHAM results and individual pmfs
  plt.rc('font',size=12)
  fig, ax = plt.subplots(constrained_layout=True)
  #for hist in histograms:
  #  ax.plot(bins, hist)
  bins, full_histogram = np.squeeze(bins), np.squeeze(full_histogram)
  #errors = np.squeeze(errors)
  #ax.fill_between(bins, -kT*np.log(full_histogram)-errors, -kT*np.log(full_histogram)+errors, alpha=0.5, edgecolor='k', facecolor='0.5', zorder=49)
  #pmfs = -kT*np.log(hists) - biases
  #for pmf in pmfs.T:
  #  ax.plot(bins, pmf)
  ax.plot(bins, -kT*np.log(full_histogram), label='{}, {:.2f}'.format(n,iter_errs[-1]), zorder=50)
  ax.set_xlabel(r'$N$')
  ax.set_ylabel(r'$\Delta F$')
  #ax.invert_xaxis()
  #ax.set_ylim(bottom=-300)
  ax.legend(frameon=False, title='it, err')
  fig.savefig(args.path + '{}_wham.svg'.format(args.output))

  return bins, full_histogram

def plot_convergence(iter_errs, args, kT):
  fig, ax = plt.subplots()
  ax.semilogy(iter_errs, 'k')
  #for err in conv:
  #  ax.semilogy(err)
  ax.set_xlabel('steps')
  ax.set_ylabel('total error')
  #ax.set_ylim(bottom=0)
  fig.tight_layout()
  fig.savefig(args.path + '{}_convergence.svg'.format(args.output))

