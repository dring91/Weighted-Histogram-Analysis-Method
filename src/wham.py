import numpy as np
from argparse import ArgumentParser
from scipy.optimize import curve_fit
from scipy.linalg import norm
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.integrate import simps

def bias(z, k, z0): return k/2*(z-z0)**2

def recombine(weights, hists, constant):
  return np.sum(weights * hists, axis=1, keepdims=True)

def weight(ns, biases, free_energies, kT):
  return ns / np.sum(ns*np.exp(-(biases-free_energies)/kT),axis=1,keepdims=True)

def wham(free_energies, grid, hists, errors, ns, biases, constant, kT):

  ## recombine histograms and unbias
  weights = weight(ns, biases, free_energies, kT)
  full_histogram = recombine(weights, hists, constant)
  integrand = full_histogram * np.exp(-biases/kT)

  ## integrate to get free energies
  integral = simps(integrand, grid, axis=0)
  free_energies = -kT * np.log(constant * integral)

  ## recalculate recombined histogram with new weights
  weights = weight(ns, biases, free_energies, kT)
  full_histogram = recombine(weights, hists, constant)
  full_errors = kT * np.sqrt(np.sum(errors**2*weights**2,axis=1,keepdims=True)) / np.abs(full_histogram)

  return free_energies, full_histogram, full_errors

def main():

  parser = ArgumentParser()
  parser.add_argument('-i','--input',nargs='+',action='append')
  parser.add_argument('-o','--output')
  parser.add_argument('-p','--path')
  parser.add_argument('-k',type=float,nargs='+')
  parser.add_argument('-x',type=float,nargs='+')
  parser.add_argument('-g','--guess',type=float,nargs='+',help='initial guess')
  args = parser.parse_args()

  ## check for meta file (describes which data sets are being whammed)
  try:
    with open(args.path + 'args','r') as file:
      for line in file:
        if line.strip() not in args.input[0]:
          print("""Input files differ from wham configuration. 
                   Change folders or delete args metafile before trying to wham again.\n
                """)
          exit()
  ## in the absence of meta file, create one and print out parser args
  except FileNotFoundError:
    with open(args.path + 'args','w') as file:
      for filename in args.input[0]:
        file.write('{}\n'.format(filename))

  ## define grid so that all histograms can be easily added together
  zmin, zmax, dz = 0, 360, 1
  grid = np.arange(zmin,zmax,dz)
  bins = np.diff(grid)/2 + grid[:-1]

  args.input = np.reshape(args.input, (-1,1))

  ## check that there are enough shifts
  if len(args.guess) == None: 
    args.guess = np.full(len(args.x), 0)
  elif len(args.guess) < len(args.x): 
    args.guess = np.pad(args.guess, (0,len(args.x)-len(args.guess)), 'constant', constant_values=(0,0))

  ns, hists, errs = [], [], []
  for reps in args.input:
    n, hist = 0, []
    for filename in reps:
      ## read in converged data (vary amount based on args.x)
      data = np.genfromtxt(filename, skip_header=30000, usemask=True)
      data = data[::50]

      ## histogram N for wham
      rep_hist, _ = np.histogram(data[:,1],density=True,bins=grid)
      n += len(data)

      ## average reps together for better statistics without more variables for WHAM
      hist.append(rep_hist)

    ## save the histograms and the number of data points for WHAM reconstruction
    hists.append(np.mean(hist,axis=0))
    errs.append(np.std(hist,axis=0))
    ns.append(n)

  hists = np.column_stack(hists)
  errs = np.column_stack(errs)

  ## check for overlap so that iteration doesn't diverge
  fig, ax = plt.subplots()
  lines = ax.plot(bins, hists)
  ax.set_xlabel(r'$N$')
  ax.set_ylabel(r'$P(N)$')
  fig.tight_layout()
  fig.savefig(args.path + 'overlap_test.png')

  ## reshape arrays so they can be broadcast together properly
  bins = np.reshape(bins, (-1,1))
  args.x, args.k, ns = np.reshape(args.x, (1,-1)), np.reshape(args.k, (1,-1)), np.reshape(ns, (1,-1))

  ## apply WHAM and iterate
  biases = bias(bins, args.k, args.x)
  guess = np.ones(ns.shape)
  #guess = np.array(args.guess)
  plt.rc('font',size=12)
  fig, ax = plt.subplots()

  iter_errs = []
  n_iters = 2000
  for n in range(n_iters):
    result, full_histogram, full_errors = wham(guess, bins, hists, errs, ns, biases, 1, 0.7)
    try:
      iter_err = norm(guess - result)
    except ValueError:
      print("Iteration Error infinite")
      iter_err = np.inf
      continue
    iter_errs.append(iter_err)
    guess = np.copy(result)
  def line(x,m,b): return m*x+b
  def N2h(x): return line(x,1/11.5,-13.7/11.5)
  lines, *_ = ax.errorbar(N2h(bins), -0.7*np.log(full_histogram), yerr=full_errors, color='k', label='{}, {:.2f}'.format(n,iter_err))
  pmfs = -0.7*np.log(hists) - biases
  #for pmf in pmfs.T:
  #  ax.plot(N2h(bins), pmf)
  xs = np.linspace(0,30,10)
  ax.plot(xs, line(xs,0.42,-7), color='k', linestyle='dashed')
  h = N2h(bins)
  mask = (h>-0.5) * (h<=2)
  ax.plot(h[mask], 0.69*np.pi*((h[mask]-3)**2-2.5**2)+10, color='k', linestyle='dashed')

  ax.set_xlabel(r'$N$')
  ax.set_xlabel(r'$h/\sigma$')
  ax.set_ylabel(r'$\Delta F$')
  #ax.invert_xaxis()
  #ax.set_ylim(bottom=-300)
  ax.legend(frameon=False, title='it, err')
  fig.tight_layout()
  fig.savefig(args.path + 'wham_results.png')

  fig, ax = plt.subplots()
  ax.semilogy(iter_errs)
  ax.set_xlabel('steps')
  ax.set_ylabel('total error')
  #ax.set_ylim(bottom=0)
  fig.tight_layout()
  fig.savefig(args.path + 'wham_convergence.png')

  with open(args.path + args.output, 'w') as file:
    np.savetxt(file, np.column_stack((bins, full_histogram, full_errors)), header='bins hist errs')

if __name__ == '__main__':
  main()
