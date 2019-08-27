import numpy as np
from sys import exit
from argparse import ArgumentParser
from scipy.optimize import curve_fit
from scipy.linalg import norm
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sklearn.utils import resample
from scipy.signal import fftconvolve

def bayesian_weights(n_samples):

  weights = np.sort(np.random.rand(n_samples - 1))
  weights = np.diff(np.concatenate(([0],weights,[1])))

  return weights

def bootstrap(n_windows, n_bootstraps, wham_opts):

  bay_weights = [bayesian_weights(n_windows) for _ in range(n_bootstraps)]

  max_iters = 2000
  histograms, iter_errs = [], []
  for w in bay_weights:
    wham_opts[4] *= w
    iter_err = []
    for n in range(max_iters):
      result, full_histogram, full_errors = wham(*wham_opts)
      #iter_err.append(norm(wham_opts[0] - result))
      wham_opts[0] = np.copy(result)
    wham_opts[4] /= w
    histograms.append(np.squeeze(full_histogram))
    iter_errs.append(iter_err)

  histograms = -wham_opts[-1] * np.log(histograms)
  shift = np.array([hist[np.isfinite(hist)][0] for hist in histograms])
  histograms -= shift[:,np.newaxis]
  average = np.mean(histograms, axis=0)
  errors = np.sqrt(np.sum((histograms - average)**2, axis=0)/(n_bootstraps - 1))

  return histograms, errors, iter_errs

def bifold(corr):
  l = len(corr) // 2
  folded = (np.flip(corr[:l+1]) + corr[l:])/2
  return folded

def autocorrelation(data):

  dxs = data[:,1] - data[:,1].mean()
  normalize = np.mean(dxs**2)
  samples = np.flip(np.arange(len(data[:,1]))+1)

  autocorr = fftconvolve(dxs, np.flip(dxs), mode='full')
  autocorr = bifold(autocorr)
  autocorr = autocorr/normalize/samples

  return autocorr

def bias(z, k, z0): return k/2*(z-z0)**2

def recombine(weights, hists, constant):
  return np.sum(weights * hists, axis=1, keepdims=True)

def weight(ns, biases, free_energies, kT, bay=1):
  #return ns*bay / np.sum(ns*bay*np.exp(-(biases-free_energies)/kT),axis=1,keepdims=True)
  return bay / np.sum(ns*bay*np.exp(-(biases-free_energies)/kT),axis=1,keepdims=True)

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
  zmin, zmax, dz = 0, 900, 1
  grid = np.arange(zmin,zmax,dz)
  bins = np.diff(grid)/2 + grid[:-1]

  args.input = np.reshape(args.input, (-1,1))
  print(args.input)

  ns, hists, ineffs, errs = [], [], [], []
  for reps in args.input:
    n, hist = 0, []
    for filename in reps:
      ## read in data
      data = np.genfromtxt(filename, usemask=True, skip_header=2000)

      ## histogram N for wham
      rep_hist, _ = np.histogram(data[::10,1],density=False,bins=grid)
      n += len(data[::10])

      ## average reps together for better statistics without more variables for WHAM
      hist.append(rep_hist)

      ## compute integrated autocorrelation time
      #autocorr = autocorrelation(data)
      #stop = np.argmax(autocorr < 0.05)
      #IACT = np.sum(autocorr[:stop])
      #ineff = 1+2*IACT
      ineff = 1

    ## save the histograms and the number of data points for WHAM reconstruction
    hists.append(np.mean(hist,axis=0))
    errs.append(np.std(hist,axis=0))
    ns.append(n/ineff)

  hists = np.column_stack(hists)
  errs = np.column_stack(errs)

  ## check for overlap so that iteration doesn't diverge
  fig, ax = plt.subplots()
  lines = ax.plot(bins, hists)
  ax.set_xlabel(r'$N$')
  ax.set_ylabel(r'$P(N)$')
  fig.tight_layout()
  fig.savefig(args.path + '{}_overlap.svg'.format(args.output))

  ## reshape arrays so they can be broadcast together properly
  bins = np.reshape(bins, (-1,1))
  args.x, args.k, ns = np.reshape(args.x, (1,-1)), np.reshape(args.k, (1,-1)), np.reshape(ns, (1,-1))

  ## apply WHAM and iterate
  biases = bias(bins, args.k, args.x)
  guess = np.ones(ns.shape)
  #guess = np.array(args.guess)
  plt.rc('font',size=12)
  fig, ax = plt.subplots(constrained_layout=True)

  ## break out wham loop into a function
  def converge_wham(max_iters, wham_opts):
    iter_errs = []
    for n in range(max_iters):
      result, full_histogram, full_errors = wham(*wham_opts)
      try:
        iter_err = norm(wham_opts[0] - result)
      except ValueError:
        print("Iteration Error infinite")
        exit()
      iter_errs.append(iter_err)
      wham_opts[0] = np.copy(result)
  
    return result, full_histogram, full_errors, iter_errs, n

  max_iters = 10000
  kT = 0.7
  result, full_histogram, full_errors, iter_errs, n = converge_wham(max_iters, [guess, bins, hists, errs, ns, biases, 1, kT])
  histograms, errors, conv = bootstrap(len(result), 5, [result, bins, hists, errs, ns, biases, 1, kT]) 

  ## plot WHAM results and individual pmfs
  #for hist in histograms:
  #  ax.plot(bins, hist)
  bins, full_histogram = np.squeeze(bins), np.squeeze(full_histogram)
  errors = np.squeeze(errors)
  ax.fill_between(bins, -kT*np.log(full_histogram)-errors, -kT*np.log(full_histogram)+errors, alpha=0.5, edgecolor='k', facecolor='0.5', zorder=49)
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

  fig, ax = plt.subplots()
  ax.semilogy(iter_errs, 'k')
  #for err in conv:
  #  ax.semilogy(err)
  ax.set_xlabel('steps')
  ax.set_ylabel('total error')
  #ax.set_ylim(bottom=0)
  fig.tight_layout()
  fig.savefig(args.path + '{}_convergence.svg'.format(args.output))

  ## concatenate bins and histograms and construct header
  out = np.concatenate((bins[:,np.newaxis],hists),axis=1)
  header = "".join(["{:5.0f}" for _ in range(hists.shape[1])])
  header = "bins " + header.format(*(np.squeeze(args.x).tolist()))
  with open(args.path + args.output + '_pmfs', 'w') as file:
    np.savetxt(file, out, header=header)
  with open(args.path + args.output + '_counts', 'w') as file:
    print(args.x.shape, ns.shape, result.shape)
    np.savetxt(file, np.column_stack((np.squeeze(args.x),np.squeeze(ns),result)), header='windows ineff fs')
    ## compute the number of points in each bin output as a 2d array

  #with open(args.path + args.output, 'w') as file:
  #  np.savetxt(file, np.column_stack((bins, full_histogram, errors)), header='bins hist errs')

if __name__ == '__main__':
  main()
