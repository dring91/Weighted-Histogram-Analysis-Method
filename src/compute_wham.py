import numpy as np
from scipy.signal import fftconvolve
from scipy.linalg import norm
from scipy.integrate import simps

def compute_inefficiency(data,method=None):
  """
  Computes the autocorrelation decay time to measure roughly the number of independent samples
  """

  if method == None or method == "constant":
    ineff = 1
  elif method == "IACT":
    # compute integrated autocorrelation time
    autocorr = autocorrelation(data)
    stop = np.argmax(autocorr < 0.05)
    IACT = np.sum(autocorr[:stop])
    ineff = 1+2*IACT
  return ineff

def autocorrelation(data):
  """
  Computes autocorrelation function using fftconvolve.
  """

  dxs = data[:,1] - data[:,1].mean()
  normalize = np.mean(dxs**2)
  samples = np.flip(np.arange(len(data[:,1]))+1)

  autocorr = fftconvolve(dxs, np.flip(dxs), mode='full')
  autocorr = bifold(autocorr)
  autocorr = autocorr/normalize/samples

  return autocorr

def bifold(corr):
  """
  Helper function for computing IACT.
  """

  l = len(corr) // 2
  folded = (np.flip(corr[:l+1]) + corr[l:])/2
  return folded

def converge_wham(max_iters, wham_opts):
  """
  Iterates wham and stops if the iteration error blows up, or if the number of max_iters is reached.
  """

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

def bias(z, k, z0): 
  """ The functional form of the biasing potential."""
  return k/2*(z-z0)**2

def recombine(weights, hists, constant):
  """ Computes the linear combination of weights and histograms. """

  return np.sum(weights * hists, axis=1, keepdims=True)

def weight(ns, biases, free_energies, kT):
  """ Functional form of the computed weights. """
  #return ns*bay / np.sum(ns*bay*np.exp(-(biases-free_energies)/kT),axis=1,keepdims=True)
  return 1 / np.sum(ns*np.exp(-(biases-free_energies)/kT),axis=1,keepdims=True)

def wham(free_energies, grid, hists, errors, ns, biases, constant, kT):
  """
  Computes a new estimate for the weights using the WHAM equations.
  """

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

