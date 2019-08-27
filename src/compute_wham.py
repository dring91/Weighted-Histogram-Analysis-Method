import numpy as np
from scipy.signal import fftconvolve
from scipy.linalg import norm
from scipy.integrate import simps

def compute_inefficiency(data,method=None):
  if method == None or method == "constant":
    ineff = 1
  elif method == "IACT":
    # compute integrated autocorrelation time
    autocorr = autocorrelation(data)
    stop = np.argmax(autocorr < 0.05)
    IACT = np.sum(autocorr[:stop])
    ineff = 1+2*IACT
  return ineff

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

def weight(ns, biases, free_energies, kT):
  #return ns*bay / np.sum(ns*bay*np.exp(-(biases-free_energies)/kT),axis=1,keepdims=True)
  return 1 / np.sum(ns*np.exp(-(biases-free_energies)/kT),axis=1,keepdims=True)
