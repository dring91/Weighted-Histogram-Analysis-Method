import numpy as np
from sklearn.utils import resample

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

def weight(ns, biases, free_energies, kT, bay=1):
  #return ns*bay / np.sum(ns*bay*np.exp(-(biases-free_energies)/kT),axis=1,keepdims=True)
  return bay / np.sum(ns*bay*np.exp(-(biases-free_energies)/kT),axis=1,keepdims=True)


