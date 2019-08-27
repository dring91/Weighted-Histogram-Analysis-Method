import numpy as np
from sys import exit
from argparse import ArgumentParser

from input_output import *
from plot_results import *
from compute_wham import *

def main():
  """
  Program to compute full free energy curves using the Weighted Histogram Analysis Method.
  """

  parser = ArgumentParser()
  parser.add_argument('-i','--input',nargs='+',action='append',help="""Data files from umbrella sampling.""")
  parser.add_argument('-o','--output',help="""Name for output files.""")
  parser.add_argument('-p','--path',help="""Path to write output files to.""")
  parser.add_argument('-k',type=float,nargs='+',help="""Biasing strength parameter.""")
  parser.add_argument('-x',type=float,nargs='+',help="""Bias location parameter.""")
  parser.add_argument('-g','--guess',type=float,nargs='+',help="""Initial guess for wham iteration.""")
  args = parser.parse_args()

  # Specify the temperature and other important parameters
  kT = 0.7
  zmin, zmax, dz = 0, 100, 1
  max_iters = 1000

  # Check for changes in wham configuration to avoid overriding a previous wham result
  check_config_change(args)

  # define one grid so that all histograms can be easily added together
  grid = np.arange(zmin,zmax,dz)
  bins = np.diff(grid)/2 + grid[:-1]

  # reshape input data, create histograms and plot the histograms to check overlap
  args.input = np.reshape(args.input, (-1,1))
  hists, errs, ns = create_histograms(args,grid)
  plot_histogram_overlap(bins, hists, args, kT)

  # reshape arrays so they can be broadcast together properly
  bins = np.reshape(bins, (-1,1))
  args.x, args.k, ns = np.reshape(args.x, (1,-1)), np.reshape(args.k, (1,-1)), np.reshape(ns, (1,-1))

  # setup WHAM inputs
  biases = bias(bins, args.k, args.x)
  guess = np.ones(ns.shape)
  
  # compute wham
  result, full_histogram, full_errors, iter_errs, n = converge_wham(max_iters, [guess, bins, hists, errs, ns, biases, 1, kT])
  #histograms, errors, conv = bootstrap(len(result), 5, [result, bins, hists, errs, ns, biases, 1, kT]) 

  # plot wham results and convergence
  bins, full_histogram = plot_wham_results(bins, full_histogram, iter_errs, n, args, kT)
  plot_convergence(iter_errs, args, kT)

  # write out the final results
  write_results(bins, hists, ns, result, args)

if __name__ == '__main__':
  main()
