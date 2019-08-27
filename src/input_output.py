import numpy as np
from compute_wham import compute_inefficiency

def write_results(bins, hists, ns, result, args):
  """
  Writes all histograms to an output text file.
  """

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

def create_histograms(args, grid):
  """
  Reads in INDUS results and constructs the appropriate histograms.
  """

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

      ineff = compute_inefficiency(data)

    ## save the histograms and the number of data points for WHAM reconstruction
    hists.append(np.mean(hist,axis=0))
    errs.append(np.std(hist,axis=0))
    ns.append(n/ineff)

  hists = np.column_stack(hists)
  errs = np.column_stack(errs)

  return hists, errs, ns

def check_config_change(args):
  """
  Writes the commandline arguments to a 'config' file if 'config' isn't found. Alternately, it checks that the commandline arguments match 'config' file and halts the program if they don't. This prevents automatically overwriting whammed configurations if the input changes. It does not protect against overwriting wham data if the progran is rerun with new wham specifications.
  """

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

