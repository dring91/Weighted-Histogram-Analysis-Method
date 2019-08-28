# Weighted Histogram Analysis Method (WHAM)

## Description

The weighted histogram analysis method is a way to stitch together umbrella sampling simulations. The probability is computed as


$$ P(\xi) = \frac{\sum_{i=1}^{N_w} g_i^{-1}h_i(\xi)}{\sum_{j=1}^{N_w}n_jg_j^{-1}exp[-\beta(w_j(\xi)-f_j)]} $$

where $N_w$ is the number of windows, $\xi$ is the reaction coordinate, $h_i$ is the histogram from the ith window, $w_i$ is the bias in the ith window, $f_j$ is the free energy offset in the ith window, $g_i$ is the statistical inefficiency in the ith window, and $n_j$ is the number of samples in the jth window. The statistical inefficiency is computed as $g_i=1+2\tau_i$ where $\tau_i$ is the integrated autocorrelation time for window i.

The wham equations are iterated by computing $P(\xi)$ using an initial guess for $\{f_i\}$ and using the equation below to compute the next set of $f_i$

$$ exp(-\beta f_j) = \int d\xi P(\xi)exp[-\beta w_j(\xi)] $$


## Run Instructions

Execute or submit the script 'submit.bash' in the top level directory.

## Dependencies

argparse, numpy, scipy
