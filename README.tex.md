# Weighted Histogram Analysis Method (WHAM)

## Description

The weighted histogram analysis method is a way to stitch together umbrella sampling simulations.


$$ P(\xi) = \frac{\sum{i=1}^{N_w} g_i^{-1}h_i(\xi)}{\sum_{j=1}^{N_w}n_jg_j^{-1}exp[-\beta(w_j(\xi)-f_j)]} $$


$$ exp(-\beta f_i) = \int d\xi exp[-\beta w_j(\xi)]P(\xi) $$


## Run Instructions

## Dependencies
