# Weighted Histogram Analysis Method (WHAM)

## Description

The weighted histogram analysis method is a way to stitch together umbrella sampling simulations. The probability is computed as


<p align="center"><img src="/tex/cba3c8088b7b1436b4a2f9bc68d40baa.svg?invert_in_darkmode&sanitize=true" align=middle width=288.6767961pt height=48.9591135pt/></p>

where <img src="/tex/4070078cd9502b3b3d8207574f65c7b2.svg?invert_in_darkmode&sanitize=true" align=middle width=23.027089799999988pt height=22.465723500000017pt/> is the number of windows, <img src="/tex/85e60dfc14844168fd12baa5bfd2517d.svg?invert_in_darkmode&sanitize=true" align=middle width=7.94809454999999pt height=22.831056599999986pt/> is the reaction coordinate, <img src="/tex/ddd3bc35b936d6a00e6a81cab0061f32.svg?invert_in_darkmode&sanitize=true" align=middle width=14.12201339999999pt height=22.831056599999986pt/> is the histogram from the ith window, <img src="/tex/c2a29561d89e139b3c7bffe51570c3ce.svg?invert_in_darkmode&sanitize=true" align=middle width=16.41940739999999pt height=14.15524440000002pt/> is the bias in the ith window, <img src="/tex/ac9424c220341fa74016e5769014f456.svg?invert_in_darkmode&sanitize=true" align=middle width=14.152495499999992pt height=22.831056599999986pt/> is the free energy offset in the ith window, <img src="/tex/681a37b53b66acbc455e39ca3e6f1c41.svg?invert_in_darkmode&sanitize=true" align=middle width=12.49148174999999pt height=14.15524440000002pt/> is the statistical inefficiency in the ith window, and <img src="/tex/54158e2c605c3ecf783cdc13e7235676.svg?invert_in_darkmode&sanitize=true" align=middle width=15.971386199999989pt height=14.15524440000002pt/> is the number of samples in the jth window. The statistical inefficiency is computed as <img src="/tex/548b958ca26501ea2d856e5ac7a59f04.svg?invert_in_darkmode&sanitize=true" align=middle width=83.59762454999999pt height=21.18721440000001pt/> where <img src="/tex/e7cdf5013524d24e01bb7ecb5878d45f.svg?invert_in_darkmode&sanitize=true" align=middle width=11.83700594999999pt height=14.15524440000002pt/> is the integrated autocorrelation time for window i.

The wham equations are iterated by computing <img src="/tex/37cc584b5a616c5d79ce00414531e92a.svg?invert_in_darkmode&sanitize=true" align=middle width=33.57029279999999pt height=24.65753399999998pt/> using an initial guess for <img src="/tex/0573c107229578771a5fb170a8d63bb6.svg?invert_in_darkmode&sanitize=true" align=middle width=29.959203449999993pt height=24.65753399999998pt/> and using the equation below to compute the next set of <img src="/tex/9b6dbadab1b122f6d297345e9d3b8dd7.svg?invert_in_darkmode&sanitize=true" align=middle width=12.69888674999999pt height=22.831056599999986pt/>

<p align="center"><img src="/tex/7f148653cc1925759e59d02105dc0280.svg?invert_in_darkmode&sanitize=true" align=middle width=264.03203805pt height=36.53007435pt/></p>


## Run Instructions

Execute or submit the script 'submit.bash' in the top level directory.

## Dependencies

argparse, numpy, scipy
