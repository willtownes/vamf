Varying-Censoring Aware Matrix Factorization
=======

This is an R package for the VAMF algorithm implementation. It is still under development. For more information see the [BioRxiv preprint](http://www.biorxiv.org/content/early/2017/07/21/166736).

## Installation

First, make sure [rstan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) is installed and working. Then you should be able to use

```r
install.packages("devtools")
devtools::install_github("willtownes/vamf")
```

## Usage

So far there is only one exported function `vamf`. Check out the help file in an interactive session or the vignette.

## Notes

The core of VAMF is written in the probabilistic programming language [Stan](http://mc-stan.org/). See `exec/vamf.stan` in the repository.

## Contact

Developer: Will Townes
