# mi-by-decoding

R scripts for estimating mutual information by decoding stochastic time series data.

## Installation

First ensure that [R](http://r-project.org) is installed on your machine. Then
from this directory either:

1. install from the pre-compiled package by opening an R session and running 
   `install.packages('release/MIdecoding_1.0-1.tar.gz')`, or
2. install from source by opening a terminal window, and then running the
   commands:
   ```bash
   > R CMD build MIdecoding
   > R CMD install MIdecoding_1.0-1.tar.gz
   ```

Then in order to use the library in your R session, you need to load it by
running `library(MIdecoding)`.

## Getting started

Check the examples in the documentation for the `MIdecoding` function. These
can be found either by running `help(MIdecoding)` from your R session or in
the PDF version of the documentation (`release/MIdecoding.pdf`). You can
easily run these examples in your R session using `example(MIdecoding)`.

The `summary` function for `MIdecoding` objects also produces useful output.
See the documentation for `summary.MIdecoding` for examples.

The licence for this software can be found in the file `MIdecoding/LICENCE`.
