# Political Science Computational Laboratory

## Description

**pscl** is an R package providing classes, methods and test data for

- Bayesian analysis of roll call data and other binary, item-response type data (e.g., from surveys or educational testing).

- elementary Bayesian statistics

- maximum likelihood estimation of zero-inflated and hurdle models for count data

- utility functions

## Historical note

The Political Science Computational Laboratory was the name of Simon Jackman's research group at Stanford University's Department of Political Science (2002-2014), where this package was first developed.

The hurdle and count data models were extensively re-written and updated by Achim Zeileis and Christian Kleiber.  

## Installation

Most users should use latest stable release of the package, which can be installed from [CRAN](https://cran.r-project.org/) by running
```R
install.packages("pscl")
```

The development version can be installed directly from [GitHub](https://github.com/atahk/pscl) by running
```R
install.packages("devtools") ## if not already installed
library(devtools)
install_github("atahk/pscl")
```
