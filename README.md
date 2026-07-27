
# GCCfactor

<!-- badges: start -->
<!-- badges: end -->

The goal of GCCfactor is to implement estimation, model selection, and inference, as developed
in "Generalised canonical correlation estimation of the multilevel factor model".

If you encounter a bug, please file an issue on the [GitHub Issues](https://github.com/rl1081/GCCfactor/issues) page.

## Installation

You can install the development version of GCCfactor from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("rl1081/GCCfactor")
```

or install from CRAN:

``` r
install.packages("GCCfactor")
```

## Example

``` r
library(GCCfactor)

panel <- UKhouse 
est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
                        overestimate = TRUE, depvar_header = "dlPrice", 
                        i_header = "Region", j_header = "LPA_Type", t_header = "Date")
```

