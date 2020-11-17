
# treeImpulses

<!-- badges: start -->
[![R-CMD-check](https://github.com/r-lib/usethis/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/usethis/actions)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->


The treeImpulses package produces and plots tree-based impulse response functions. 
Users may choose from two random forest variants, implement instrumental variable, allow state-dependence, and choose lags based on information criteria or forecasting accuracy. 

---

## Installation
The package may be downlaoded or cloned from github:  

    devtools::install_github('tylerJPike/treeImpulses')

## Usage example 
This short example shows how to calculate impulse responses based on a system of three simulated variables. While the user may choose between trees based on [Athey and Wager (2019)](https://arxiv.org/abs/1510.04342) or [Oprescu et al (2019)](http://proceedings.mlr.press/v97/oprescu19a.html), the following impulse responses are calculated using the honest forest of Wager and Athey (2019).

    # import library
    library(treeImpulses)

    # generate data
    n = 300
    x = rnorm(n)
    y = rpois(n)
    z = 0.8*x + 0.1*y + 0.1*rnorm(n)
    Data = data.frame(x, y, z)

    # estimate impulse response  
    # (using a honest forest)
    impulseResponses = 
        localProjectionIRF_VAR(
            data = Data,                 
            shocks = c('x','y','z'),                
            targets = c('x','y','z'),
            horizons = 12,
            lags = 4,
            lags_max = 12,
            standardize = FALSE,       
            confidence = 0.1,          
            states = NULL,             
            NW = FALSE,                
            NW_lags = NULL,            
            NW_prewhite = NULL,        
            engine = 'AW',              
            treeBag = 5,               
            cores = -1,                
            seed = 1116,               
            IV = TRUE,                 
            Cholesky = FALSE)

    # plot the impulse responses
    plot(impulseResponses)

## Runtime notes
1. It is advised to run the routine in parallel, controled by the `cores` command
2. `NW` commands are only available with the to OLS estimation engine
