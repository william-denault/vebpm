## Overview

Perform empirical Bayes inference on $\mu_i$ in the model 
$$x_i|\mu_i\sim \text{Poisson}(s_i\exp(\mu_i)),$$
$$\mu_i\sim g(\cdot).$$

## Installation

```r
# install.packages('devtools')

## Make sure R package ebnm is installed
# devtools::install_github('stephenslab/ebnm')

devtools::install_github('DongyueXie/vebpm')
```
