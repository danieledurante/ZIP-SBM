Tutorial
================
This tutorial contains guidelines and code to perform the analyses for the second scenario (`simulated_network_2.RData`) considered in the simulation study of the article [**Zero-Inflated Stochastic Block Modeling of Efficiency-Security Tradeoffs in Weighted Criminal Networks**](https://arxiv.org/abs/2410.23838). In particular, you will find a detailed step-by-step guide and `R` code to **implement the Gibbs sampler presented in the article** for ZIP-SBM and to **reproduce the results for the second scenario** presented in Table 1 of the article. For implementation purposes, please **execute the code below considering the same order in which is presented**.

Import the data
================
To start the analysis, **set the working directory** in this `Tutorial` folder (where `simulated_network_2.RData` are located). Once this has been done, **clean the workspace, and load the data along with useful** `R` **packages**.

``` r
rm(list=ls())
source("zip_sbm_source.R")
Rcpp::sourceCpp('stirling.cpp')

library(reshape)
library(gdata)
library(igraph)
library(mcclust.ext)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(coda)
library(dummies)
library(randnet)
library(greed)
library(LaplacesDemon)

load("simulated_network_2.RData")
V <- dim(Y)[1]

# note that Y must have diagonal equal to 0
diag(Y)
```

As discussed in the article, the network under analysis has *V=80* nodes and *5* groups displaying heterogenous block architectures beyond simple communities, combined with a higher prevalence of covert structures [see second matrix of Figure 3 in the article]. 

Posterior computation for ZIP-SBM and its competitors (ESBM, P-SBM)
================
This section contains the code to **implement the Gibbs sampler for ZIP-SBM** [function `esbm_zip()`], and its competitors ESBM [function `esbm_binary()`] and P-SBM [function `esbm_count()`]. See the source code `zip_sbm_source.R` for a detailed description of the inputs and the outputs of these functions.

ESBM (first competitor)
------------------
Let us first perform posterior computation under the **ESBM**, which provides a competitor of ZIP-SBM. This alternative is originally developed for binary edges and, hence, is applied to a dichotomized version of the observed network. To perform posterior computation under ESBM, execute the code below.

``` r
N_iter <- 20000
burn_in <- 10000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

# define the vector with node attributes
my_x <- z_noise

# set hyperparameters for the Dirichlet-Multinomial cohesion function (see the article)
my_alpha_xi <- rep(1,5)

# run Gibbs sampler to sample the cluster indicators
Z_binary <- esbm_binary(Y_bin=(Y>0)*1, my_seed, N_iter, my_z, a=1, b=1, gamma_GN = 0.3,x = my_x, alpha_xi = my_alpha_xi)

# compute a point estimate of z
c_Z_binary <- pr_cc(Z_binary[,(burn_in+1):N_iter])
memb_Z_binary <- minVI(c_Z_binary,method="avg",max.k=20)
hat_Z_binary <- memb_Z_binary$cl

# run Gibbs sampler to sample the block probabilities given hat_z_binary
pi_binary <- sample_pi(z=hat_Z_binary,Y_bin=(Y>0)*1, my_seed, N_iter, a=1, b=1)
```

Once the above steps have been done, **save the output** in the file `binary_SBM.RData`.

``` r
# save the output
save(Z_binary,pi_binary,file="binary_SBM.RData")
rm(Z_binary,pi_binary)
```

P-SBM (second competitor)
------------------
Let us now perform posterior computation under the **P-SBM**, which provides another relevant competitor of ZIP-SBM. This alternative is developed for counts edges (but does not consider zero-inflation) and, hence, is applied directly to the observed network. To perform posterior computation under P-SBM, execute the code below.

``` r
N_iter <- 20000
burn_in <- 10000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

# define the vector with node attributes
my_x <- z_noise

# set hyperparameters for the Dirichlet-Multinomial cohesion function (see the article)
my_alpha_xi <- rep(1,5)

# run Gibbs sampler to sample the cluster indicators
Z_count <- esbm_count(Y_count=Y, my_seed, N_iter, my_z, a_1=1, a_2=1, gamma_GN = 0.3,x = my_x, alpha_xi = my_alpha_xi)

# compute a point estimate of z
c_Z_count <- pr_cc(Z_count[,(burn_in+1):N_iter])
memb_Z_count <- minVI(c_Z_count,method="avg",max.k=20)
hat_Z_count <- memb_Z_count$cl

# run Gibbs sampler to sample the block lambdas given hat_z_lambda
lambda_count <- sample_lambda(z=hat_Z_count,Y_count=Y, my_seed, N_iter, a_1=1, a_2=1)
```

Once the above steps have been done, **save the output** in the file `count_SBM.RData`.

``` r
# save the output
save(Z_count,lambda_count,file="count_SBM.RData")
rm(Z_count,lambda_count)
```

ZIP-SBM (proposed model)
------------------
Let us finally perform posterior computation under the proposed **ZIP-SBM** applied directly to the observed network. To perform posterior computation under ZIP-SBM, execute the code below.

``` r
N_iter <- 20000
burn_in <- 10000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

# define the vector with node attributes
my_x <- z_noise

# set hyperparameters for the Dirichlet-Multinomial cohesion function (see the article)
my_alpha_xi <- rep(1,5)

# run Gibbs sampler to sample the cluster indicators
start_time <- Sys.time()

Z_zip <- esbm_zip(Y=Y, my_seed, N_iter, my_z, a=1, b=9, a_1=1, a_2=1, gamma_GN = 0.3,x = my_x, alpha_xi = my_alpha_xi)

end_time <- Sys.time()
(end_time - start_time)/20
#7.4 seconds

# compute a point estimate of z
c_Z_zip <- pr_cc(Z_zip[,(burn_in+1):N_iter])
memb_Z_zip <- minVI(c_Z_zip,method="avg",max.k=20)
hat_Z_zip <- memb_Z_zip$cl

# run Gibbs sampler to sample the block lambdas given hat_z_lambda
pi_lambda_zip <- sample_lambda_pi(z=hat_Z_zip,Y=Y, my_seed, N_iter, a=1, b=9, a_1=1, a_2=1)
```

Once the above steps have been done, **save the output** in the file `zip_SBM.RData`.

``` r
save(Z_zip,pi_lambda_zip,file="zip_SBM.RData")
rm(Z_zip,pi_lambda_zip)
```

Posterior inference under ZIP-SBM and its competitors (ESBM and P-SBM) [Table 1: Scenario 2]
================
This section contains the **code to perform estimation, uncertainty quantification and model selection for ZIP-SBM and its competitors (ESBM and P-SBM)** leveraging the samples from the previous Gibbs samplers. In particular, we **reproduce the analyses in Table 1 of the article**, for **scenario 2**. To accomplish this goal let us first **clear the workspace** and **load useful libraries and data**.

``` r
rm(list=ls())
source("zip_sbm_source.R")
Rcpp::sourceCpp('stirling.cpp')

library(reshape)
library(gdata)
library(igraph)
library(mcclust.ext)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(coda)
library(dummies)
library(randnet)
library(greed)
library(aricode)

load("simulated_network_2.RData")
V <- dim(Y)[1]

# note that Y must have diagonal equal to 0
diag(Y)```

Once this has been done, the performance measures discussed in the article for ESBM, P-SBM and the proposed ZIP-SBM (under scenario 2) can be obtained by executing the code below.

``` r
############################################################
# ------------------------------------
# ESBM
# ------------------------------------
############################################################
load("binary_SBM.RData")

N_iter <- 20000
burn_in <- 10000
V <- dim(Y)[1]

# posterior mean of VI distance from z_0
VI(z_0,t(Z_binary[,(burn_in+1):N_iter]))

# compute a point estimate of z
c_Z_binary <- pr_cc(Z_binary[,(burn_in+1):N_iter])
memb_Z_binary <- minVI(c_Z_binary,method="avg",max.k=20)
hat_Z_binary <- memb_Z_binary$cl

# estimated H
max(hat_Z_binary )

# horizontal bound of the credible ball
credibleball(memb_Z_binary$cl,t(Z_binary[,(burn_in+1):N_iter]))[[5]]

# VI distance of point estimate from z_0
VI(z_0,t(hat_Z_binary))

# Normalized mutual information
NMI(z_0,hat_Z_binary)

  
############################################################
# ------------------------------------
# P-SBM
# ------------------------------------
############################################################
load("count_SBM.RData")

N_iter <- 20000
burn_in <- 10000
V <- dim(Y)[1]

# posterior mean of VI distance from z_0
VI(z_0,t(Z_count[,(burn_in+1):N_iter]))

# compute a point estimate of z
c_Z_count <- pr_cc(Z_count[,(burn_in+1):N_iter])
memb_Z_count <- minVI(c_Z_count,method="avg",max.k=20)
hat_Z_count <- memb_Z_count$cl

# estimated H
max(hat_Z_count)

# horizontal bound of the credible ball
credibleball(memb_Z_count$cl,t(Z_count[,(burn_in+1):N_iter]))[[5]]
#[1] 0.6956093

# VI distance of point estimate from z_0
VI(z_0,t(hat_Z_count))

# Normalized mutual information
NMI(z_0,hat_Z_count)


############################################################
# ------------------------------------
# ZIP-SBM
# ------------------------------------
############################################################
load("zip_SBM.RData")

N_iter <- 20000
burn_in <- 10000
V <- dim(Y)[1]

# posterior mean of VI distance from z_0
VI(z_0,t(Z_zip[,(burn_in+1):N_iter]))

# compute a point estimate of z
c_Z_zip <- pr_cc(Z_zip[,(burn_in+1):N_iter])
memb_Z_zip <- minVI(c_Z_zip,method="avg",max.k=20)
hat_Z_zip <- memb_Z_zip$cl

# estimated H
max(hat_Z_zip)

# horizontal bound of the credible ball
credibleball(memb_Z_zip$cl,t(Z_zip[,(burn_in+1):N_iter]))[[5]]

# VI distance of point estimate from z_0
VI(z_0,t(hat_Z_zip))

# Normalized mutual information
NMI(z_0,hat_Z_zip)
 ```
