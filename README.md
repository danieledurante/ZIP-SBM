# ZIP-SBM: Zero-inflated Poisson stochastic block model

This repository is associated with the article **Zero-Inflated Stochastic Block Modeling of Efficiency-Security Tradeoffs in Weighted Criminal Networks**, and aims at providing detailed tutorials and codes to implement the **ZIP-SBM** presented in the article*.

The documentation is organized as described below.  

- `zip_sbm_source.R`.  It contains all the **source** `R` **functions** which are required to perform posterior computation and inference under the proposed **ZIP-SBM**  and its relevant competitors (i.e., ESBM and P-SBM).

- `Tutorial`. It contains data [see `simulated_network_2.RData`], source codes [see `zip_sbm_source.R` and  `stirling.cpp`], and a step-by-step tutorial [see `scenario_2.md`]  to implement the proposed **ZIP-SBM**  and its relevant competitors (ESBM and P-SBM), with a focus on reproducing the results for **scenario 2** in Table 1 in the article.

The analyses are performed with a  **MacBook Air (M1, 2020), CPU 8–core and 8GB RAM (macOS Monterey, version 12.5)**, using the `R` version 4.2.2

All the above functions rely on a basic and reproducible `R` implementation, mostly meant to provide a clear understanding of the computational routines associated with the proposed model. Optimized computational routines relying on C++ coding can be easily considered. 

A seed is set at the beginning of each code to ensure full reproducibility. Slight changes in the final numbers (if any) depend on which version of the external `R` libraries employed has been used in the implementation of the codes. This is due to internal changes of certain functions when versions of some packages have been updated. However, the magnitude of these minor variations (if any) is negligible and does not affect the final conclusions.