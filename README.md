# __prosail__ <img src="man/figures/logo.png" align="right" alt="" width="200" />

# An R package for the simulation of canopy reflectance using the model PROSAIL (PROSPECT+SAIL).

[![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)](https://www.r-project.org/Licenses/GPL-3)
[![Build Status](https://gitlab.com/jbferet/prosail/badges/master/pipeline.svg)](https://gitlab.com/jbferet/prosail/pipelines/latest)

# 1 Install

After installing package `devtools`, the package `prosail` can be installed with the following command line in R session:
```
devtools::install_gitlab('jbferet/prosail')
```

... if you are already on this webpage, but `prosail` is still not publicly available... Lucky you!!!
then install the `getPass` package, and run this command line:

```
devtools::install_git('https://gitlab.com/jbferet/prosail',credentials = git2r::cred_user_pass('Your_Gitlab_UserName',getPass::getPass())) 
```


# 2 Tutorial

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- ```{r include = FALSE} -->
<!-- knitr::opts_chunk$set( -->
<!--   collapse = TRUE, -->
<!--   comment = "#>", -->
<!--   fig.path = "man/figures/README-", -->
<!--   out.width = "100%" -->
<!-- ) -->
<!-- ``` -->

A tutorial vignette is available [here](https://jbferet.gitlab.io/prosail/articles/prosail.html).

# 3 Citation

If you use **prosail**, please cite the following references:

## PROSPECT
Féret J-B, Gitelson AA, Noble SD & Jacquemoud S, 2017. PROSPECT-D: Towards modeling leaf optical properties through a complete lifecycle. Remote Sensing of Environment, 193, 204–215. https://doi.org/10.1016/j.rse.2017.03.004
Féret J-B, Berger K, de Boissieu F & Malenovský Z, 2020. Estimation of leaf protein and carbon-based constituent content from optical properties with the PROSPECT-PRO model. ArXiv200311961 Q-Bio.

## 4SAIL & 4SAIL2
Verhoef W & Bach H, 2007. Coupled soil–leaf-canopy and atmosphere radiative transfer modeling to simulate hyperspectral multi-angular surface reflectance and TOA radiance data. Remote Sensing of Environment, 109:166-182. doi:10.1016/j.rse.2006.12.013
Verhoef W, Jia L, Xiao Q & Su Z, 2007. Unified optical-thermal four-stream radiative transfer theory for homogeneous vegetation canopies. IEEE Transactions in Geosciences and Remote Sensing, 45:1808–1822. https://doi.org/10.1109/TGRS.2007.895844

## PROSAIL
Jacquemoud S, Verhoef W, Baret F, Bacour C, Zarco-Tejada PJ, Asner GP, François C & Ustin SL, 2009. PROSPECT+ SAIL models: A review of use for vegetation characterization. Remote Sensing of Environment, 113:S56–S66. https://doi.org/doi:10.1016/j.rse.2008.01.026
Berger K, Atzberger C, Danner M, D’Urso G, Mauser W, Vuolo F & Hank T 2018. Evaluation of the PROSAIL Model Capabilities for Future Hyperspectral Model Environments: A Review Study. Remote Sensing, 10:85. https://doi.org/10.3390/rs10010085
