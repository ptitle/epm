  <!-- badges: start -->
  [![R-CMD-check](https://github.com/ptitle/epm/workflows/R-CMD-check/badge.svg)](https://github.com/ptitle/epm/actions)
  [![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/epm)](https://github.com/r-hub/cranlogs.app)
  [![cran version](https://www.r-pkg.org/badges/version/epm)](https://cran.r-project.org/package=epm)
  [![DOI](https://zenodo.org/badge/348488706.svg)](https://zenodo.org/badge/latestdoi/348488706)
  <!-- badges: end -->
# epm
## EcoPhyloMapper R package

This R package facilitates the aggregation of species' geographic ranges from vector or raster spatial data, and enables the calculation of various morphological and phylogenetic community metrics across geography, with a particular focus on handling of geometric morphometric data.

The following R package dependencies are required:

`terra, sf, ape, viridisLite, pbapply, Rcpp`

Additionally, the R packages `exactextractr`, `tmap`, `data.table` and 
`spdep` provide 
additional functionality or speed improvements, and are recommended.

The current release of epm can be [downloaded from CRAN]( https://cran.r-project.org/package=epm).

To install the latest development version from GitHub, run the following in R: `remotes::install_github('ptitle/epm')`

[Check out the wiki for explanations and discussions of the available options.](https://github.com/ptitle/epm/wiki)

## Citation
If you use this R package, please cite the accompanying article: <br>
Pascal O. Title, Donald L. Swiderski and Miriam L. Zelditch. 2022. EcoPhyloMapper: an R package for integrating geographic ranges, phylogeny, and morphology. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.13914


