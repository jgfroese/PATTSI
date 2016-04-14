# PATTSI
Spatial PATTern SuItability analysis

© 2016, Jens G. Froese (jens.froese@uq.net.au)

## Description
This R scipt describes a methodology for measuring and relating spatial patterns of habitat variables to their suitability for a given species response (e.g. presence, abundance, habitat suitability etc.)

---
PATTSI requires an **input raster layer (.TIFF)** describing the spatial distribution of a given habitat variable.

---
Spatial patterns of habitat variables are measured using moving window analysis, whereby the value of a focal raster pixel is computed by summarizing the values of all neighbouring raster pixels contained within an analysis window. This window is incrementally moved across the study area, centering on each raster pixel.

Currently, three measurements of spatial patterns (landscape metrics) are currently implemented:
* maximum distance-weighted value of a numerical habitat variable (using function `focal {raster}`)
* average value of a numerical habitat variable (using the `ESRI ArcGIS Focal Statistics` tool)
* average distance-weighted value of a numerical habitat variable (using the `ESRI ArcGIS Focal Statistics` tool)

External `ESRI ArcGIS` software was used to circumvent technical issues with treatment of `NA` / `0` values in `focal {raster}, fun = mean`, which should be resolved in future. PATTSI could also be adapted to other landscape metrics of numerical as well as categorical habitat variables.

---
PATTSI also requires an **input CSV file** relating a given landscape metric to a numerical suitability index, i.e. dscribing the species response to spatial patterns of a habitat variable. 

Currently, these responses have been elicited from experts with knowledge of the study system. However, they could equally be parameterized from empirical data.

---
PATTSI returns an **output raster layer (.TIFF)** describing the spatial distribution of a suitability index for a given habitat variable.

---
A worked example is also included: `PATTSI-example.R`.

## References
* Bates, D. and Maechler, M. 2015. Package 'Matrix': sparse and dense matrix classes and methods. URL http://Matrix.R-forge.R-project.org/.
* Hijmans, R.J. 2015. Package 'raster': geographic data analysis and modeling. URL http://cran.r-project.org/web/packages/raster/.
* RCoreTeam 2015. R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
* Scroggie, M. 2012. Applying a circular moving window filter to raster data in R. URL https://scrogster.wordpress.com/2012/10/05/applying-a-circular-moving-window-filter-to-raster-data-in-r/.
