# Getting a server mode

## Solution 1 : Shiny
### pros:
 * Use R code and all R libs
 * simple web interface with shiny lib

### cons:
 * Require a [shiny server](https://www.rstudio.com/products/shiny/shiny-server/): the open source server as limited functions (no authentification, no ssl, no multiple process, no monitoring)

### implementation
 * upload file with [fileInput](http://shiny.rstudio.com/gallery/file-upload.html)
 * start tami in background using batch processing on cluster machine with R lib (see available scheduler at [HPC & R](https://cran.r-project.org/web/views/HighPerformanceComputing.html) :
  - [flowr](https://cran.r-project.org/web/packages/flowr/index.html)
  - [batchjobs](https://github.com/tudo-r/BatchJobs)
  - [batch](https://cran.r-project.org/web/packages/batch/index.html)
 * return data with a download link and/or table ([example](http://shiny.rstudio.com/gallery/file-download.html) for download link + table function)

## Installation
 * for [shinyserver](https://www.rstudio.com/products/shiny/download-server/) see [installation reference](http://docs.rstudio.com/shiny-server/#installation)
  - should be in a docker + local port and using nginx as proxy
  - with flowr or batchjobs require cluster access
