library(mixOmics)
library(lmms)
library(reshape2)
library(knitr)
library(parallel)
library(nlme)
library(lmeSplines)

library(dynOmics)
library(cluster)
library(mclust)
library(edci)
library(clValid)
library(devtools)
library(dynOmics)

sourceFolder = function(folder, recursive = FALSE, ...)
{
  files = list.files(folder, pattern = "[.][rR]$",
                     full.names = TRUE, recursive = recursive)
  if (!length(files))
    stop(simpleError(sprintf('No R files in folder "%s"', folder)))
  src = invisible(lapply(files, source, ...))
  message(sprintf('%s files sourced from folder "%s"', length(src), folder))
}
sourceFolder("/home/antoine/Documents/timeOmics_dev/R/Code-MAW/")
