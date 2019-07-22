# Integration of longitudinal microbiome studies 

We propose a generic data-driven framework to integrate different types of longitudinal data measured on the same biological specimens with microbial communities data, and select key temporal features with strong associations within the same sample group. The framework ranges from filtering and modelling, to integration using smoothing splines and multivariate dimension reduction methods to address some of the analytical challenges of microbiome-derived data. We illustrate our framework on different types of multi-omics case studies in bioreactor experiments as well as human studies.

More details about the integrative framework can be found in https://www.biorxiv.org/content/10.1101/585802v1

To reproduce the examples, the user can compile the `.Rmd` files in the `Example/` folder. 
The list of required packages and their versions is available in the file "needed_packages.txt".  
Please install them before compiling the source files.

This repository is structured as follows:

## Data

Raw data used to conduct examplar analysis.
Files are loaded when the `.Rmd` files are compiled.

## Examples

Folder that contains the Rmarkdown scripts to re-run our 2 examples.

* **Baby gut microbiota developpement**: studies microbiome dynamic during the first 100 days of infant early life.
* **Waste degradation**: integrates microbiome, metabolomic and a continious response of interest.

This folder also contains the pdf output files.


## Rscripts

Contains the R files used for time-course microbiome clustering and integration.



