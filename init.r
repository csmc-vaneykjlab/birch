## Install R packages
#packages <- c("shiny", "shinyjs", "shinyalert", "shinyvalidate", "shinythemes", "shinycssloaders", "rsconnect", "DT", "tidyr", "textreadr")
#packages <- c("optparse", "dplyr", "tibble", "ggplot2", "assertthat", "readr", "purrr", "tidyverse")
packages <- c("ggpubr", "RColorBrewer", "gtable", "grid", "gridExtra", "tools", "missRanger", "plotly")
install_if_missing <- function(p) {
    print(p)
    if (!p %in% rownames(installed.packages())) {
        install.packages(p, clean=TRUE, quiet=TRUE)
    }
}
invisible(sapply(packages, install_if_missing))

packages_bio <- c("pcaMethods", "SummarizedExperiment", "MSnbase")
install_if_missing2 <- function(h) {
#    if (!require("BiocManager", quietly = TRUE)) {
#        install.packages("BiocManager")
    BiocManager::install(version = "3.10")
    if (!h %in% rownames(installed.packages())) {
        print(h)
        BiocManager::install(h, clean=TRUE, quiet=TRUE)
    }
#   }
}
invisible(sapply(packages_bio, install_if_missing2))


#if (!require("BiocManager", quietly = TRUE)) {
#install.packages("BiocManager")
install.packages(c("localpkgs/rjson_0.2.15.tar.gz", "localpkgs/GetoptLong.tar.gz", "localpkgs/ComplexHeatmap.tar.gz", "localpkgs/DEP.tar.gz", "localpkgs/proBatch.tar.gz"), repos=NULL, type="source")
#}


