## Install R packages
#packages <- c("shiny", "shinyjs", "shinyalert", "shinyvalidate", "shinythemes", "shinycssloaders", "rsconnect", "DT", "tidyr", "textreadr")
#packages <- c("optparse", "dplyr", "tibble", "ggplot2", "assertthat", "readr", "purrr", "tidyverse")
#packages <- c("ggpubr", "RColorBrewer", "gtable", "grid", "gridExtra", "sqldf", "tools", "missRanger", "plotly")
#install_if_missing <- function(p) {
#    print(p)
#    if (!p %in% rownames(installed.packages())) {
#        install.packages(p, clean=TRUE, quiet=TRUE)
#    }
#}
#invisible(sapply(packages, install_if_missing))

#packages_bio <- c("pcaMethods")
packages_bio <- c("DEP", "proBatch", "ComplexHeatmap", "SummarizedExperiment")
install_if_missing2 <- function(h) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
        if (!h %in% rownames(installed.packages())) {
            print(h)
            BiocManager::install(h, clean=TRUE, quiet=TRUE)
        }
    }
}
invisible(sapply(packages_bio, install_if_missing2))

#install.packages(c("localpkgs/DEP.tar.gz", "localpkgs/proBatch.tar.gz"), repos=NULL, type="source")

