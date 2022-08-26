## Install R packages
install.packages('shinyFeedback', dependencies = TRUE, repos='http://cran.rstudio.com/')

#install.packages('textreadr', dependencies = TRUE, repos='http://cran.rstudio.com/')
#install.packages('pdftools', dependencies = TRUE, repos='http://cran.rstudio.com/')

#install.packages('RSQLite', version="2.2.7", dependencies = FALSE, repos='http://cran.rstudio.com/')
#install.packages('RSQLite', version="2.2.7", dependencies = FALSE,repos='http://cran.us.r-project.org', type="binary")
#remove.packages("BiocManager")
#install.packages("DBI", version="1.1.2", dependencies=FALSE,repos='http://cran.rstudio.com/')
#install.packages("rlang", version="0.4.11", dependencies=FALSE,repos='http://cran.rstudio.com/')
#install.packages("bit", version="4.0.4", dependencies=FALSE,repos='http://cran.rstudio.com/')
#install.packages("bit64", version="4.0.5", dependencies=FALSE,repos='http://cran.rstudio.com/')
#install.packages("blob", version="1.2.2", dependencies=FALSE,repos='http://cran.rstudio.com/')
#install.packages("Rcpp", version="1.0.6", dependencies=FALSE,repos='http://cran.rstudio.com/')
#install.packages("RSQLite", version="2.2.15", dependencies=FALSE,repos='http://cran.rstudio.com/')
#install.packages("rmarkdown", version="2.11", dependencies=FALSE,repos='http://cran.rstudio.com/')
#install.packages("BiocManager", version="3.10", dependencies=FALSE,repos='http://cran.rstudio.com/')
#BiocManager::install("S4Vectors", version = "0.24.4")
#BiocManager::install("IRanges", version = "2.20.2")
#BiocManager::install("BioGenerics", version = "0.32.0")
#BiocManager::install("Biobase", version = "2.46.0")
#BiocManager::install("AnnotationDbi", version = "1.48.0")




#packageurl <- "other/RSQLite_2.2.7.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

#packageurl <- "other/DEP.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

#packageurl <- "other/proBatch.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

#packages <- c("shiny", "shinyjs","shinyalert", "shinyvalidate", "shinythemes", "shinycssloaders", "DT", "tidyr", "textreadr", "optparse", "dplyr", "tibble", "ggplot2", "assertthat", "readr", "purrr", "tidyverse", "ggpubr", "RColorBrewer", "gtable", "grid", "gridExtra", "sqldf", "tools", "missRanger", "plotly","rsconnect")
install_if_missing <- function(p) {
    if (!p %in% rownames(installed.packages())) {
        install.packages(p)
    }
}
#invisible(sapply(packages, install_if_missing))

#packages_bio <- c("DEP", "proBatch", "pcaMethods", "ComplexHeatmap", "SummarizedExperiment")
#remove.packages(packages_bio)
#install.packages('DEP', version = "0.9.1", dependencies = TRUE, repos='http://cran.rstudio.com/') 
install_if_missing2 <- function(h) {
    if (!require("BiocManager", quietly = TRUE)) {}
        #install.packages("BiocManager")
        if (!h %in% rownames(installed.packages())) {
            BiocManager::install(h)
        }
    }

#invisible(sapply(packages_bio, install_if_missing2))

