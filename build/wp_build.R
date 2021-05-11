#install.packages("devtools")
library(devtools)

pkg=file.path("~/git", "WadingPool")
setwd(pkg)


#use_description()


#### Assembling data ####
usethis::use_testthat()               # unit testing
usethis::use_vignette("vignette")  # vignette documentation
usethis::use_data_raw()               # raw data


#### Building ####
#Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/MacOS/pandoc")
usethis::use_package("assertthat")
usethis::use_package("GenomicRanges")
usethis::use_package("IRanges")
usethis::use_package("S4Vectors")
usethis::use_package("GenomeInfoDb")
usethis::use_package("VariantAnnotation")
usethis::use_package("stats")
usethis::use_package("utils")
usethis::use_package("dplyr")
usethis::use_package("purrr")
usethis::use_package("DescTools")
usethis::use_package("ggplot2")
usethis::use_package("reshape2")
usethis::use_package("gridExtra")
usethis::use_package("depmixS4")


devtools::load_all()
devtools::document(pkg)
devtools::check(pkg)
devtools::build_vignettes(pkg)

devtools::build(pkg)
devtools::install(pkg)
detach("package:WadingPool", unload=TRUE)
library(WadingPool)
# devtools::install_github("quevedor2/aneuploidy_score")
devtools::install_github("quevedor2/WadingPool", ref = "dev")

# library(AneuploidyScore)
# AneuploidyScore::listData()
