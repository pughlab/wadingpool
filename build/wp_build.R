#install.packages("devtools")
library(devtools)

pkg=file.path("~/git", "WadingPool")
setwd(pkg)


#use_description()

#### Assembling data ####
usethis::use_data_raw()
devtools::load_all()

#### Building ####
#Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/MacOS/pandoc")

usethis::use_package("GenomicRanges")
usethis::use_package("IRanges")
usethis::use_package("S4Vectors")
usethis::use_package("GenomeInfoDb")
usethis::use_package("VariantAnnotation")
usethis::use_package("stats")
usethis::use_package("utils")


devtools::document(pkg)
devtools::check(pkg)
devtools::build_vignettes(pkg)

devtools::build(pkg)
devtools::install(pkg)
# devtools::install_github("quevedor2/aneuploidy_score")
devtools::install_github("quevedor2/WadingPool", ref = "dev")

# library(AneuploidyScore)
# AneuploidyScore::listData()
