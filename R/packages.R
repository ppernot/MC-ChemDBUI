
# Local install ----

# # CRAN Libraries ----
# libs <- c(
#   "devtools", "boot",
#   "shiny", "shinyFiles","shinycssloaders",
#   "data.table","CHNOSZ"
#
# )
# for (lib in libs) {
#   if (!require(lib, character.only = TRUE, quietly = TRUE))
#     install.packages( lib,
#       # repos = "https://cran.univ-paris1.fr",
#       dependencies = FALSE
#     )
#   library(lib, character.only = TRUE, quietly = TRUE)
# }
# # ## Other libraries ----
# devtools::install_github("ppernot/ErrViewLib")

# Cloud deployment ----

library("devtools")
library("shiny")
library("shinythemes")
library("shinyFiles")
library("shinyAce")
library("shinycssloaders")
library("DT")
library("boot")
library("ErrViewLib")
library("data.table")
library("CHNOSZ")
library("stringr")
library("ape")
library("bibtex")
