
# Local install ----

# # CRAN Libraries ----
# libs <- c(
#   "devtools", "boot",
#   "shiny", "shinyFiles","shinycssloaders",
#   "data.table"
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
# library(lib,character.only = TRUE)

# Cloud deployment ----

library("devtools")
library("shiny")
library("shinyFiles")
library("shinyAce")
library("shinycssloaders")
library("boot")
library("ErrViewLib")
library("data.table")
