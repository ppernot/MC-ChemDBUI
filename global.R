version = 0.0

# Options ####
Sys.setlocale(category = "LC_NUMERIC", locale = "C")
options(
  # shiny.maxRequestSize=30*1024^2, # Increase max loading size to 30 Mo
  width = 70,
  warn = 0,
  stringsAsFactors = FALSE
)

# Load packages ####
source("R/packages.R")

# Proportions of Side/Main Panels ####
sideWidth  <- 3
mainWidth  <- 12 - sideWidth
plotHeight <- 550
plotWidth  <- 550

# Graphical parameters ####
gPars = ErrViewLib::setgPars('shiny')

# Fine tune graph. params
gPars$cex = 1.5
gPars$mar[3] = 2

# Data paths ####
neutralsSource = file.path('..','MC-ChemDB','Neutrals','Source')

# Load data and functions ####
source('R/massCalc.R')
elements      = unlist(read.csv("data/elements.csv", header = FALSE))
massElem      = CHNOSZ::mass(elements)
dummySpecies  = unlist(read.csv("data/dummySpecies.csv", header = FALSE))
stoechFilters = read.csv("data/stoechFilters.csv", header = FALSE)

