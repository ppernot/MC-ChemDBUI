version = 0.2

# Options ####
Sys.setlocale(category = "LC_NUMERIC", locale = "C")
options(
  # shiny.maxRequestSize=30*1024^2, # Increase max loading size to 30 Mo
  width = 70,
  warn = 0,
  stringsAsFactors = FALSE
)

source_ui <- function(...) {
  source(
    file.path("ui_files", ...),
    local = TRUE
  )$value
}

# Load packages ####
source("R/packages.R")

# Proportions of Side/Main Panels ####
sideWidth  <- 3
mainWidth  <- 12 - sideWidth
plotHeight <- 600
plotWidth  <- 550

# Graphical parameters ####
gPars = ErrViewLib::setgPars('shiny')
col2tr =function(x,alpha=80){
  rgb(unlist(t(col2rgb(x))),alpha=alpha,maxColorValue = 255)
}

# Fine tune graph. params
gPars$cex = 1.5
gPars$mar[3] = 2

# Data paths ####
neutralsSource = file.path('..','MC-ChemDB','Neutrals','Source')
neutralsPublic = file.path('..','ChemDBPublic','Neutrals')
ionsSource     = file.path('..','MC-ChemDB','Ions','Source')
ionsTmp        = file.path('..','MC-ChemDB','Ions','Tmp')
ionsPublic     = file.path('..','ChemDBPublic','Ions')

# Load data and functions ####
source('R/massCalc.R')
maxReacts     = 3     # Max number of reactants slots in generated dBases
maxProds      = 4     # Max number of product slots in generated dBases
elTable       = read.csv(file.path('data','elements.csv'), header = FALSE)
elements      = unlist(elTable[1,])
numElecElem   = as.numeric(unlist(elTable[2,]))
massElem      = CHNOSZ::mass(elements)
dummySpecies  = unlist(read.csv(file.path('data','dummySpecies.csv'), header = FALSE))
stoechFilters = read.csv(file.path('data','stoechFilters.csv'), header = FALSE, 
                         allowEscapes = TRUE)
tabNeuFiles   = read.csv(file.path('data','neutralsDBFiles.csv'),header = FALSE)

source('R/rateFormulas.R')

ionsRateParKwdList = c('ALPHA','BETA','GAMMA') 
ionsReacTypes = c('dr','kooij','ionpol1','ionpol2')
source('R/ionsFunctions.R')

# Biblio for ions
bibFile = file.path('..','MC-ChemDB','Doc','refsDR.bib')
bibProc = file.path('data','bib.Rdata')
if(!file.exists(bibProc)) {
  cat('*** Processing .bib file\n')
  bib = bibtex::read.bib(file=bibFile)
  save(bib, file=bibProc)
} else {
  sourceTime = file.info(bibFile)["mtime"]
  bibTime    = file.info(bibProc)["mtime"]
  if(sourceTime > bibTime) {
    cat('*** Processing .bib file\n')
    bib = bibtex::read.bib(file=bibFile)
    save(bib, file=bibProc)    
  } else {
    cat('*** Loading  processed .bib file\n')
    load(bibProc)
  }
}



