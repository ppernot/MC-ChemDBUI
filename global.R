version = 0.3

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
sideWidth  = 3
mainWidth  = 12 - sideWidth
plotHeight = 600
plotWidth  = 550

# Max. nb. of options in drop-down menus 
# (upper limit to nb. of reactions)
maxOptions = 10000 

# Graphical parameters ####
gPars = ErrViewLib::setgPars('shiny')
col2tr =function(x,alpha=80){
  rgb(unlist(t(col2rgb(x))),alpha=alpha,maxColorValue = 255)
}

# Fine tune graph. params
gPars$cex = 1.5
gPars$mar[3] = 2

# Data paths ####
neutralsSource = file.path('..','MC-ChemDB','Neutrals')
neutralsPublic = file.path('..','ChemDBPublic','Neutrals')
ionsSource     = file.path('..','MC-ChemDB','Ions')
ionsTmp        = file.path('..','MC-ChemDB','Tmp')
ionsPublic     = file.path('..','ChemDBPublic','Ions')
photoSource    = file.path('..','MC-ChemDB','PhotoProcs')
photoPublic    = file.path('..','ChemDBPublic','PhotoProcs')

# Load data and functions ####
source('R/massCalc.R')
maxReacts     = 3     # Max number of reactants slots in generated dBases
maxProds      = 4     # Max number of product slots in generated dBases
elTable       = read.csv(file.path('data','elements.csv'), header = FALSE)
elements      = unlist(elTable[1,])
numElecElem   = as.numeric(unlist(elTable[2,]))
massElem      = CHNOSZ::mass(elements)
dummySpecies  = unlist(read.csv(file.path('data','dummySpecies.csv'), 
                                header = FALSE))
stoechFilters = read.csv(file.path('data','stoechFilters.csv'), 
                         header = FALSE, 
                         allowEscapes = TRUE)

neutralsRateParKwdList = c('A1','B1','C1','F1','G1',
                           'A2','B2','C2','F2','G2',
                           'A3','B3','C3','F3','G3',
                           'FC') 
neutralsReacTypes  = c('kooij','assocMD','assocVV')
source('R/rateFormulas.R')

ionsRateParKwdList = c('ALPHA','BETA','GAMMA') 
ionsReacTypes = c('dr','kooij','ionpol1','ionpol2')
source('R/ionsFunctions.R')

photoKwdList = c('CHANNEL','XS_SOURCE','XS_F','BR_SOURCE','OUTPUT')
photoXSSources = c('Leiden','SWRI','Hebrard')
photoBRSources = c('SWRI','Plessis')
photoXSResolutions = c(1,0.1) # nm
photoDefaultuF = 1.2
photoRuBRN  = 0.2  # relative uncertainty for Neutral branching ratios
photoRuBRI  = 0.2  # relative uncertainty for Ionic branching ratios
photoRuBRNI = 0.03 # relative uncertainty on Ionic vs. Neutral channels
photoEps    = 5e-3 # threshold fo zero in compositions
source('R/photoFunctions.R')

# Bibliography
bibFile = file.path('..','MC-ChemDB','Doc','refsDR.bib')
bibProc = file.path('..','ChemDBPublic','bib.Rdata')
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



