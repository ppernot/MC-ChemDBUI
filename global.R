version = '1.2'

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
# (upper limit to nb. of reactions in db)
maxOptions = 10000 

# Graphical parameters ####
gPars = ErrViewLib::setgPars('shiny')
col2tr =function(x,alpha=80){
  rgb(unlist(t(col2rgb(x))),alpha=alpha,maxColorValue = 255)
}

# Fine tune graph. params
gPars$cex = 1.5
gPars$mar[3] = 2

# Truncate standard normal random draws to +/- truncFactor
# to avoid outlying samples
truncFactor = 3

# Data paths ####
neutralsSource = file.path('..','MC-ChemDB','Neutrals')
neutralsPublic = file.path('..','ChemDBPublic','Neutrals')
ionsSource     = file.path('..','MC-ChemDB','Ions')
ionsTmp        = file.path('..','MC-ChemDB','Tmp')
ionsPublic     = file.path('..','ChemDBPublic','Ions')
photoSource    = file.path('..','MC-ChemDB','PhotoProcs')
photoPublic    = file.path('..','ChemDBPublic','PhotoProcs')
beamSource     = file.path('..','MC-ChemDB','BeamSpectrumFiles')
beamPublic     = file.path('..','ChemDBPublic')

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

# Neutrals ####
neutralsRateParKwdList = c('A1','B1','C1','F1','G1',
                           'A2','B2','C2','F2','G2',
                           'A3','B3','C3','F3','G3',
                           'FC')

# For neutrals, each reaction type (X) should have a corresponding k_X() 
# function in rateFormulas.R
neutralsReacTypes  = c('kooij','assocmd','assocvv','assoc0','assoctroe')
source('R/rateFormulas.R')
for (type in neutralsReacTypes)
  if(!exists(paste0('k_',type)))
    stop(paste0('Missing function: ',paste0('k_',type)))

# Ions ####
# For ions, the rate functions are hard coded in scripts
# TBD: create functions in rateFormulas.R
ionsRateParKwdList = c('ALPHA','BETA','GAMMA') 
ionsReacTypes = c('dr','kooij','ionpol1','ionpol2')
source('R/ionsFunctions.R')

# Photo-processes ####
photoKwdList = c('CHANNEL','XS_SOURCE','XS_F','BR_SOURCE')
photoXSSources = c('leiden','swri','hebrard','vulcan')
photoBRSources = c('swri','vulcan','plessis')
photoXSResolutions = c(1,0.1) # nm
photoDefaultuF = 1.2
photoRuBRN  = 0.20 # relative uncertainty for Neutral branching ratios
photoRuBRI  = 0.10 # relative uncertainty for Ionic branching ratios
photoRuBRNI = 0.03 # relative uncertainty on Ionic vs. Neutral channels
photoEps    = 5e-3 # threshold for zero-rounding in compositions
source('R/photoFunctions.R')

# Bibliography ####
bibFile = file.path('..','MC-ChemDB','Doc','refsDR.bib')
bibProc = file.path('..','ChemDBPublic','bib.Rdata')
if(!file.exists(bibProc)) {
  cat('*** Processing .bib file\n')
  bib = bibtex::read.bib(file = bibFile)
  save(bib, file = bibProc)
  cat('*** Processed .bib file saved to ChemDBPublic\n')
} else {
  sourceTime = file.info(bibFile)["mtime"]
  bibTime    = file.info(bibProc)["mtime"]
  if (sourceTime > bibTime) {
    cat('*** Processing .bib file\n')
    bib = bibtex::read.bib(file = bibFile)
    save(bib, file = bibProc)
    cat('*** Processed .bib file saved to ChemDBPublic\n')
  } else {
    cat('*** Loading processed .bib file\n')
    load(bibProc)
    cat('*** Processed .bib loaded from ChemDBPublic\n')
  }
}


