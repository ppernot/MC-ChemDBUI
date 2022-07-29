filterFormula <- function (sp) {
  # Normalize chemical formula for mass estimation
  if(is.null(stoechFilters))
    stop(">>> Need stoechFilters in global Env. !")
  sp1 = sp
  for(i in 1:nrow(stoechFilters))
    sp1 = sub(stoechFilters[i,1], stoechFilters[i,2], sp1)
  return(sp1)
}
calcAtoms = function(formula) {
  # Split formula into chemical composition vector
  compo = matrix(0, nrow = 1, ncol = length(elements))
  colnames(compo) = elements
  if (!(formula %in% c("E","HV"))) { # Mass = 0
    atoms = suppressWarnings(CHNOSZ::i2A(formula))
    compo[1, colnames(atoms)] = atoms
  }
  return(compo[1, ])
}
get.atoms <- function (sp) {
  # Normalize formula and extract composition vector
  # Rq: stoechFilters is a global variable
  sp1 = filterFormula(sp)
  tryCatch(
    calcAtoms(sp1),
    error = function(x)
      rep(NA, length(elements))
  )
}
massFormula = function(sto) {
  # Compute mass from composition vector
  sum(sto * massElem)
}
getMassList = function (species, excludeList = 'Products') {
  # Compute mass for set of formulae
  if (any(species %in% excludeList)) {
    mass = NA
  } else {
    mass = sum(
      sapply(
        X = species,
        FUN = function(x){massFormula(get.atoms(x))}
      ),
      na.rm = TRUE
    )
  }
  return(mass)
}
checkBalance = function(reactants, products) {
  reac = unlist(reactants)
  prod = unlist(products)
  massReacs = getMassList(reac, excludeList = dummySpecies)
  massFrags = getMassList(prod, excludeList = dummySpecies)
  reacTags = paste0(
    paste0(reac,collapse = ' + '),
    ' --> ',
    paste0(prod,collapse = ' + ')
  )
  msg = NULL
  if(!is.na(massReacs) & !is.na(massFrags))
    if(abs(massFrags-massReacs) > 0.01) 
      msg = paste0(
        'Unbalanced masses: \n',
        reacTags,'\n',
        massReacs, ' /= ', massFrags
      ) 
  return(msg)
}
numElec <- function(sto) {
  # Calculate electron number of composition
  sum(sto * numElecElem)
}
spCharge = function(species) {
  charge = rep(0,length(species))
  ions   = grepl("\\+$",species)
  charge[ions] = 1
  charge[which(species == 'E')] = -1
  return(charge)
}
selectSpecies <- function(species, categs) {
  compo   = t(apply(as.matrix(species,ncol=1),1,get.atoms)) # Enforce matrix
  colnames(compo) = elements
  
  # Remove dummy species
  sel0 = ! species %in% dummySpecies
  
  if ('all' %in% categs) {
    return ( sel0)
    
  } else {
    
    # Charge
    charge = rep(0,length(species))
    ions   = grepl("\\+$",species)
    charge[ions] = 1
    neus = !ions
    sel  = rep(FALSE,length(species))
    for (cl in categs) {
      if (cl == "neutrals") {
        sel = sel | neus
      } else if (cl == "ions") {
        sel = sel | ions
      }
    }
    selCharge = sel
    
    # Radicals
    radic = (apply(compo, 1, numElec) - charge) %% 2
    for (cl in categs) {
      if (cl == "radicals") {
        sel = radic
      }
    }
    selRadic = sel
    
    # Composition Elementale
    azot = grepl("N",species)
    oxy  = grepl("O",species)
    sel  = rep(FALSE,length(species))
    for (cl in categs) {
      if (cl == "hydrocarbons") {
        sel = sel | (!azot & !oxy)
      } else if (cl == "N-bearing") {
        sel = sel | azot
      } else if (cl == "O-bearing") {
        sel = sel | oxy
      }
    }
    selElem = sel
    
    # Heavy elements
    nHeavy = nbHeavyAtoms(species)
    sel    = rep(FALSE,length(species))
    for (cl in categs) {
      if (cl == "C0") {
        sel = sel | nHeavy == 0
      } else if (cl == "C1") {
        sel = sel | nHeavy == 1
      } else if (cl == "C2") {
        sel = sel | nHeavy == 2
      } else if (cl == "C3") {
        sel = sel | nHeavy == 3
      } else if (cl == "C4") {
        sel = sel | nHeavy == 4
      } else if (cl == "C5") {
        sel = sel | nHeavy == 5
      } else if (cl == "C6") {
        sel = sel | nHeavy == 6
      } else if (cl == "Cmore") {
        sel = sel | nHeavy > 6
      }
    }
    selHeavy = sel
    
    return(sel0 & selCharge & selRadic & selElem & selHeavy)
    
  }
  
}
nbHeavyAtoms <- function(spList) {
  # Number of non-hydrogen atoms in formula
  compo = t(sapply(spList, get.atoms))
  nHeavy = rowSums(compo) - compo[, 1]
  return(nHeavy)
}
