filterFormula <- function (sp, stoechFilters = NULL) {
  # Normalize chemical formula for mass estimation
  if(is.null(stoechFilters))
    stop(">>> Need stoechFilters !")
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
get.atoms <- function (sp, stoechFilters = NULL) {
  # Normalize formula and extract composition vector
  sp1 = filterFormula(sp, stoechFilters)
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
getMassList = function (species,excludeList = 'Products', stoechFilters = NULL) {
  # Compute mass for set of formulae
  if (any(species %in% excludeList)) {
    mass = NA
  } else {
    mass = sum(
      sapply(
        X = species,
        FUN = function(x){massFormula(get.atoms(x,stoechFilters))}
      ),
      na.rm = TRUE
    )
  }
  return(mass)
}
checkBalance = function(reactants, products, stoechFilters = NULL) {
  reac = unlist(reactants)
  prod = unlist(products)
  massReacs = getMassList(reac, excludeList = dummySpecies, 
                          stoechFilters = stoechFilters)
  massFrags = getMassList(prod, excludeList = dummySpecies, 
                          stoechFilters = stoechFilters)
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
