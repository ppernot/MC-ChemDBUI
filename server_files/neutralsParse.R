reacScheme = reactiveVal(NULL)

# id = shiny::showNotification(
#   h4(paste0('Processing ',filename)),
#   closeButton = FALSE,
#   duration = 5,
#   type = 'message'
# )


observeEvent(
  input$neutralsParseBtn,
  {
    req(neutralsVersion())
    
    nbReac = 0
    reactants = products = params = type = orig =  list()
    
    filename = 'Titan - Réactions bimoléculaires.csv'
    scheme  = read.csv(
      file = file.path(neutralsSource,neutralsVersion(),filename),
      header = FALSE
    )
    comments = scheme[, ncol(scheme)]
    scheme  = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
    for (i in 1:nrow(scheme)) {
      if(substr(scheme[i,1],1,1)=='#') next
      nbReac = nbReac + 1
      terms = scheme[i, 1:3]
      reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 4:8]
      products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 9:13]
      params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
      params[[nbReac]][6] = 'kooij'
      type[[nbReac]]      = 'kooij'
      orig[[nbReac]]      = filename
    }
    
    filename = 'Titan - Réactions trimoléculaires.csv'
    scheme  = read.csv(
      file = file.path(neutralsSource,neutralsVersion(),filename),
      header = FALSE
    )
    comments = c(comments, scheme[, ncol(scheme)])
    scheme  = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
    for (i in 1:nrow(scheme)) {
      if(substr(scheme[i,1],1,1)=='#') next
      nbReac = nbReac + 1
      terms = scheme[i, 1:3]
      reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 4:8]
      products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 9:24]
      params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
      params[[nbReac]][17]= 'assocMD'
      type[[nbReac]]      = 'assocMD'
      orig[[nbReac]]      = filename
    }
    
    ## Additional bimolecular data from misc. sources
    filename = 'bimol_supp.csv'
    scheme  = read.csv(
      file = file.path(neutralsSource,neutralsVersion(),filename),
      header = FALSE
    )
    comments = c(comments, scheme[, ncol(scheme)])
    scheme  = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
    for (i in 1:nrow(scheme)) {
      if(substr(scheme[i,1],1,1)=='#') next
      nbReac = nbReac + 1
      terms = scheme[i, 1:3]
      reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 4:8]
      products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 9:13]
      params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
      params[[nbReac]][6] = 'kooij'
      type[[nbReac]]      = 'kooij'
      orig[[nbReac]]      = filename
    }
    
    ## Additional bimolecular data from misc. sources
    filename = 'trimol_supp.csv'
    scheme  = read.csv(
      file = file.path(neutralsSource,neutralsVersion(),filename),
      header = FALSE
    )
    comments = c(comments, scheme[, ncol(scheme)])
    scheme  = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
    for (i in 1:nrow(scheme)) {
      if(substr(scheme[i,1],1,1)=='#') next
      nbReac = nbReac + 1
      terms = scheme[i, 1:3]
      reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 4:8]
      products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 9:24]
      params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
      params[[nbReac]][17]= 'assocMD'
      type[[nbReac]]      = 'assocMD'
      orig[[nbReac]]      = filename
    }
    
    ## Additional trimolecular data from Vuitton2019
    ## (With specific parameterization)
    
    filename = 'trimol_VV.csv'
    scheme  = read.csv(
      file = file.path(neutralsSource,neutralsVersion(),filename),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    comments = c(comments, scheme[, ncol(scheme)])
    scheme  = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
    for (i in 1:nrow(scheme)) {
      if(substr(scheme[i,1],1,1)=='#') next
      nbReac = nbReac + 1
      terms = scheme[i, 1:3]
      reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 4:8]
      products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
      terms = scheme[i, 9:24]
      params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
      params[[nbReac]][17]= 'assocVV'
      type[[nbReac]]      = 'assocVV'
      orig[[nbReac]]      = filename
    }
    
    # Auxilliary data
    species = sort(unique(unlist(c(reactants, products))))
    compo = 
      t(
        apply(
          X = as.matrix(species,ncol=1), 
          MARGIN = 1,
          FUN = function(x){get.atoms(x,stoechFilters)}
        )
      )
    colnames(compo)=elements
    mass  = apply(compo,1,massFormula)
    names(mass) = species
    dummyMass   = round(max(mass,na.rm=TRUE)+2)
    mass[dummySpecies] = dummyMass    
    
    if(any(is.na(mass))) {
      sp = names(mass[is.na(mass)])
      id = shiny::showNotification(
        h4(
          paste0(paste0(sp,collapse=', '),
                 ' has/have no mass !')
        ),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
    }
    
    reacScheme(
      list(
        nbReac    = nbReac,
        reactants = reactants,
        products  = products,
        params    = params,
        type      = type,
        orig      =  orig,
        species   = species,
        mass      = mass
      )
    )
  }
)
output$tabScheme = DT::renderDT(
  {
    req(reacScheme())
    
    nbReac    = reacScheme()$nbReac
    reactants = reacScheme()$reactants
    products  = reacScheme()$products
    params    = reacScheme()$params
    type      = reacScheme()$type
    orig      = reacScheme()$orig
    
    formatReac <- function(reactants, products, params, type) {
      react = paste0(reactants, collapse = ' + ')
      prods = paste0(products , collapse = ' + ')
      pars  = paste0(params, collapse = ', ')
      typ   = unlist(type)
      
      return(data.frame(
        Reactants = react,
        Products  = prods,
        Params    = pars,
        Type      = typ
      ))
    }
    
    dat = data.frame(
      Reactants = NA,
      Products = NA,
      Params = NA,
      Type = NA
    )
    if (is.na(input$targetSpecies) |
        input$targetSpecies == "") {
      for (i in 1:nbReac) {
        msg = checkBalance(reactants[[i]],products[[i]],
                           stoechFilters = stoechFilters)
        if(!is.null(msg))
          id = shiny::showNotification(
            h4(msg),
            closeButton = TRUE,
            duration = NULL,
            type = 'error'
          )
        dat = rbind(dat, formatReac(
          reactants[[i]],products[[i]],params[[i]],type[[i]]
        ))
        
      }
      
    } else {
      # Species-specific reaction list
      nOK = 0
      for (i in 1:nbReac)
        if (input$targetSpecies %in% reactants[[i]] ||
            input$targetSpecies %in% products[[i]]) {
          dat = rbind(dat, formatReac(
            reactants[[i]],products[[i]],params[[i]],type[[i]]
          ))
          nOK = nOK +1
        }
      if (nOK == 0)
          id = shiny::showNotification(
            h4(paste0('Species not in scheme: ', input$targetSpecies)),
            closeButton = TRUE,
            duration = NULL,
            type = 'warning'
          )
    }
    
    return(dat[-1, ])
  },
  rownames = FALSE,
  extensions = c('Scroller'),
  options = list(
    dom         = 'Btip',
    deferRender = TRUE,
    scrollY     = 600,
    scroller    = TRUE
  )
)

output$massScheme = DT::renderDT(
  {
    req(reacScheme())
    
    species = reacScheme()$species
    masssp  = reacScheme()$mass

    dat = data.frame(
      Mass = NA,
      Species = NA
    )
    mu = sort(unique(masssp))
    for (i in seq_along(mu)) {
      m = mu[i]
      spm = names(masssp)[masssp==m]
      dat = rbind(
        dat,
        data.frame(
          Mass = mu[i], 
          Species = paste0(spm, collapse = ', ')
        )
      )
    }
    return(dat[-1,])
  },
  rownames = FALSE,
  extensions = c('Scroller'),
  options = list(
    dom         = 'Btip',
    deferRender = TRUE,
    scrollY     = 600,
    scroller    = TRUE
  )
)
