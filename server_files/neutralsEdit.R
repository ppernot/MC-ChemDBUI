neutralsReacs         = shiny::reactiveVal()
neutralsRateMask      = shiny::reactiveVal()
neutralsSimulSamples  = shiny::reactiveVal()
neutralsReacsFiltered = shiny::reactiveVal()
neutralsReacID        = shiny::reactiveVal()

# Manage reacs list ####
observeEvent(
  input$neutralsReacSelInit,
  {
    shiny::updateTextInput(inputId = "neutralsReacSel", value = "")
    shiny::updateRadioButtons(inputId = "neutralsReacSelKind", selected = "Both")
  })

observe({
  req(neutralsDB())
  reacs = neutralsDB()$REACTANTS
  prods = neutralsDB()$PRODUCTS
  tag   = neutralsDB()$TAG       # "tags" is reserved by shiny...

  if(input$neutralsReacSel != "" & !is.null(input$neutralsReacSel)){
    if(input$neutralsReacSelKind == "Reactant") {
      sel = sapply(
        reacs,
        FUN = function(x) input$neutralsReacSel %in% getSpecies(x)
      )
    } else if(input$neutralsReacSelKind == "Product"){
      sel = sapply(
        prods,
        FUN = function(x) input$neutralsReacSel %in% getSpecies(x)
      )
    } else if(input$neutralsReacSelKind == "Both"){
      sel = apply(
        cbind(reacs,prods),
        1,
        FUN = function(x) 
          input$neutralsReacSel %in% getSpecies(paste0(x[1],' + ',x[2]))
      ) 
    }
    if(sum(sel) == 0) {
      id = shiny::showNotification(
        strong(
          paste0(
            'Unknown species: ',input$neutralsReacSel,
            ' in context: ',input$neutralsReacSelKind)
        ),
        closeButton = TRUE,
        duration = 5,
        type = 'error'
      )
      tag = NULL
    } else {
      tag = tag[sel]
    }
  }
  req(tag)
  
  nums = 1:length(tag)  
  
  if(max(nums) > maxOptions)
    id = shiny::showNotification(
      strong(
        paste0(
          'Nb. of reactions (',max(nums),') exceeds limit ',
          'of drop-down menu (',maxOptions,'). Increase variable ',
          'maxOptions in MC-ChemDBUI/global.R.'
        )
      ),
      closeButton = TRUE,
      duration = 5,
      type = 'warning'
    )

  # Update selector
  names(nums) = tag
  shiny::updateSelectizeInput(
    session  = session,
    inputId  = "neutralsReaction", 
    choices  = nums,
    selected = 1,
    server   = TRUE
  )
  
  neutralsReacsFiltered(nums) # Partial list used by other parts
})

observeEvent(
  input$neutralsMinus,
  {
    iReac = as.numeric(input$neutralsReaction)
    if(iReac > 1) {
      iReac = iReac - 1
      shiny::updateSelectizeInput(
        session = session,
        "neutralsReaction",
        selected = iReac
      )      
    }
  })

observeEvent(
  input$neutralsPlus,
  {
    req(neutralsReacsFiltered())
    
    reacs = neutralsReacsFiltered()
    iReac = as.numeric(input$neutralsReaction)
    if(iReac < length(reacs)) {
      iReac = iReac + 1
      shiny::updateSelectizeInput(
        session = session,
        "neutralsReaction",
        selected = iReac
      )      
    }
  })

# Parse data file ####
shiny::observe({
  req(neutralsDB())
  req(input$neutralsReaction != "0")
  
  # Entry in table
  isel  = as.numeric(input$neutralsReaction)
  if(input$neutralsReacSel != ""){
    tag = names(neutralsReacsFiltered())[isel]
  } else {
    tag = neutralsDB()$TAG[isel]
  }
  req(tag)
  
  # Absolute index of selected reac
  iReac = which(neutralsDB()$TAG == tag)
  
  neutralsReacID(neutralsDB()$ID[iReac]) # Set ID for save
  
  # (Re-)Init samples for plots
  neutralsSimulSamples(NULL)
  
  # Format data for masks
  mask = list()
  reactants = getSpecies(neutralsDB()$REACTANTS[iReac])
  mask[['REACTANTS']] = reactants
  # massReactants = getMassList(reactants, excludeList = dummySpecies)
  
  products = getSpecies(neutralsDB()$PRODUCTS[iReac])
  mask[['PRODUCTS']] = products
  # massProducts = getMassList(products, excludeList = dummySpecies)
  
  mask[['TYPE']] = neutralsDB()$TYPE[iReac]
  
  for (kwd in neutralsRateParKwdList)
    mask[[kwd]] = neutralsDB()[[kwd]][iReac]
  
  mask[['REFS']] = neutralsDB()$REFS[iReac]
  mask[['RQ']] = neutralsDB()$COMMENTS[iReac]
  mask[['TIMESTAMP']] = neutralsDB()$TIMESTAMP[iReac]
  neutralsRateMask(mask)

})

output$neutralsRateMask = shiny::renderUI({
  req(neutralsRateMask())
  
  mask = neutralsRateMask()
  type = tolower(mask[['TYPE']])
  
  list(
    br(),
    fluidRow(
      # All built on width 11 instead of 12...
      column(
        4,
        shiny::textInput(
          "neutralsReacReactants",
          "Reactants",
          value = paste0(mask[['REACTANTS']],collapse = ' + ')
        )
      ),
      column(
        4,
        shiny::textInput(
          "neutralsReacProducts",
          "Products",
          value = paste0(mask[['PRODUCTS']],collapse = ' + ')
        )
      ),
      column(
        3,
        shiny::selectInput(
          "neutralsReacTYPE",
          "Reaction type",
          neutralsReacTypes,
          selected = type
        )
      )
    ),
    fluidRow(
      column(
        1,
        strong(ifelse(type == 'assocvv','kInf','k0'))
      ),
      column(
        2,
        shiny::textInput("neutralsReacA1", "A1", value = mask[['A1']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacB1", "B1", value = mask[['B1']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacC1", "C1", value = mask[['C1']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacF1", "F1", value = mask[['F1']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacG1", "G1", value = mask[['G1']])
      )
    ),
    fluidRow(
      column(
        1,
        strong(ifelse(type == 'assocvv','k0','kInf'))
      ),
      column(
        2,
        shiny::textInput("neutralsReacA2", "A2", value = mask[['A2']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacB2", "B2", value = mask[['B2']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacC2", "C2", value = mask[['C2']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacF2", "F2", value = mask[['F2']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacG2", "G2", value = mask[['G2']])
      )
    ),
    fluidRow(
      column(
        1,
        strong('kr')
      ),
      column(
        2,
        shiny::textInput("neutralsReacA3", "A3", value = mask[['A3']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacB3", "B3", value = mask[['B3']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacC3", "C3", value = mask[['C3']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacF3", "F3", value = mask[['F3']])
      ),
      column(
        2,
        shiny::textInput("neutralsReacG3", "G3", value = mask[['G3']])
      )
    ),
    fluidRow(
      column(
        8,
        shiny::textAreaInput(
          "neutralsReacRQ",
          "Comments",
          value = mask[['RQ']],
          height = '200px'
        )
      ),
      column(
        3,
        shiny::textInput(
          "neutralsReacFC", 
          "Fc", 
          value = mask[['FC']]
        ),
        shiny::textInput(
          "neutralsReacREF",
          "References",
          value = mask[['REFS']]
        ),
        shiny::textInput(
          "neutralsReacTIMESTAMP",
          "Time Stamp",
          value = mask[['TIMESTAMP']]
        )
      )
    )
  )
  
})
outputOptions(output, "neutralsRateMask", suspendWhenHidden = FALSE)


# Save ####
observeEvent(
  input$neutralsParseSave,
  {
    req(neutralsDB())

    reactants = rep(NA, 3)
    elts = trimws(unlist(strsplit(input$neutralsReacReactants,'+',fixed = TRUE)))
    reactants[1:length(elts)] = elts
   
    products = rep(NA, 5)
    elts = trimws(unlist(strsplit(input$neutralsReacProducts,'+',fixed = TRUE)))
    products[1:length(elts)] = elts

    msg = checkBalance(reactants,products)
    if(!is.null(msg))
      id = shiny::showNotification(
        strong(msg),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
    req(is.null(msg))
    
    if(input$neutralsParseComment)
      reactants[1] = paste0('#',reactants[1])
    
    ratePars = rep(NA,length(neutralsRateParKwdList))
    names(ratePars) = neutralsRateParKwdList
    for (kwd in neutralsRateParKwdList) 
      ratePars[kwd] = as.numeric(input[[paste0('neutralsReac',kwd)]])

    id = neutralsReacID()

    line = c(
      reactants,
      products,
      paste0(trimws(input$neutralsReacTYPE)),
      ratePars,
      paste0(input$neutralsReacREF),
      paste0(input$neutralsReacRQ),
      paste0(Sys.time()),
      id
    )
    line = matrix(line,nrow=1)
    line = capture.output(
      write.table(line,sep=';',row.names = FALSE, col.names = FALSE)
    )
    # line = paste0(line,collapse = ";")
    
    data = neutralsDB()
    tag  = makeTag(
      c(id,
        input$neutralsReacReactants,
        input$neutralsReacProducts)
    )
    iReac = which(data$TAG == tag)
    
    # Update editor's content
    dataEditor = neutralsEditDBText()
    dl = length(dataEditor)
    
    if(length(iReac) != 0) {
      # Tag exists: replace in situ
      dataEditor[id+1] = line   # Header counts as first line...
    } else {
      # Tag does not exist: set new ID
      line = sub(paste0('^',id),paste0(dl),line)
      dataEditor = c(dataEditor,line)
    }
    
    neutralsEditDBText(dataEditor)
    
  }
)



# Sampling ####
topow = function(x,p) {
  if(x==0)
    return(0.0)
  else
    return(x^p)
}
oneSampleNeutralsPars = function(i, pars, type) {
  p = pars = as.numeric(pars)
  rmax = ifelse(type == 'kooij', 1, 3)
  for (r in 1:rmax) {
    rnd = ifelse(
      i==0, 
      0, 
      truncnorm::rtruncnorm(1,-truncFactor,truncFactor,0,1)) 
    i0 = (r-1)*5
    p[i0+4] = topow(pars[i0+4], rnd)
    p[i0+5] = pars[i0+5] * rnd
  }
  return(p)
}
sampleNeutralsPars = function(nMC, pars, type) {
  tab = matrix(NA, nrow = nMC+1, ncol = length(pars))
  for (i in 0:nMC)
    tab[i+1, ] = oneSampleNeutralsPars(i, pars, type)
  return(tab)
}
neutralsParsSampling = reactive({
    req(neutralsRateMask())
    req(input$neutralsReacTIMESTAMP) # Ensure that dynamic interface is built

    nMC = as.numeric(input$neutralsSimulateSize)

    # Rate parameters
    pars = rep(0,length(neutralsRateParKwdList))
    names(pars) = neutralsRateParKwdList
    for (kwd in neutralsRateParKwdList)
      pars[kwd] = as.numeric(input[[paste0('neutralsReac',kwd)]])

    type = tolower(input$neutralsReacTYPE)

    sampleRateParams = sampleNeutralsPars(nMC, pars, type)
    colnames(sampleRateParams) = neutralsRateParKwdList

    return(sampleRateParams)

  }
)

# Plot rates ####
output$plotNeutralsRate = renderPlot({
  req(neutralsRateMask())
  
  sample = neutralsParsSampling()
  
  T0 = as.numeric(input$T0Plot)
  tRange = seq(input$neutralsTempRangePlot[1],
               input$neutralsTempRangePlot[2],5)
  M0  = 10^as.numeric(input$M0Plot)
  mRange = 10^seq(input$densRangePlot[1],input$densRangePlot[2],0.5)

  type = tolower(input$neutralsReacTYPE)
  
  krateT = matrix(NA, ncol = nrow(sample),nrow= length(tRange))
  krateM = matrix(NA, ncol = nrow(sample),nrow= length(mRange))
  
  for (i in 1:nrow(sample)) {
    pars = sample[i,]
    # fixed M, T varies
    krateT[,i] = switch(
      type,
      kooij     = k_kooij(pars, tempRange = tRange),
      assocmd   = k_assocmd(pars, tempRange = tRange, M0),
      assocvv   = k_assocvv(pars, tempRange = tRange, M0),
      assoc0    = k_assoc0(pars, tempRange = tRange, M0),
      assoctroe = k_assoctroe(pars, tempRange = tRange, M0),
      rep(0,length(tRange))
    )
    # fixed T, M varies
    krateM[,i] = switch(
      type,
      kooij     = k_kooij(pars, tempRange = T0),
      assocmd   = k_assocmd(pars, tempRange = T0, mRange),
      assocvv   = k_assocvv(pars, tempRange = T0, mRange),
      assoc0    = k_assoc0(pars, tempRange = T0, mRange),
      assoctroe = k_assoctroe(pars, tempRange = T0, mRange),
      rep(0,length(mRange))
    )
  }
  
  # Plot 
  par(
    mfrow = c(2, 1),
    mar = c(3,3,1,1),
    mgp = gPars$mgp,
    tcl = gPars$tcl,
    lwd = gPars$lwd,
    pty = 'm',
    cex = 1
  )
  
  tempRange = tRange
  if (diff(range(krateT)) != 0) {
    matplot(
      tempRange,
      krateT,
      type = 'l',
      lty = 1,
      col = gPars$cols_tr[6],
      lwd = 1.5*gPars$lwd,
      log = 'y',
      xlab = 'T [K]',
      ylab = 'Rate constant [cm^3.s^-1]',
      main = ''
    )
    grid()
    lines(tempRange, krateT[, 1], col = gPars$cols[2],lwd = 1.5*gPars$lwd)
    legend('bottomright',title = paste0('M = ',M0,'cm^-3'), legend = NA, bty='n')
    box()
  }
  
  # P-dep
  if (diff(range(krateM)) != 0) {
    matplot(
      mRange,
      krateM,
      type = 'l',
      lty = 1,
      col = gPars$cols_tr[6],
      lwd = 1.5*gPars$lwd,
      log = 'xy',
      xlab = 'M [cm^-3]',
      ylab = 'Rate constant [cm^3.s^-1]',
      main = ''
    )
    lines(mRange, krateM[, 1], col = gPars$cols[2], lwd = 1.5*gPars$lwd)
    grid()
    legend('bottomright',title = paste0('T = ',T0,' K'), legend = NA, bty='n')
    box()
  }  
  
},
height = plotHeight)

# Biblio ####
output$neutralsBiblio = shiny::renderUI({
  req(neutralsRateMask())
  req(input$neutralsReacREF)
  
  k = input$neutralsReacREF
  keys = unlist(str_split(k,';'))
  keys = trimws(sort(unique(keys)))

  refs = '<H4>References</H4><DL>'
  if(length(keys) != 0) {
    for (key in keys)
      refs = paste0(
        refs,
        '<DT>[',key,']</DT><DD>',
        format(bib[key],style="html"),
        '</DD>')
  }
  refs = paste0(refs,'</DL>')
  list(
    HTML(refs)
  )
})

# output$checkSpecies <- renderText({ 
#   req(input$ace_cursor)
#   sp = input$ace_selection
#   compo = get.atoms(sp, stoechFilters = stoechFilters)
#   names(compo) = elements
#   mass  = massFormula(compo)
#   paste0(
#     "Selection: \"", sp, "\"\n",
#     "Mass= ", mass
#   )
# })

# shiny::observeEvent(
#   input$neuSave,
#   {
#     req(neutralsCopyVersion())
#     req(neutralsFile())
#     
#     # Make copy of source directory
#     neutralsOrigDir = file.path(neutralsSource,neutralsOrigVersion())
#     neutralsCopyDir = file.path(neutralsSource,neutralsCopyVersion())
#     if(!dir.exists(neutralsCopyDir)) {
#       dir.create(neutralsCopyDir)
#       files = list.files(path = neutralsOrigDir)
#       for(file in files)
#         file.copy(
#           from = file.path(neutralsOrigDir,file),
#           to   = file.path(neutralsCopyDir,file)
#         )
#       id = shiny::showNotification(
#         h4(paste0('Created new version: ', neutralsCopyVersion())),
#         closeButton = FALSE,
#         duration = 5
#       )
#     }
#     
#     # Save modified file to target version
#     data = isolate(input$ace)
#     writeLines(
#       data,
#       con = file.path(neutralsCopyDir,neutralsFile())
#     )
#     id = shiny::showNotification(
#       h4(paste0('Saved file: ', neutralsFile())),
#       closeButton = FALSE,
#       duration = 5
#     )
#     
#   }
# )
