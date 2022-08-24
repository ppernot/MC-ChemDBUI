neutralsDB            = shiny::reactiveVal()
neutralsReacs         = shiny::reactiveVal()
neutralsRateMask      = shiny::reactiveVal()
neutralsSimulSamples  = shiny::reactiveVal()
neutralsReacsFiltered = shiny::reactiveVal()
neutralsReacID        = shiny::reactiveVal()

makeTag = function(x) 
  trimws(paste0(x[1],': ',x[2],' -> ',x[3]))

# Manage reacs list ####
shiny::observe({
  req(input$aceNeutralsDB)

  # Get all data in (do not use comment.char here because it mixes up with id.)
  data = read.table(
    header = TRUE, 
    text = input$aceNeutralsDB, 
    sep=';',
    comment.char = "")
  
  # Filter out comment lines
  sel = substr(data$R1,1,1) != '#'
  data = data[sel,]
  
  data[['REACTANTS']] = apply(data[,paste0("R",1:3)], 1, 
                              function(x) paste0(x[!is.na(x) & x != ""],
                                                 collapse = ' + '))
  data[['PRODUCTS']]  = apply(data[,paste0("P",1:5)], 1, 
                              function(x) paste0(x[!is.na(x) & x != ""],
                                                 collapse = ' + '))
  # Build unique tag
  data[['TAG']]       = apply(data[,c("ID","REACTANTS","PRODUCTS")], 1, 
                              makeTag)
  
  neutralsDB(data)
})

observeEvent(
  input$neutralsReacSelInit,
  {
    shiny::updateTextInput(inputId = "neutralsReacSel", value = "")
    shiny::updateRadioButtons(inputId = "neutralsReacSelKind", selected = "Both")
  })

output$selNeutralsReac = shiny::renderUI({
  req(neutralsDB())
  reacs = neutralsDB()$REACTANTS
  prods = neutralsDB()$PRODUCTS
  tag   = neutralsDB()$TAG       # "tags" is reserved by shiny...

  if(input$neutralsReacSel != ""){
    if(input$neutralsReacSelKind == "Reactant") {
      sel = grepl(input$neutralsReacSel,reacs) 
    } else if(input$neutralsReacSelKind == "Product"){
      sel = grepl(input$neutralsReacSel,prods)
    } else if(input$neutralsReacSelKind == "Both"){
      sel = grepl(input$neutralsReacSel,tag)    
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
  names(nums) = tag
  neutralsReacsFiltered(nums) # Partial list used by other parts

  list(
    fluidRow(
      column(
        8,
        shiny::selectInput(
          "neutralsReaction",
          "Reactions",
          choices  = nums,
          selected = 1
        )
      ),
      column(
        4,
        fluidRow(
          column(
            6,
            actionButton(
              "neutralsMinus",
              "",
              icon = icon('angle-down',verify_fa = FALSE)
            ),
            tags$style(
              type='text/css',
              "#neutralsMinus { width:100%; margin-top: 30px;}"
            )
          ),
          column(
            6,
            actionButton(
              "neutralsPlus",
              "",
              icon = icon('angle-up',verify_fa = FALSE)
            ),
            tags$style(
              type='text/css',
              "#neutralsPlus { width:100%; margin-top: 30px;}"
            )
          )
        )
      )
    )
  )
})

observeEvent(
  input$neutralsMinus,
  {
    iReac = as.numeric(input$neutralsReaction)
    if(iReac > 1) {
      iReac = iReac - 1
      shiny::updateSelectInput(
        session=session,
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
      shiny::updateSelectInput(
        session=session,
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
  massReactants = getMassList(reactants, excludeList = dummySpecies)
  
  products = getSpecies(neutralsDB()$PRODUCTS[iReac])
  mask[['PRODUCTS']] = products
  massProducts = getMassList(products, excludeList = dummySpecies)
  
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
  list(
    br(),
    fluidRow(
      column(
        2,
        shiny::textInput(
          "neutralsReacReactants",
          "Reactants",
          value = paste0(mask[['REACTANTS']],collapse = ' + ')
        )
      ),
      column(
        2,
        shiny::textInput(
          "neutralsReacProducts",
          "Products",
          value = paste0(mask[['PRODUCTS']],collapse = ' + ')
        )
      ),
      column(
        2,
        shiny::selectInput(
          "neutralsReacTYPE",
          "Reaction type",
          neutralsReacTypes,
          selected = mask[['TYPE']]
        )
      )
    ),
    fluidRow(
      column(
        1,
        strong('k0')
      ),
      column(
        2,
        shiny::textInput("neutralsReacA1", "A1", value = mask[['A1']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacB1", "B1", value = mask[['B1']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacC1", "C1", value = mask[['C1']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacF1", "F1", value = mask[['F1']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacG1", "G1", value = mask[['G1']])
      )
    ),
    fluidRow(
      column(
        1,
        strong('kInf')
      ),
      column(
        2,
        shiny::textInput("neutralsReacA2", "A2", value = mask[['A2']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacB2", "B2", value = mask[['B2']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacC2", "C2", value = mask[['C2']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacF2", "F2", value = mask[['F2']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacG2", "G2", value = mask[['G2']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacFC", "Fc", value = mask[['FC']])
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
        1,
        shiny::textInput("neutralsReacB3", "B3", value = mask[['B3']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacC3", "C3", value = mask[['C3']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacF3", "F3", value = mask[['F3']])
      ),
      column(
        1,
        shiny::textInput("neutralsReacG3", "G3", value = mask[['G3']])
      )
    ),
    fluidRow(
      column(
        5,
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
          "neutralsReacREFS",
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

    ratePars = rep(NA,length(neutralsRateParKwdList))
    names(ratePars) = neutralsRateParKwdList
    for (kwd in neutralsRateParKwdList) 
      ratePars[kwd] = as.numeric(input[[paste0('neutralsReac',kwd)]])

    id = neutralsReacID()

    line = c(
      id,
      reactants,
      products,
      trimws(input$neutralsReacTYPE),
      ratePars,
      gsub(';',',',input$neutralsReacREFS),
      input$neutralsReacRQ,
      paste0(Sys.time())
    )
    line = paste0(line,collapse = ";")

    
    data = neutralsDB()
    tag  = makeTag(
      c(id,
        input$neutralsReacReactants,
        input$neutralsReacProducts)
    )
    iReac = which(neutralsDB()$TAG == tag)
    
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
observeEvent(
  input$neutralsSimulateBtn,
  {
    req(neutralsRateMask())
    
    nMC = as.numeric(input$neutralsSimulateSize)
    
    # Sanity checks !!!!!
    # TBD...
    # * Mass consistency
    
    # Rate parameters
    sampleRateParams = matrix(
      NA, 
      nrow = nMC,
      ncol = length(neutralsRateParKwdList))
    colnames(sampleRateParams) = neutralsRateParKwdList
    rateParDistStrings = rep(NA,length(neutralsRateParKwdList))
    names(rateParDistStrings) = neutralsRateParKwdList
    for (kwd in neutralsRateParKwdList) {
      # stringDist = neutralsRateMask()[[kwd]]
      stringDist = input[[paste0('neutralsReac',kwd)]] # Enable user mod
      rateParDistStrings[kwd]=stringDist
      sampleDist = sampleDistString(stringDist, nMC)
      sampleRateParams[1:nMC, kwd] = sampleDist
    }
    
   
    
    neutralsSimulSamples(
      list(
        sampleSize       = nMC,
        sampleRateParams = sampleRateParams,
        rateParDistStrings = rateParDistStrings
      )
    )
    
  }
)

# Plot rates ####
output$plotRate = renderPlot({
  req(input$reacNbPlot)
  iReac = as.numeric(input$reacNbPlot)
  
  req(reacScheme())
  reac0   = reacScheme()$reactants[[iReac]]
  prod0   = reacScheme()$products[[iReac]]
  params0 = reacScheme()$params[[iReac]]
  typ0    = reacScheme()$type[[iReac]]
  note    = reacScheme()$notes[[iReac]]
  
  tag  = paste0(
    'Reac. ',iReac, ': ',
    paste0(reac0, collapse = ' + '),' --> ',
    paste0(prod0, collapse = ' + ')
  )
  legText = paste0(tag, '\n', 'Rate law: ', typ0)
  legText = switch(
    typ0,
    kooij    = paste0(
      legText,
      '\n',
      'Parameters: ',
      paste0(params0[1:5],
             collapse = ' / '),
      '\n'
    ),
    assocMD  = paste0(
      legText,
      '\n',
      'Parameters Fc : ', params0[16],'\n',
      'Parameters k0 : ',
      paste0(params0[1:5], collapse = ' / '),
      '\n',
      'Parameters kInf : ',
      paste0(params0[6:10], collapse = ' / '),
      '\n',
      'Parameters kr : ',
      paste0(params0[11:15], collapse = ' / '),
      '\n'
    ),
    assocVV  = paste0(
      legText,
      '\n',
      'Parameters Fc : ', params0[16],'\n',
      'Parameters kInf : ',
      paste0(params0[1:5], collapse = ' / '),
      '\n',
      'Parameters k0 : ',
      paste0(params0[6:10], collapse = ' / '),
      '\n',
      'Parameters kR : ',
      paste0(params0[11:15], collapse = ' / '),
      '\n'
    ),
    legText
  )
  legText = paste0(legText,note)
  
  # Samples directory
  neutralsSampleDir = paste0(neutralsPublic,'_',neutralsVersion())
  req(dir.exists(neutralsSampleDir))
  samplesList = list.files(
    path = neutralsSampleDir,
    pattern = 'run_',
    full.names = TRUE
  )
  
  # Generate curves
  T0 = as.numeric(input$T0Plot)
  tRange = seq(input$tempRangePlot[1],input$tempRangePlot[2],5)
  M0  = 10^as.numeric(input$M0Plot)
  mRange = 10^seq(input$densRangePlot[1],input$densRangePlot[2],0.5)
  
  irun = 0
  krateT = matrix(NA,ncol = length(samplesList),nrow= length(tRange))
  krateM = matrix(NA,ncol = length(samplesList),nrow= length(mRange))
  for (file in samplesList) {
    irun = irun + 1
    
    # Get params and generate tags
    scheme  = read.csv(
      file = file,
      skip = iReac-1,
      nrows = 1,
      header = FALSE,
      sep = ';'
    )
    scheme  = gsub(" ", "", scheme)
    params = scheme[8:23]
    pars   = as.numeric(params)
    type   = scheme[ 24]
    
    krateT[,irun] = switch(
      type,
      kooij    = kooij(pars,tRange,M0),
      assocMD  = k3body(pars,tRange,M0),
      assocVV  = kEq18VV(pars,tRange,M0),
      rep(0,length(tRange))
    )
    
    # fixed T, M varies
    krateM[,irun] = switch(
      type,
      kooij    = kooij(pars,T0,mRange),
      assocMD  = k3body(pars,T0,mRange),
      assocVV  = kEq18VV(pars,T0,mRange),
      rep(0,length(mRange))
    )
    
    # if (sum(k <= 0) != 0)
    #   alerts = c(alerts, paste0('Null RC: ', tag[[i]], '\n'))
  }
  
  # Plot 
  
  trBlue = col2tr('blue', 60)
  
  par(
    mfrow = c(1, 2),
    mar = c(3,3,9,2),
    mgp = gPars$mgp,
    tcl = gPars$tcl,
    lwd = gPars$lwd,
    pty = 's',
    cex = 1.5
  )
  
  tempRange = tRange
  if (diff(range(krateT)) != 0) {
    matplot(
      tempRange,
      krateT,
      type = 'l',
      lty = 1,
      col = gPars$cols_tr2[6],
      lwd = 1.5*gPars$lwd,
      log = 'y',
      xlab = 'T [K]',
      ylab = 'Rate constant [cm^3.s^-1]',
      main = ''
    )
    grid()
    lines(tempRange, krateT[, 1], col = gPars$cols[2],lwd = 1.5*gPars$lwd)
    legend('top',title = paste0('M = ',M0,'cm^-3'), legend = NA, bty='n')
    mtext(
      legText,
      side = 3,
      cex = 1.5,
      adj = 0,
      line = 8,
      padj = 1,
      col = gPars$cols[1]
    )
    box()
  }
  
  # P-dep
  if (diff(range(krateM)) != 0) {
    matplot(
      mRange,
      krateM,
      type = 'l',
      lty = 1,
      col = gPars$cols_tr2[6],
      lwd = 1.5*gPars$lwd,
      log = 'xy',
      xlab = 'M [cm^-3]',
      ylab = 'Rate constant [cm^3.s^-1]',
      main = ''
    )
    lines(mRange, krateM[, 1], col = gPars$cols[2], lwd = 1.5*gPars$lwd)
    grid(col = 'darkgray')
    legend('top',title = paste0('T = ',T0,' K'), legend = NA, bty='n')
    box(lwd = 4)
  }  
  
},
height = plotHeight)

# Biblio ####
output$ionsBiblio = shiny::renderUI({
  req(ionsRateMask())
  req(ionsBRMask())
  
  keys = c()
  bibKwd = paste0('REF_',c(ionsRateParKwdList,'BR'))
  for (elt in bibKwd) {
    k = ionsRateMask()[[elt]]
    if(k != "NA" & !is.na(k) & k!="")
      keys = c(keys,unlist(str_split(k,';')))
  }
  keys = sort(unique(keys))
  
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
