ionsDB            = shiny::reactiveVal()
ionsReacs         = shiny::reactiveVal()
ionsRateMask      = shiny::reactiveVal()
ionsBRMask        = shiny::reactiveVal()
ionsSimulSamples  = shiny::reactiveVal()
ionsReacsFiltered = shiny::reactiveVal()
ionsReacID        = shiny::reactiveVal()

# Manage reacs list ####
observeEvent(
  input$ionsReacSelInit,
  {
    shiny::updateTextInput(inputId = "ionsReacSel", value = "")
    shiny::updateRadioButtons(inputId = "ionsReacSelKind", selected = "Both")
  })

observe({
  req(ionsDB())
  reacs = ionsDB()$REACTANTS
  prods = ionsDB()$STRINGBR
  tag   = ionsDB()$TAG       # "tags" is reserved by shiny...

  if(input$ionsReacSel != ""){
    if(input$ionsReacSelKind == "Reactant") {
      sel = sapply(
        reacs,
        FUN = function(x) 
          input$ionsReacSel %in% getSpecies(x)
      )
    } else if(input$ionsReacSelKind == "Product"){
      sel = sapply(
        prods,
        FUN = function(x) 
          input$ionsReacSel %in% getSpecies(x)
      )
    } else if(input$ionsReacSelKind == "Both"){
      sel = apply(
        cbind(reacs,prods),
        1,
        FUN = function(x) 
          input$ionsReacSel %in% getSpecies(paste0(x[1],' + ',x[2]))
      ) 
    }
    if(sum(sel) == 0) {
      id = shiny::showNotification(
        strong(
          paste0(
            'Unknown species: ',input$ionsReacSel,
            ' in context: ',input$ionsReacSelKind)
        ),
        closeButton = TRUE,
        duration = NULL,
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
  ionsReacsFiltered(nums) # Partial list used by other parts
  
  shiny::updateSelectInput(
    session  = session,
    "ionsReaction",
    choices  = nums,
    selected = 1
  )    
  
})

observeEvent(
  input$ionsMinus,
  {
    iReac = as.numeric(input$ionsReaction)
    if(iReac > 1) {
      iReac = iReac - 1
      shiny::updateSelectInput(
        session = session,
        "ionsReaction",
        selected = iReac
      )      
    }
  })

observeEvent(
  input$ionsPlus,
  {
    req(ionsReacsFiltered())
    
    reacs = ionsReacsFiltered()
    iReac = as.numeric(input$ionsReaction)
    if(iReac < length(reacs)) {
      iReac = iReac + 1
      shiny::updateSelectInput(
        session = session,
        "ionsReaction",
        selected = iReac
      )      
    }
  })

# Parse data file ####
shiny::observe({
  req(ionsDB())
  req(input$ionsReaction != "0")
  
  # Entry in table
  isel  = as.numeric(input$ionsReaction)
  if(input$ionsReacSel != ""){
    tag = names(ionsReacsFiltered())[isel]
  } else {
    tag = ionsDB()$TAG[isel]
  }
  req(tag)
  
  # Absolute index of selected reac
  iReac = which(ionsDB()$TAG == tag)
  
  ionsReacID(ionsDB()$ID[iReac])
  
  # (Re-)Init samples for plots
  ionsSimulSamples(NULL)
  
  # Format data for masks
  mask = list()
  reactants = getSpecies(ionsDB()$REACTANTS[iReac])
  mask[['REACTANTS']] = reactants
  # massReactants = getMassList(reactants, excludeList = dummySpecies)
  
  mask[['TYPE']] = ionsDB()$TYPE[iReac]
  for (kwd in ionsRateParKwdList)
    mask[[kwd]] = ionsDB()[[kwd]][iReac]
  bibKwd = paste0('REF_',c(ionsRateParKwdList,'BR'))
  for (kwd in bibKwd)
    mask[[kwd]] = ionsDB()[[kwd]][iReac]
  mask[['RQ']] = ionsDB()$COMMENTS[iReac]
  mask[['TIMESTAMP']] = ionsDB()$TIMESTAMP[iReac]
  ionsRateMask(mask)
  
  # BRs
  ionsBRMask(
    list(
      nBR      = ionsDB()$NBR[iReac],
      StringBR = ionsDB()$STRINGBR[iReac],
      REF_BR   = mask[['REF_BR']]
    )
  )
  
  # Check mass balance
  stringDist = ionsDB()$STRINGBR[iReac]
  stringDist = gsub('\n','',stringDist)
  stringDist = gsub('\t','',stringDist)
  tags = getTagsFromTaggedDist(stringDist)
  for(i in 1:length(tags)) {
    products = trimws(strsplit(tags[i],' \\+ ')[[1]])
    msg = checkBalance(reactants,products)
    if(!is.null(msg))
      id = shiny::showNotification(
        strong(paste0('Reac ',iReac,', Channel ',i,': ',msg)),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
    req(is.null(msg))
  }
})

output$ionsRateMask = shiny::renderUI({
  req(ionsRateMask())
  
  mask = ionsRateMask()
  list(
    br(),
    fluidRow(
      column(
        6,
        shiny::textInput(
          "ionsReacReactants",
          "Reactants:",
          value = paste0(mask[['REACTANTS']],collapse = ' + ')
        ),
        shiny::textInput(
          "ionsReacALPHA",
          "ALPHA dist.:",
          value = mask[['ALPHA']]
        ),
        shiny::textInput(
          "ionsReacBETA",
          "BETA dist.",
          value = mask[['BETA']]
        ),
        shiny::textInput(
          "ionsReacGAMMA",
          "GAMMA dist.",
          value = mask[['GAMMA']]
        ),
        shiny::textAreaInput(
          "ionsReacRQ",
          "Comments",
          value = mask[['RQ']],
          height = '200px'
        )
      ),
      column(
        6,
        shiny::selectInput(
          "ionsReacTYPE",
          "Reaction type:",
          ionsReacTypes,
          selected = mask[['TYPE']]
        ),
        shiny::textInput(
          "ionsReacREF_ALPHA",
          "Refs ALPHA",
          value = mask[['REF_ALPHA']]
        ),
        shiny::textInput(
          "ionsReacREF_BETA",
          "Refs BETA",
          value = mask[['REF_BETA']]
        ),
        shiny::textInput(
          "ionsReacREF_GAMMA",
          "Refs GAMMA",
          value = mask[['REF_GAMMA']]
        ),
        shiny::textInput(
          "ionsReacTIMESTAMP",
          "Time Stamp",
          value = mask[['TIMESTAMP']]
        )
      )
    )
  )
})
outputOptions(output, "ionsRateMask", suspendWhenHidden = FALSE)

output$ionsBRMask = shiny::renderUI({
  req(ionsBRMask())
  
  list(
    br(),
    shiny::textAreaInput(
      "ionsStringDist",
      "StringBR",
      width  = '400px',
      height = '250px',
      value  = formatStringBR(ionsBRMask()[['StringBR']])  
    ),
    fluidRow(
      column(
        6,
        shiny::textInput(
          "ionsReacREF_BR",
          "Refs BR",
          value = ionsBRMask()[['REF_BR']]
        )
      ),
      column(
        6,
        shiny::textInput(
          "ionsReacNBR",
          "Nb Channels",
          value = ionsBRMask()[['nBR']]
        )
      )
    )
  )
})
outputOptions(output, "ionsBRMask", suspendWhenHidden = FALSE)

# Save ####
observeEvent(
  input$ionsParseSave,
  {
    req(ionsDB())
    
    # Rate parameters
    rateParDistStrings = rep(NA,length(ionsRateParKwdList))
    names(rateParDistStrings) = ionsRateParKwdList
    for (kwd in ionsRateParKwdList) 
      rateParDistStrings[kwd] = input[[paste0('ionsReac',kwd)]]
    
    # Branching ratios
    req(input$ionsStringDist)
    stringDist = input$ionsStringDist # Enable user mod
    stringDist = gsub('\n','',stringDist)
    stringDist = gsub('\t','',stringDist)
    tags = getTagsFromTaggedDist(stringDist)
    
    if(length(tags) != input$ionsReacNBR)
      id = shiny::showNotification(
        h4(paste0('>>> Pb. with channels number:',
                  paste0(tags,collapse = ';'))),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
    
    # Sanity checks !!!!!
    # TBD...
    # * Correct distributions for Pars (Delta, Logu, Logn, Unif...)
    # * Correct distributions for BRs (Diri, Diun, Dirg, Mlgn...)
    
    # Biblio
    bibKwd = paste0('REF_',c(ionsRateParKwdList,'BR'))
    refBib = rep(NA, length(bibKwd))
    names(refBib) = bibKwd
    for (kwd in bibKwd)
      refBib[kwd] = input[[paste0('ionsReac',kwd)]]
    
    reactants = trimws(input$ionsReacReactants)
    if(input$ionsParseComment)
      reactants = paste0('#',reactants)
    
    id = ionsReacID()
    
    line = c(
      reactants,
      trimws(input$ionsReacTYPE),
      rateParDistStrings,
      length(tags),
      trimws(stringDist),
      refBib,
      input$ionsReacRQ,
      paste0(Sys.time()),
      id
    )
    line = matrix(line,nrow=1)
    line = capture.output(
      write.table(line,sep=';',row.names = FALSE, col.names = FALSE)
    )
    
    # Update editor's content
    data  = ionsDB()
    tag = paste0(id,': ',input$ionsReacReactants)
    iReac = which(data$TAG == tag)
    
    # Update editor's content
    dataEditor = ionsEditDBText()
    dl = length(dataEditor)
    
    if(length(iReac) != 0) {
      # Tag exists: replace in situ
      dataEditor[id+1] = line   # Header counts as first line...
    } else {
      # Tag does not exist: set new ID
      line = sub(paste0('^',id),paste0(dl),line)
      dataEditor = c(dataEditor,line)
    }
    
    ionsEditDBText(dataEditor)
    
  }
)


# Sampling ####
ionsSimulate = reactive({
  req(ionsRateMask())
  req(ionsBRMask())
  
  nMC = as.numeric(input$ionsSimulateSize)
  
  # Sanity checks !!!!!
  # TBD...
  # * Correct distributions for Pars (Delta, Logu, Logn, Unif...)
  # * Correct distributions for BRs (Diri, Diun, Dirg, Mlgn...)
  
  # Rate parameters
  sampleRateParams = matrix(
    NA, 
    nrow = nMC,
    ncol = length(ionsRateParKwdList))
  colnames(sampleRateParams) = ionsRateParKwdList
  rateParDistStrings = rep(NA,length(ionsRateParKwdList))
  names(rateParDistStrings) = ionsRateParKwdList
  for (kwd in ionsRateParKwdList) {
    # stringDist = ionsRateMask()[[kwd]]
    stringDist = input[[paste0('ionsReac',kwd)]] # Enable user mod
    rateParDistStrings[kwd]=stringDist
    sampleDist = sampleDistString(stringDist, nMC)
    sampleRateParams[1:nMC, kwd] = sampleDist
  }
  
  # Branching ratios
  # stringDist = ionsBRMask()$StringBR
  req(input$ionsStringDist)
  stringDist = input$ionsStringDist # Enable user choice
  stringDist = gsub('\n','',stringDist)
  stringDist = gsub('\t','',stringDist)
  tags = getTagsFromTaggedDist(stringDist)
  
  if(length(tags) != input$ionsReacNBR)
    id = shiny::showNotification(
      h4(paste0('>>> Pb. with tags:',
                paste0(tags,collapse = ';'))),
      closeButton = TRUE,
      duration = NULL,
      type = 'error'
    )
  
  if(length(tags) > 1) {
    
    # Generate BR sample #
    stringBR = getDistFromTaggedDist(stringDist)
    sampleBR = nds(nMC,stringBR)
    
    # Build Newick string for tree plotting 
    newickBR = getNewickFromTaggedDist(stringDist)
    mytree <- ape::read.tree(text = newickBR)
    
    nodeTags = getNodesFromTaggedDist(stringDist)
    edgeTags = NULL # TBD
    treeDepth= max(ape::node.depth(mytree))
    
  } else {
    # Single pathway with BR=1
    sampleBR = matrix(1,ncol=1,nrow=nMC)
    # tags     = stringDist
    nodeTags = NULL
    edgeTags = NULL
    mytree   = NULL
    treeDepth= NULL
    
  }
  
  return(
    list(
      sampleSize       = nMC,
      sampleRateParams = sampleRateParams,
      rateParDistStrings = rateParDistStrings,
      sampleBR         = sampleBR,
      nodeTags         = nodeTags,
      edgeTags         = edgeTags,
      mytree           = mytree,
      treeDepth        = treeDepth,
      tags             = tags
    )
  )
  
}
)

# Plot rates ####
output$plotIonsParsSample = renderPlot({
  # req(ionsSimulSamples())
  req(ionsRateMask())
  
  mask = isolate(ionsRateMask())
  reacType  = tolower(mask$TYPE)
  reactants = mask$REACTANTS
  
  ionsSimulSamples = ionsSimulate()
  
  sampleSize         = ionsSimulSamples$sampleSize
  sampleRateParams   = ionsSimulSamples$sampleRateParams
  rateParDistStrings = ionsSimulSamples$rateParDistStrings
  
  meanPars = rep(NA, ncol(sampleRateParams))
  names(meanPars) = colnames(sampleRateParams)
  sigPars = rep(NA, ncol(sampleRateParams))
  names(sigPars) = colnames(sampleRateParams)
  for (kwd in ionsRateParKwdList) {
    sample = sampleRateParams[, kwd]
    if(substr(rateParDistStrings[kwd],1,5)!='Delta') {
      meanPars[kwd] = exp(mean(log(sample)))
      sigPars[kwd]  = exp(sd(log(sample)))
    } else {
      meanPars[kwd] = mean(sample)
      sigPars[kwd]  = 1
    }
  }
  
  split.screen(c(2, 1))
  split.screen(c(1,3),screen = 1)
  split.screen(c(1,2),screen = 2)
  iscreen=2
  for (kwd in ionsRateParKwdList) {
    iscreen = iscreen + 1
    screen(iscreen)
    if (!is.finite(sigPars[kwd]) | sigPars[kwd] == 1)
      next
    par(mar = c(4, 5, 3, 1),
        mgp = gPars$mgp,
        tcl = gPars$tcl,
        lwd = gPars$lwd,
        pty = 's',
        cex = 1)
    d = density(sampleRateParams[, kwd])
    plot(
      d,
      xlim = range(sampleRateParams[, kwd]),
      xaxs = 'i',
      xlab = kwd, #paste0(kwd, '~', rateParDistStrings[kwd]),
      yaxs = 'i',
      ylim = c(0,1.1*max(d$y)),
      col = gPars$cols_tr2[5],
      main = ''
    )
    polygon(d$x, d$y, col = gPars$cols_tr[5], border = NA)
    # grid()
    abline(v = meanPars[kwd],
           col = gPars$cols[2],
           lwd = 2)
    abline(
      v = ErrViewLib::vhd(sampleRateParams[, kwd],p = c(0.05,0.95)),
      col = gPars$cols[2],
      lty = 2,
      lwd = 1
    )
    box()
  }
  
  temp = seq(input$ionsTempRangePlot[1],
             input$ionsTempRangePlot[2],
             10)
  nt = length(temp)
  if (reacType == 'kooij' | reacType == 'dr') {
    if ('E' %in% reactants) {
      rateFun = function(t, pars)
        pars['ALPHA'] * (300 / t) ^ pars['BETA']
    } else {
      rateFun = function(t, pars)
        pars['ALPHA'] * (t / 300) ^ pars['BETA'] * exp(-pars['GAMMA'] / t)
    }
  } else {
    if (reacType == 'ionpol1') {
      # KIDA::ionpol1 type
      # BEWARE: notation change because here branching ratios are stored in BR
      rateFun = function(t, pars)
        pars['ALPHA'] * (0.62 + 0.4767 * pars['BETA'] * (300 / t) ^ 0.5)
    } else {
      # KIDA::ionpol2 type
      # BEWARE: notation change because here branching ratios are stored in BR
      rateFun = function(t, pars)
        pars['ALPHA'] * (1 + 0.0967 * pars['BETA'] * (300 / t) ^ 0.5 +
                           pars['BETA'] ^ 2 * 300 / (10.526 * t))
    }
  }
  
  np = min(sampleSize, 500) # nb of plotted samples
  Y = matrix(NA, nrow = np, ncol = nt)
  for (ip in 1:np)
    Y[ip, 1:nt] = rateFun(temp, sampleRateParams[ip, ])
  yMean = 10^apply(log10(Y),2,mean)
  yF    = 10^apply(log10(Y),2,sd)
  screen(6)
  par(mar = c(4, 5, 1, 1),
      mgp = gPars$mgp,
      tcl = gPars$tcl,
      lwd = gPars$lwd,
      pty = 'm',
      cex = 1)
  
  matplot(
    temp, t(Y),
    type = 'l', 
    lty = 1,
    col = gPars$cols_tr[5],
    lwd = 4,
    log = 'y',
    xlab = 'T / K',
    xaxs = 'i',
    ylab = 'Rate constant / cm^3.s^-1',
    ylim = c(min(Y) / 1.5, max(Y) * 1.5),
    main = ''
  )  
  grid()
  lines(temp,
        yMean, #rateFun(temp, meanPars),
        col = gPars$cols[2],
        lwd = 2)
  box()
  
  screen(7)
  par(mar = c(4, 5, 1, 1),
      mgp = gPars$mgp,
      tcl = gPars$tcl,
      lwd = gPars$lwd,
      pty = 'm',
      cex = 1)
  plot(
    temp, yF,
    type = 'l',
    lty = 1,
    col = gPars$cols[2],
    lwd = 3,
    xlab = 'T / K',
    xaxs = 'i',
    ylab = 'Uncert. factor',
    main = ''
  )  
  grid()
  box()
  close.screen(all = TRUE)
},
height = plotHeight, width = 1.2*plotWidth)
# Plot BRs ####
## Sample ####
output$plotIonsBRSample = renderPlot({
  req(ionsSimulate())
  req(ionsRateMask())
  
  ionsSimulSamples = ionsSimulate()
  tags = ionsSimulSamples$tags
  
  req(!is.null(ionsSimulSamples$mytree))
  
  sampleSize       = ionsSimulSamples$sampleSize
  sampleBR         = ionsSimulSamples$sampleBR
  # tags             = ionsSimulSamples$tags
  
  meanBR    = colMeans(sampleBR)
  meanBR    = meanBR / sum(meanBR)
  sigBR     = apply(sampleBR, 2, sd)
  
  np = min(sampleSize, 500) # nb of plotted samples
  nt = length(tags)
  if (nt >= 2) {
    m0= max(1,40-3*nt)
    par(
      mar = c(m0, 15, 4, 1), 
      cex.main = 1)
    matplot(
      t(sampleBR)[nt:1, 1:np],
      1:nt,
      type = 'l',
      lty = 1,
      col = gPars$cols_tr[5],
      lwd = 5,
      yaxt = 'n',
      ylab = '',
      yaxs = 'i',
      xlim = c(0, 1),
      xlab = '',
      xaxt = 'n',
      xaxs = 'i',
      main = 'Branching Ratios'
    )
    lines(meanBR, nt:1, lwd = 3, col = gPars$cols[2])
    lines(cumsum(meanBR), nt:1, lwd = 3, lty = 2, col = gPars$cols[2])
    text(
      x = seq(0, 1, by = 0.1),
      y = nt,
      labels = seq(0, 1, by = 0.1),
      xpd = TRUE,
      cex = 1,
      pos = 3
    )
    text(
      x = -0.02,
      y = nt:1,
      labels = tags,
      srt = 0,
      adj = 1,
      xpd = TRUE,
      cex = 1
    )
    # grid()
    abline(
      h = 1:length(tags),
      col = 'darkgray',
      lwd = 1,
      lty = 3
    )
    abline(
      v = seq(0.1,1,by=0.1),
      col = 'darkgray',
      lwd = 1,
      lty = 3
    )
    legend('topright', bty = 'n',
           legend = c('Sample', 'Mean', 'Cumul.'),
           lty = c(1,1,2), 
           lwd = 3, pch = NA,
           col = gPars$cols[c(5,2,2)])
    box()
  }
},
height = plotHeight, width = plotWidth)
## Tree ####
output$plotIonsBRTree = renderPlot({
  
  ionsSimulSamples = ionsSimulate()
  req(!is.null(ionsSimulSamples$mytree))
  
  sampleSize       = ionsSimulSamples$sampleSize
  sampleBR         = ionsSimulSamples$sampleBR
  nodeTags         = ionsSimulSamples$nodeTags
  edgeTags         = ionsSimulSamples$edgeTags
  mytree           = ionsSimulSamples$mytree
  tags             = ionsSimulSamples$tags
  depth            = ionsSimulSamples$treeDepth
  
  meanBR    = colMeans(sampleBR)
  meanBR    = meanBR / sum(meanBR)
  sigBR     = apply(sampleBR, 2, sd)
  
  # Branching ratios
  nt = length(tags)
  if (nt >= 2) {
    tagStat = tags
    for (ip in 1:nt) {
      tagStat[ip] = paste0(' ',tags[ip], ' (',
                           signif(meanBR[ip], 2), ' +/- ',
                           signif(sigBR[ip], 1), ')')
    }
    m0 = max(1,40-3*nt)
    m1 = max(1,40-3*depth)
    par(
      mar = c(m0, 1, 1, m1), 
      cex = 1)
    mytree$tip.label = tagStat
    mytree$edge.length = rep(1, dim(mytree$edge)[1])
    plot(
      mytree,
      type = 'clado',
      y.lim = c(length(tags), 1), # Reverse axis
      show.tip.label = TRUE,
      tip.color = gPars$cols[1],
      use.edge.length = TRUE,
      root.edge = TRUE,
      edge.width = 2,
      edge.color = gPars$cols[3],
      align.tip.label = TRUE,
      font = 1,
      main = ''
    )
    mynodelabels(nodeTags, bg = 'gold')
    # myedgelabels(edgeTags) # Not defined yet...
    
  }
},
height = plotHeight, width = 1.5*plotWidth)

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
