photoReacs         = shiny::reactiveVal()
photoXSMask        = shiny::reactiveVal()
photoBRMask        = shiny::reactiveVal()
photoSimulSamples  = shiny::reactiveVal()
photoReacsFiltered = shiny::reactiveVal()
photoReacID        = shiny::reactiveVal()

# Manage reacs list ####
observeEvent(
  input$photoReacSelInit,
  {
    shiny::updateTextInput(inputId = "photoReacSel", value = "")
  })

observe({
  req(photoDB())
  reacs = photoDB()$REACTANTS
  tag   = reacs[ photoDB()$CHANNEL == 1 ] 
  
  if(input$photoReacSel != ""){
    sel = grepl(input$photoReacSel,reacs) & photoDB()$CHANNEL == 1
    if(sum(sel) == 0) {
      id = shiny::showNotification(
        strong(
          paste0('Unknown reactant: ',input$photoReacSel)
        ),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
      tag = NULL
    } else {
      tag = reacs[sel]
    }
  }
  req(tag)

  nums = 1:length(tag)  
  names(nums) = tag
  photoReacsFiltered(nums) # Partial list used by other parts
  
  shiny::updateSelectInput(
    session  = session,
    "photoReaction",
    choices  = nums,
    selected = 1
  )    
  
})

observeEvent(
  input$photoMinus,
  {
    iReac = as.numeric(input$photoReaction)
    if(iReac > 1) {
      iReac = iReac - 1
      shiny::updateSelectInput(
        session = session,
        "photoReaction",
        selected = iReac
      )      
    }
  })

observeEvent(
  input$photoPlus,
  {
    req(photoReacsFiltered())
    
    reacs = photoReacsFiltered()
    iReac = as.numeric(input$photoReaction)
    if(iReac < length(reacs)) {
      iReac = iReac + 1
      shiny::updateSelectInput(
        session = session,
        "photoReaction",
        selected = iReac
      )      
    }
  })

# Parse data file ####
shiny::observe({
  req(photoDB())
  req(input$photoReaction != "0")
  
  reacs = photoDB()$REACTANTS
  tag   = reacs[ photoDB()$CHANNEL == 1 ] 
  
  # Entry in table
  isel  = as.numeric(input$photoReaction)
  if(input$photoReacSel != ""){
    tag = names(photoReacsFiltered())[isel]
  } else {
    tag = tag[isel]
  }
  req(tag)
  
  # Absolute index of selected reac
  iReac = which(photoDB()$REACTANTS == tag & photoDB()$CHANNEL == 1)
  
  photoReacID(photoDB()$ID[iReac])
  
  # (Re-)Init samples for plots
  photoSimulSamples(NULL)
  
  # Format data for masks
  mask = list()
  reactants = getSpecies(photoDB()$REACTANTS[iReac])
  mask[['REACTANTS']] = reactants

  products = getSpecies(photoDB()$PRODUCTS[iReac])
  mask[['PRODUCTS']] = products

  for (kwd in photoKwdList)
    mask[[kwd]] = photoDB()[[kwd]][iReac]
  
  mask[['REFS']] = photoDB()$REFS[iReac]
  mask[['COMMENTS']] = photoDB()$COMMENTS[iReac]
  mask[['TIMESTAMP']] = photoDB()$TIMESTAMP[iReac]

  photoXSMask(mask)
  
  # # BRs
  # photoBRMask(
  #   list(
  #     nBR      = ionsDB()$NBR[iReac],
  #     StringBR = ionsDB()$STRINGBR[iReac],
  #     REF_BR   = mask[['REF_BR']]
  #   )
  # )
})

output$photoXSMask = shiny::renderUI({
  req(photoXSMask())
  
  mask = photoXSMask()
  
  uncF = mask[['XS_F']]
  if(is.null(uncF) | is.na(uncF) | uncF == "")
    uncF = photoDefaultuF
  
  list(
    br(),
    fixedRow(
      column(
        12,
        shiny::textInput(
          "photoReactants",
          "Reactants:",
          value = paste0(mask[['REACTANTS']],collapse = ' + ')
        ),
        fluidRow(
          column(
            6,
            shiny::selectInput(
              "photoXSSource",
              "Data source:",
              photoXSSources,
              selected = mask[['XS_SOURCE']]
            )
          ),
          column(
            6,
            conditionalPanel(
              condition = "input.photoXSSource != 'Leiden'",
              shiny::textInput(
                "photoXS_F",
                "Unc. Factor (>1) :",
                value = uncF
              )
            )
          )
        ),
        shiny::textInput(
          "photoReacREF",
          "Refs",
          value = mask[['REFS']]
        ),
        shiny::textAreaInput(
          "photoRQ",
          "Comments",
          value = mask[['COMMENTS']],
          height = '100px'
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
outputOptions(output, "photoXSMask", suspendWhenHidden = FALSE)

# output$ionsBRMask = shiny::renderUI({
#   req(ionsBRMask())
#   
#   list(
#     br(),
#     shiny::textAreaInput(
#       "ionsStringDist",
#       "StringBR",
#       width  = '400px',
#       height = '250px',
#       value  = formatStringBR(ionsBRMask()[['StringBR']])  
#     ),
#     fluidRow(
#       column(
#         6,
#         shiny::textInput(
#           "ionsReacREF_BR",
#           "Refs BR",
#           value = ionsBRMask()[['REF_BR']]
#         )
#       ),
#       column(
#         6,
#         shiny::textInput(
#           "ionsReacNBR",
#           "Nb Channels",
#           value = ionsBRMask()[['nBR']]
#         )
#       )
#     )
#   )
# })
# outputOptions(output, "ionsBRMask", suspendWhenHidden = FALSE)

# Save ####
# observeEvent(
#   input$ionsParseSave,
#   {
#     req(ionsDB())
#     
#     # Rate parameters
#     rateParDistStrings = rep(NA,length(ionsRateParKwdList))
#     names(rateParDistStrings) = ionsRateParKwdList
#     for (kwd in ionsRateParKwdList) 
#       rateParDistStrings[kwd] = input[[paste0('ionsReac',kwd)]]
#     
#     # Branching ratios
#     req(input$ionsStringDist)
#     stringDist = input$ionsStringDist # Enable user mod
#     stringDist = gsub('\n','',stringDist)
#     stringDist = gsub('\t','',stringDist)
#     tags = getTagsFromTaggedDist(stringDist)
#     
#     if(length(tags) != input$ionsReacNBR)
#       id = shiny::showNotification(
#         h4(paste0('>>> Pb. with channels number:',
#                   paste0(tags,collapse = ';'))),
#         closeButton = TRUE,
#         duration = NULL,
#         type = 'error'
#       )
#     
#     # Sanity checks !!!!!
#     # TBD...
#     # * Mass consistency
#     # * Correct distributions for Pars (Delta, Logu, Logn, Unif...)
#     # * Correct distributions for BRs (Diri, Diun, Dirg, Mlgn...)
#     
#     # Biblio
#     bibKwd = paste0('REF_',c(ionsRateParKwdList,'BR'))
#     refBib = rep(NA, length(bibKwd))
#     names(refBib) = bibKwd
#     for (kwd in bibKwd)
#       refBib[kwd] = input[[paste0('ionsReac',kwd)]]
#     
#     reactants = trimws(input$ionsReacReactants)
#     if(input$ionsParseComment)
#       reactants = paste0('#',reactants)
#     
#     id = ionsReacID()
#     line = c(
#       id,
#       reactants,
#       trimws(input$ionsReacTYPE),
#       rateParDistStrings,
#       length(tags),
#       trimws(stringDist),
#       refBib,
#       input$ionsReacRQ,
#       paste0(Sys.time())
#     )
#     line = matrix(line,nrow=1)
#     line = capture.output(
#       write.table(line,sep=';',row.names = FALSE, col.names = FALSE)
#     )
#     
#     # Update editor's content
#     data  = ionsDB()
#     tag = paste0(id,': ',input$ionsReacReactants)
#     iReac = which(data$TAG == tag)
#     
#     # Update editor's content
#     dataEditor = ionsEditDBText()
#     dl = length(dataEditor)
#     
#     if(length(iReac) != 0) {
#       # Tag exists: replace in situ
#       dataEditor[id+1] = line   # Header counts as first line...
#     } else {
#       # Tag does not exist: set new ID
#       line = sub(paste0('^',id),paste0(dl),line)
#       dataEditor = c(dataEditor,line)
#     }
#     
#     ionsEditDBText(dataEditor)
#     
#   }
# )


# Sampling ####
photoXSSimulate = reactive({
  req(photoXSMask())

  nMC  = as.numeric(input$photoSimulateSize)
  reso = as.numeric(input$photoXSReso)
  type = input$photoXSSource
  
  sp   = photoXSMask()[['REACTANTS']][1]
  
  # Get data
  source = file.path(photoSource,photoEditOrigVersion(),'Data', type)
  
  if(type == 'Leiden') {
    xsl  = getXShdf5( sp, source_dir = source)
    if(is.null(xsl)) {
      id = shiny::showNotification(
        strong(paste0('No data for:', sp,' in ',type)),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
      return(NULL)
    }
    wl   = xsl$wavelength
    xs   = xsl$photoabsorption
    uF   = xsl$uncF
    
  } else if (type == 'SWRI') {
    file = file.path(source,paste0(sp, '.dat'))
    if (!file.exists(file)) {
      id = shiny::showNotification(
        strong(paste0('No data for:', sp,' in ',type)),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
      return(NULL)
    }
    S = read.table(file = file,
                   header = TRUE,
                   check.names = FALSE)
    
    wl = S[, 1] / 10 # Convert A to nm
    xs = S[, 2]
    if(!is.null(input$photoXS_F))
      uF = as.numeric(input$photoXS_F)
    else
      uF = photoDefaultuF
    
  } else if (type == 'Hebrard') {
    file = file.path(source,paste0('se', sp, '.dat'))
    if (!file.exists(file)) {
      id = shiny::showNotification(
        strong(paste0('No data for:', sp,' in ',type)),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
      return(NULL)
    }
    S = read.table(file = file,
                   header = FALSE,
                   check.names = FALSE,
                   fill = TRUE)
    wl  = S[, 1]
    xs  = S[, 2]
    if(!is.null(input$photoXS_F))
      uF = as.numeric(input$photoXS_F)
    else
      uF = photoDefaultuF
    
  }
  
  # Interpolate on regular grid
  xsl  = downSample(wl, xs, reso = reso)
  wl   = xsl$wl
  xs   = xsl$xs
  
  # Remove tailing zeroes
  first = which(xs != 0)[1]
  last  = length(xs) - which(rev(xs) != 0)[1] + 1
  wl   = wl[first:last]
  xs   = xs[first:last]

  # Systematic perturbation
  sampleXS = matrix(NA, nrow = nMC + 1, ncol = length(wl))
  sampleXS[1, ] = xs
  for (iMC in 1:nMC) {
    rnd =  truncnorm::rtruncnorm(1,-3,3,0,1) # Avoid outliers
    # rnd = rlnorm(1, meanlog = 0, sdlog = log(uF))
    sampleXS[1 + iMC, ] = xs * exp( log(uF) * rnd )
  }
  
  return(list(
    sampleSize  = nMC,
    sampleXS    = sampleXS,
    sampleWl    = wl,
    sampleTitle = paste0(sp,' / ',type,' / ',
                         reso,' nm / uF = ',uF) 
  ))
  
})

# Plot XS ####
output$plotPhotoXSSample = renderPlot({
  req(photoXSMask())
  
  photoSimulSamples = photoXSSimulate()
  
  sampleXS    = photoSimulSamples$sampleXS
  sampleWl    = photoSimulSamples$sampleWl
  sampleTitle = photoSimulSamples$sampleTitle
  
  par(mar = c(4, 4, 2, 1),
      mgp = gPars$mgp,
      tcl = gPars$tcl,
      lwd = gPars$lwd,
      pty = 'm',
      cex = 1.5)
  matplot(
    sampleWl, t(sampleXS),
    type = 'l', 
    lwd  = 3,
    lty  = 1,
    xaxs = 'i',
    xlab = 'Wavelength [nm]',
    xlim = input$photoWLPlotRange,
    yaxs = 'i',
    ylim = c(0, 1.1*max(sampleXS[1,])),
    ylab = expression(paste('Cross-section [', cm ^ 2, ']')),
    col = gPars$cols_tr[5],
    main = sampleTitle
  )
  grid()
  lines(sampleWl, sampleXS[1,], lwd = 2, col= gPars$cols[2])
  box()
  
},
height = plotHeight, width = plotWidth)
# Plot BRs ####
## Sample ####
output$plotIonsBRSample = renderPlot({
  # req(ionsSimulSamples())
  ionsSimulSamples = ionsSimulate()
  
  req(!is.null(ionsSimulSamples$mytree))
  
  sampleSize       = ionsSimulSamples$sampleSize
  sampleBR         = ionsSimulSamples$sampleBR
  tags             = ionsSimulSamples$tags
  
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
