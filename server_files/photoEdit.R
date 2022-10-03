photoReacs         = shiny::reactiveVal()
photoXSMask        = shiny::reactiveVal()
photoBRMask        = shiny::reactiveVal()
photoReacsFiltered = shiny::reactiveVal()
photoReacID        = shiny::reactiveVal()

# Manage reacs list ####
observeEvent(
  input$photoReacSelInit,
  shiny::updateTextInput(inputId = "photoReacSel", value = "")
)

shiny::observe({
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

shiny::observeEvent(input$photoMinus,
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

shiny::observeEvent(input$photoPlus,
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
  
  # BRs
  chans  = which(photoDB()$REACTANTS == tag)
  nBR    = length(chans)
  source = photoDB()[['BR_SOURCE']][iReac] # Assumed identical for all channels
  
  photoBRMask(
    list(
      nBR      = nBR,
      channels = chans,
      source   = source
    )
  )
  
})

output$photoXSMaskUI = shiny::renderUI({
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
shiny::outputOptions(output, "photoXSMaskUI", suspendWhenHidden = FALSE)

output$photoBRMaskUI = shiny::renderUI({
  req(photoBRMask())
  req(photoDB())

  mask = photoBRMask()
  
  out = list()
  out[[1]] = br() 
  for( i in 1:mask$nBR) {
    iReac = mask$channels[i]
    out[[i+1]] = fluidRow(
      column(
        2,
        shiny::textInput(
          paste0("photoBRChan_",i),
          label = NULL,
          value = i
        )
      ),
      column(
        1,
        shiny::checkboxInput(
          paste0("photoBROut_",i),
          label = NULL,
          value = photoDB()[['OUTPUT']][iReac]
        )
      ),
      column(
        6,
        shiny::textInput(
          paste0("photoBRProds_",i), 
          label = NULL,
          value = photoDB()[['PRODUCTS']][iReac]
        )
      ),
      if (i==1)
        column(
          3,
          shiny::selectInput(
            "photoBRSource",
            label = NULL,
            choices = photoBRSources,
            selected = mask$source
          )
        )
    )
  }

  return(out)
})
shiny::outputOptions(output, "photoBRMaskUI", suspendWhenHidden = FALSE)

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
photoXSSimulate = shiny::reactive({
  req(photoXSMask())

  nMC  = as.numeric(input$photoSimulateSize)
  reso = as.numeric(input$photoXSReso)
  type = input$photoXSSource
  req(type)
  
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

photoBRSimulate = shiny::reactive({
  req(photoDB()) 
  req(photoXSMask())
  req(photoBRMask())
  
  nMC  = as.numeric(input$photoSimulateSize)
  reso = as.numeric(input$photoXSReso)
  type = photoBRMask()$source
  req(type)

  sp   = photoXSMask()[['REACTANTS']][1]
  nBR  = photoBRMask()$nBR
  
  # Get data
  source = file.path(photoSource,photoEditOrigVersion(),'Data', type)
  
  if (type == 'SWRI') {
    file = file.path(source,paste0(sp, '.dat'))
    if (!file.exists(file)) {
      id = shiny::showNotification(
        strong(paste0('No BR data for:', sp,' in ',type)),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
      return(NULL)
    }
    S = read.table(file = file,
                   header = TRUE,
                   check.names = FALSE)
    wl  = S[, 1] / 10 # Convert A to nm
    xs  = S[, 2]
    xsl = downSample(wl, xs, reso = reso)
    wl1  = xsl$wl
    xs1  = xsl$xs
    # Remove tailing zeroes
    first = which(xs1 != 0)[1]
    last  = length(xs1)-which(rev(xs1) !=0)[1] + 1
    wl1    = wl1[first:last]
    
    # Quantum yields
    np = ncol(S) - 2
    if (np == 0) {
      # Single channel: compute qy from total xs => qy=1
      np = 1
      i0 = 1
    } else {
      # Normal case: several channels
      i0 = 2
    }
    qy = matrix(0, ncol = np, nrow = length(wl1))
    for (i in 1:np) {
      br = downSample(wl, S[, i+i0]/xs, reso = reso)$xs
      qy[ ,i] = br[first:last] 
    }
    wl = wl1
    
  #   if(!is.null(input$photoXS_F))
  #     uF = as.numeric(input$photoXS_F)
  #   else
  #     uF = photoDefaultuF
  #   
  } else if (type == 'Plessis') {
    
    pattern = paste0('qy', sp)
    files = list.files(source, pattern)
    if (length(files) == 0) {
      id = shiny::showNotification(
        strong(paste0('No BR data for:', sp,' in ',type)),
        closeButton = TRUE,
        duration = NULL,
        type = 'error'
      )
      return(NULL)
    }

    S = read.table(file = file.path(source,files[1]),
                   header = TRUE,
                   check.names = FALSE)
    wl  = S[, 1]
    xs  = S[, 2]
    xsl = downSample(wl, xs, reso = reso)
    wl  = xsl$wl
    qy = matrix(0, ncol = length(files), nrow = length(wl))
    
    for (i in seq_along(files)) {
      S = read.table(file = file.path(source,files[i]),
                     header = TRUE,
                     check.names = FALSE)
        wl  = S[, 1]
        xs  = S[, 2]
        xsl = downSample(wl, xs, reso = reso)
        wl  = xsl$wl
        xs  = xsl$xs
        qy[ ,i] = xs
    }
    
  }

  # Sample
  
  if(nBR == 1) {
    # A single channel: no uncertainty in BR
    qySample = array(
      data = 1,
      dim = c(nMC, ncol(qy), nBR)
    )

  } else {

    qy = qy / rowSums(qy)

    # Generate ordered Diri-based samples at each wavelength
    
    ionic = c()
    for (i in 1:nBR) {
      iReac = photoBRMask()$channels[i]
      ionic[i] = 'E' %in% photoDB()$PRODUCTS[iReac]
    }
    
    if (sum(ionic) * sum(!ionic) == 0) {
      # Diri sampling
      qySample = diriSample(
        qy,
        ru = ifelse( sum(ionic) == 0, photoRuBRN, photoRuBRI),
        nMC,
        photoEps)
    } else {
      # Nested sampling
      qySample = hierSample(
        qy, ionic,
        ru = c(photoRuBRNI, photoRuBRN, photoRuBRI),
        nMC,
        photoEps)
    }
  }
  
  return(list(
    sampleSize  = nMC,
    sampleBR0   = qy,
    sampleBR    = qySample,
    sampleWl    = wl,
    sampleTitle = paste0(sp,' / ',type,' / ',reso) 
  ))
  
})


# Plot XS ####
output$plotPhotoXSSample = shiny::renderPlot({
  req(photoXSMask())
  req(input$photoXSSource)
  
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
output$plotPhotoBRSample = shiny::renderPlot({
  req(photoBRMask())
    photoSimulSamples = photoBRSimulate()

  nMC         = photoSimulSamples$sampleSize
  sampleBR0   = photoSimulSamples$sampleBR0
  sampleBR    = photoSimulSamples$sampleBR
  sampleWl    = photoSimulSamples$sampleWl
  sampleTitle = photoSimulSamples$sampleTitle
  nBR         = photoBRMask()$nBR
 
  par(mar = c(4, 4, 2, 1),
      mgp = gPars$mgp,
      tcl = gPars$tcl,
      lwd = gPars$lwd,
      pty = 'm',
      cex = 1.5)
  
  cols    = rep(gPars$cols,2)
  cols_tr = rep(gPars$cols_tr,2)
  lty     = c(rep(1,length(gPars$cols)),rep(2,length(gPars$cols)))
  
  matplot(
    sampleWl, sampleBR[1,,],
    type = 'l', 
    lwd  = 3,
    xaxs = 'i',
    xlab = 'Wavelength [nm]',
    xlim = input$photoWLPlotRange,
    yaxs = 'i',
    ylim = c(-0.01, 1.01),
    ylab = paste('Branching ratios'),
    col  = cols_tr,
    lty  = lty,
    main = sampleTitle
  )
  grid()

  for(iMC in 2:nMC) {
    matlines(
      sampleWl, sampleBR[iMC,,],
      lwd  = 3,
      col  = cols_tr,
      lty  = lty
    )
  }
  
  matlines(
    sampleWl, 
    sampleBR0, 
    lwd = 2, 
    col = cols, 
    lty = lty
  )
  
  legend(
    'right', bty = 'n',
    legend = 1:nBR,
    col = cols,
    lty = lty,
    pch = NULL
  )
  
  box()
  
},
height = plotHeight, width = plotWidth)


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
