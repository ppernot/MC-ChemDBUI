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
        7,
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
observeEvent(
  input$photoParseSave,
  {
    req(photoDB())
    req(photoBRMask())
    
    nBR = photoBRMask()$nBR
    Lines = c()
    
    # ID;R1;R2;P1;P2;P3;P4;CHANNEL;XS_SOURCE;XS_F;BR_SOURCE;REFS;COMMENTS;TIMESTAMP
    reactants = getSpecies(trimws(input$photoReactants))
    R1 = reactants[1]
    R2 = reactants[2]
    if(input$ionsParseComment)
      R1 = paste0('#',R1)
    
    XS_SOURCE = input$photoXSSource
    XS_F      = input$photoXS_F
    
    id = photoReacID()
    
    P1 = P2 = P3 = P4 = ""
    for( i in 1:nBR) {
      products = getSpecies(trimws(input[[paste0('photoBRProds_',i)]]))
      P1 = products[1]
      if(length(products) >=2)
        P2 = products[2]
      if(length(products) >=3)
        P3 = products[3]
      if(length(products) ==4)
        P4 = products[4]
      
      line = c(
        id + i -1,
        R1, R2, P1, P2, P3, P4, i,
        input$photoXSSource, 
        input$photoXS_F, 
        input$photoBRSource, 
        input$photoReacREF, 
        input$photoRQ,
        paste0(Sys.time())
      )
      line = matrix(line,nrow=1)
      Lines[i] = capture.output(
        write.table(line,sep=';',row.names = FALSE, col.names = FALSE)
      )
    }
    
    # Update editor's content
    dataEditor = photoEditDBText()
    # dl = length(dataEditor)
    dataEditor[(id+1):(id+nBR)] = Lines
    
    # if(length(iReac) != 0) {
    #  # Tag exists: replace in situ
    #   dataEditor[id+1] = line   # Header counts as first line...
    # } else {
    #   # Tag does not exist: set new ID
    #   line = sub(paste0('^',id),paste0(dl),line)
    #   dataEditor = c(dataEditor,line)
    # }

    photoEditDBText(dataEditor)

  }
)


# Sampling ####
## XS ####
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
## BRs ####
photoBRSimulate = shiny::reactive({
  req(photoDB()) 
  req(photoXSMask())
  req(photoBRMask())
  
  nMC  = as.numeric(input$photoSimulateSize)
  reso = as.numeric(input$photoXSReso)
  type = photoBRMask()$source
  req(type)
  
  sort = FALSE # input$photoBRSort

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
    qy0 = matrix(0, ncol = np, nrow = length(wl1))
    for (i in 1:np) {
      br = downSample(wl, S[, i+i0]/xs, reso = reso)$xs
      qy0[ ,i] = br[first:last] 
    }
    wl = wl1
    
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
    qy0 = matrix(0, ncol = length(files), nrow = length(wl))
    
    for (i in seq_along(files)) {
      S = read.table(file = file.path(source,files[i]),
                     header = TRUE,
                     check.names = FALSE)
        wl  = S[, 1]
        xs  = S[, 2]
        xsl = downSample(wl, xs, reso = reso)
        wl  = xsl$wl
        xs  = xsl$xs
        qy0[ ,i] = xs
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
    
    qy = qy0 / rowSums(qy0)
    
    if(input$photoEditGP_Fit) {
      
      nw = length(wl)
      
      qySample = array(data = NA,
                       dim = c(nMC, nw, nBR))
      ## Apply threshold
      mask = qy < photoEps
      qy[mask] = photoEps
      
      # Transform to logit-space
      qy1   = log( qy/(1-qy) )
      noise = 0.2 * sqrt( qy*(1-qy) ) # Distrib of uncert peaks at 0.06
      unc1  = noise/( qy*(1-qy) )
      
      ##  GP
      gp = list()
      kernel = c('gauss', 'matern5_2', 'matern3_2', 'exp')[2]
      coef.var = input$photoGPCorVar
      coef.cov = input$photoGPCorLen
      sel = seq(1, nw, length.out = 30)
      for (i in 1:nBR) {
        gp[[i]] = DiceKriging::km(
          design   = data.frame(x=wl[sel]),
          response = data.frame(y=qy1[sel,i]),
          covtype  = kernel,
          coef.trend = mean(qy1[,i]),
          coef.cov = coef.cov,
          coef.var = coef.var,
          noise.var= unc1[sel,i]^2
        )
      }
      
      for (iMC in 1:nMC) {
        pred = matrix(nrow = nw, ncol = nBR)
        for (i in 1:nBR) {
          p = DiceKriging::simulate(
            gp[[i]],
            newdata = data.frame(x = wl),
            cond = TRUE)
          pred[, i] = p
        }
        ## back-transform
        brPred = 1 / (1 + exp(-pred))
        brPred[mask] = 0 # Restore true 0s
        qySample[iMC, , ] = brPred / rowSums(brPred)
      }

    } else {
      # Generate ordered Diri-based samples at each wavelength
      
      ionic = c()
      for (i in 1:nBR) {
        iReac = photoBRMask()$channels[i]
        ionic[i] = 'E' %in% getSpecies(photoDB()$PRODUCTS[iReac])
      }
      
      if (input$flatTree) {
        
        uqy =  matrix(0, ncol =ncol(qy), nrow = nrow(qy))
        if(input$useDirg) {
          # Uncertainties on BRs
          r = rep(photoRuBRN,nBR)
          if(sum(ionic) != 0)
            r[ionic] = photoRuBRI 
          for(i in 1:nrow(qy0))
            uqy[i,] = fuBr(qy0[i,],r*qy0[i,])
        }
          
        qySample = flatSample(
          qy, uqy,
          ru = photoRuBRN,
          nMC = nMC,
          eps = photoEps,
          useDirg = input$useDirg,
          newDiri = input$newDiri,
          newDirg = input$newDirg)
        
      } else {
        
        # Nested sampling
        qySample = hierSample(
          qy0, 
          ionic = ionic,
          ru = c(photoRuBRNI, photoRuBRN, photoRuBRI),
          nMC = nMC,
          eps = photoEps,
          useDirg = input$useDirg,
          newDiri = input$newDiri,
          newDirg = input$newDirg)
      }
      
      # Rearrange samples for better wavelength continuity
      if(input$photoBRArrange)
        qySample = arrangeSample(
          qySample, 
          useRanks = input$photoBRUseRanks
        )
      
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
  
  log = ''
  ylim = c(0, 1.1*max(sampleXS[1,]))
  if(input$photoXSLog){
    log = 'y'  
    ylim = NULL
  }

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
    log  = log,
    ylim = ylim,
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
  req(photoDB())
  req(photoBRMask())
  
  photoSimulSamples = photoBRSimulate()
  nMC         = photoSimulSamples$sampleSize
  sampleBR0   = photoSimulSamples$sampleBR0
  sampleBR    = photoSimulSamples$sampleBR
  sampleWl    = photoSimulSamples$sampleWl
  sampleTitle = photoSimulSamples$sampleTitle
  nBR         = photoBRMask()$nBR
  channels    = photoBRMask()$channels
  
 
  par(mar = c(4, 4, 2, 1),
      mgp = gPars$mgp,
      tcl = gPars$tcl,
      lwd = gPars$lwd,
      pty = 'm',
      cex = 1.5)
  
  cols    = rep(gPars$cols,2)
  cols_tr = rep(gPars$cols_tr,2)
  lty     = c(rep(1,length(gPars$cols)),rep(2,length(gPars$cols)))
  
  ylab = 'Branching ratios'
  ylim = c(-0.01, 1.4)
  
  if (input$photoEditBRDisplay == 0 | nBR == 1) {
    qy0 = sampleBR0
    qy  = sampleBR
    leg = photoDB()$PRODUCTS[channels]
    
    
  } else if (input$photoEditBRDisplay == 1) {
    ionic = c()
    for (i in 1:nBR)
      ionic[i] = 'E' %in% getSpecies(photoDB()$PRODUCTS[channels[i]])
    qy0 = matrix(data = 0, nrow = length(sampleWl), ncol = 2)
    qy0[, 1] = rowSums(sampleBR0[, !ionic, drop = FALSE])
    qy0[, 2] = rowSums(sampleBR0[,  ionic, drop = FALSE])
    qy = array(data = 0, dim  = c(nMC, length(sampleWl), 2))
    qy[, , 1] = rowSums(sampleBR[, , !ionic, drop = FALSE], dims = 2)
    qy[, , 2] = rowSums(sampleBR[, ,  ionic, drop = FALSE], dims = 2)
    leg = c("Neutrals", "Ions")
    
  } else if (input$photoEditBRDisplay == 2) {
    qy0 = rowSums(sampleBR0)
    qy = array(data = 0, dim  = c(nMC, length(sampleWl), 1))
    qy[, , 1] = rowSums(sampleBR, dims = 2)
    leg = c("Sum-to-one")
    
  } else if (input$photoEditBRDisplay == 3) {
    mu  = apply(sampleBR, c(2,3), mean, na.rm = TRUE)
    sig = apply(sampleBR, c(2,3), sd  , na.rm = TRUE)
    mu[mu==0] = 1
    qy0 = sig / mu
    leg = photoDB()$PRODUCTS[channels]
    ylab = 'Relative uncertainty'
    ylim = c(-0.01,1.4*max(qy0, na.rm = TRUE))
    
  } else if (input$photoEditBRDisplay == 4) {
    ionic = c()
    for (i in 1:nBR)
      ionic[i] = 'E' %in% getSpecies(photoDB()$PRODUCTS[channels[i]])
    
    mu  = apply(sampleBR, c(2,3), mean, na.rm = TRUE)
    sig = apply(sampleBR, c(2,3), sd  , na.rm = TRUE)
    mu[mu==0] = 1
    qy = sig / mu
    nzMean = function(X) {
      x = X[X != 0]
      if(length(x) == 0) 
        return(0)
      else
        return( mean(x) )
    }
    qy0 = matrix(data = 0, nrow = length(sampleWl), ncol = 3)
    qy0[,1] = apply(qy[, !ionic, drop = FALSE], 1, nzMean)
    qy0[,2] = apply(qy[,  ionic, drop = FALSE], 1, nzMean)
    qy0[,3] = apply(qy, 1, nzMean)
    leg = c("Neutrals", "Ions", "Global")
    ylab = 'Mean relative uncertainty'
    ylim = c(-0.01,1.4*max(qy0, na.rm = TRUE))
  }
  
  matplot(
    sampleWl, qy0,
    type = 'l', 
    lwd  = 3,
    xaxs = 'i',
    xlab = 'Wavelength [nm]',
    xlim = input$photoWLPlotRange,
    yaxs = 'i',
    ylim = ylim,
    ylab = ylab,
    col  = cols,
    lty  = lty,
    main = sampleTitle
  )
  grid()

  if (input$photoEditBRDisplay < 3) {
    for(iMC in 1:nMC) {
      matlines(
        sampleWl, qy[iMC,,],
        lwd  = 3,
        col  = cols_tr,
        lty  = lty
      )
    }
    matlines(
      sampleWl, qy0, 
      lwd = 4, 
      col = cols, 
      lty = lty
    )
  }
  
  legend(
    'top', ncol =2, box.col = "white",
    legend = leg,  cex = 0.75,
    col = cols,
    lty = lty,
    lwd = 4,
    pch = NULL
  )
  
  box()
  
},
height = plotHeight, width = plotWidth)


# Biblio ####
output$photoBiblio = shiny::renderUI({
  req(photoXSMask())
  req(input$photoReacREF)
  
  k = input$photoReacREF
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

