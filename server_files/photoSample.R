photoSampleReacsFiltered = shiny::reactiveVal()

# Sample ####
observeEvent(
  input$photoSampleBtn,
  {
    req(photoDB())
    
    nMC = as.numeric(input$photoSampleSize) 
    
    sortBR = FALSE #input$photoBRSampleSort
    
    reacs = photoDB()$REACTANTS
    
    resos = input$photoSampleReso
    if(resos == "All")
      resos = photoXSResolutions
    else
      resos = as.numeric(resos)
    
    photoSampleDir = paste0(photoPublic, '_', photoEditOrigVersion())
    if (!dir.exists(photoSampleDir))
      dir.create(photoSampleDir)
    
    for (reso in resos){
      target = file.path(photoSampleDir, paste0(reso,'nm'))
      if(!dir.exists(target)) {
        dir.create(target)
        
      } else {
        # Clean
        if(!input$photoSampleCheck) {
          fileList = list.files(path = target, full.names = TRUE)
          file.remove(fileList)
        }
      }
    }
    
    photoSchemeFile = file.path(photoSampleDir,'PhotoScheme.dat')
    schemeTab = data.frame(
      R1 = photoDB()$R1,
      R2 = photoDB()$R2,
      P1 = photoDB()$P1,
      P2 = photoDB()$P2,
      P3 = photoDB()$P3,
      P4 = photoDB()$P4,
      CHANNEL = photoDB()$CHANNEL
    )
    Lines <- with(
      schemeTab,
      sprintf(
        "%-11s%-11s%-11s%-11s%-11s%-11s%-11.2i", 
        R1, R2, P1, P2, P3, P4, CHANNEL
      )
    )
    writeLines(Lines,photoSchemeFile)
    
    photoSourceFile = file.path(photoSampleDir,'provenance.txt')
    if(file.exists(photoSourceFile))
      file.remove(photoSourceFile)
    
    # Scale for progress / factor 2 is for XS and BRs
    len = 2 * sum(photoDB()$CHANNEL == 1) * length(resos) 
    shiny::withProgress(
      message = 'Sampling : ', 
      {
        for (reso in resos){
          
          target = file.path(photoSampleDir, paste0(reso,'nm'))

          for (iReac in seq_along(photoDB()$REACTANTS)) {
            

            # Loop on reactant species only
            if(photoDB()$CHANNEL[iReac] != 1)
              next
    
            ## Get XS ####                    
            
            # Get descriptors
            reactants = getSpecies(photoDB()$REACTANTS[iReac])
            sp = reactants[1]
            type     = photoDB()$XS_SOURCE[iReac]
            uF0      = photoDB()$XS_F[iReac]
            refs     = photoDB()$REFS[iReac]
            comments = photoDB()$COMMENTS[iReac]
            
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
              wl0 = xsl$wavelength
              xs0 = xsl$photoabsorption
              uF0 = xsl$uncF # Priority on DB parameter
              file = file.path(source, 'cross_sections', 
                               sp, paste0( sp, '.hdf5') )
              
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
              
              wl0 = S[, 1] / 10 # Convert A to nm
              xs0 = S[, 2]
              
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
              wl0 = S[, 1]
              xs0 = S[, 2]
            }
            
            if(!is.null(uF0))
              uF = as.numeric(uF0)
            else
              uF = photoDefaultuF
            
            sink(file = photoSourceFile, append = TRUE)
            cat(sp,'\n','  XS :',file,'; Unc. F = ',uF,'\n')
            sink()
            
            # Interpolate on regular grid
            xsl  = downSample(wl0, xs0, reso = reso)
            wl   = xsl$wl
            xs   = xsl$xs
            
            # Remove tailing zeroes
            first = which(xs != 0)[1]
            last  = length(xs) - which(rev(xs) != 0)[1] + 1
            wavlXS = wl[first:last]
            valXS  = xs[first:last]
          
            # print(paste0('Got XS for ',sp))
            # print(str(wavlXS))
            
            # Get BRs ####
            # Get descriptors
            type  = photoDB()$BR_SOURCE[iReac]
            chans = which(photoDB()$REACTANTS == photoDB()$REACTANTS[iReac])
            nBR   = length(chans)

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
              wl1 = xsl$wl
              xs1 = xsl$xs
              # Remove tailing zeroes
              first = which(xs1 != 0)[1]
              last  = length(xs1)-which(rev(xs1) !=0)[1] + 1
              wl1   = wl1[first:last]
              
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
              for (i in 1:np)
                qy0[ ,i] = downSample(wl, S[, i+i0]/xs, 
                                     reso = reso)$xs[first:last]
              wavlBR = wl1
              files = paste0(sp, '.dat')
              
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
                             header = TRUE, check.names = FALSE)
              wl  = S[, 1]
              xs  = S[, 2]
              xsl = downSample(wl, xs, reso = reso)
              wl1 = xsl$wl
              xs1 = xsl$xs
              # Remove tailing zeroes
              first = which(xs1 != 0)[1]
              last  = length(xs1)-which(rev(xs1) !=0)[1] + 1
              wl1    = wl1[first:last]
              
              qy0  = matrix(0, ncol = length(files), nrow = length(wl1))
              for (i in seq_along(files)) {
                S = read.table(file = file.path(source,files[i]),
                               header = TRUE, check.names = FALSE)
                wl  = S[, 1]
                xs  = S[, 2]
                qy0[ ,i] = downSample(wl, xs, reso = reso)$xs[first:last]
              }
              wavlBR = wl1
            }
            # print(paste0('Got BRs for ',sp))
            # print(tail(qy))
            
            sink(file = photoSourceFile, append = TRUE)
            cat(paste0('   BR : ',source,'/',files, collapse = '\n'),'\n')
            sink()
            
            ## Define common wavl range for XS and qy
            ## BRs outside of XS range are useless...
            lims   = range(wavlXS)
            selBR  = wavlBR >= lims[1] & wavlBR <= lims[2]
            wavlBR = wavlBR[selBR]
            qy0    = qy0[selBR,]
            # print(str(wavlBR))
            # print(rowSums(qy))
            
            # Expand BR range to XS if necessary
            lims = range(wavlBR)
            sel  = wavlXS >= lims[1] & wavlXS <= lims[2] 
            if(sum(!sel) != 0 ) {
              qy1 = matrix(0, ncol = ncol(qy0), nrow = length(wavlXS))
              qy1[sel,] = qy0
              # print(sel)
              if(!sel[1]) {
                # print('Extrap. left')
                # Extrapolate towards short wavl
                i0 = which(sel)[1]
                for(i in 1:(i0-1))
                  qy1[i,] = qy0[1,]
              }
              if(!sel[length(sel)]) {
                # print('Extrap. right')
                # Extrapolate towards long wavl
                i0 = length(sel) - which(rev(sel))[1] + 1
                for(i in (i0+1):nrow(qy1))
                  qy1[i,] = qy0[nrow(qy0),]
              }
              # print(qy1)
              qy0 = qy1
              wavlBR = wavlXS
            }
            # print(rowSums(qy))
            
            
            # Sample XS ####
            incProgress(1/len, 
                        detail = paste0(' XS ',sp,'/',type,'/',reso,' nm'))
            if(!input$photoSampleCheck) { 
              for (iMC in 0:nMC) {
                prefix = paste0(sprintf('%04i', iMC), '_')
                # Systematic perturbation of XS
                if(iMC == 0) {
                  # Nominal run
                  Frnd = 1
                } else {
                  rnd  =  truncnorm::rtruncnorm(1,-3,3,0,1) # Avoid outliers
                  Frnd = exp( log(uF) * rnd )
                }
                write.table(
                  cbind(wavlXS, valXS * Frnd),
                  sep = ' ',
                  row.names = FALSE,
                  col.names = FALSE,
                  file = gzfile(
                    file.path(target, paste0(prefix, 'se', sp, '.dat.gz'))
                  )
                )
              }
            }
           
            # Sample BRs ####
            incProgress(1/len, 
                        detail = paste0(' BRs ',sp,'/',type,'/',reso,' nm'))
            if(!input$photoSampleCheck) { 
              if(nBR == 1) {
                # A single channel: no uncertainty in BR
                qySample = array(
                  data = 1,
                  dim = c(1+nMC, ncol(qy), nBR)
                )
                
              } else {
                
                qy = qy0 / rowSums(qy0)
                
                # Generate Diri-based samples at each wavelength
                ionic = c()
                for (i in 1:nBR)
                  ionic[i] = 'E' %in% getSpecies(photoDB()$PRODUCTS[chans[i]])
                
                if (input$photoSampleFlatTree) {
                  
                  qySample = flatSample(
                    qy0, 
                    ionic = ionic,
                    ru  = c(photoRuBRN,photoRuBRI),
                    nMC = nMC,
                    eps = photoEps)
                  
                } else {
                  
                  # Nested sampling
                  qySample = hierSample(
                    qy0, 
                    ionic = ionic,
                    ru = c(photoRuBRNI, photoRuBRN, photoRuBRI),
                    nMC = nMC,
                    eps = photoEps)
                }
                
                # Rearrange samples for better wavelength continuity
                if(input$photoSampleBRArrange)
                  qySample = arrangeSample(qySample)
                
              }
              
              # Nominal run
              iMC = 0
              prefix = paste0(sprintf('%04i', iMC), '_')
              # save
              for (i in 1:nBR) {
                fileOut = file.path(
                  target, 
                  paste0(prefix,'qy',sp,'_',i,'.dat.gz')
                )
                write.table(
                  cbind(wavlBR, qy[, i]),
                  sep = ' ',
                  row.names = FALSE,
                  col.names = FALSE,
                  file = gzfile(fileOut)
                )
              }
              for(iMC in 1:nMC) {
                prefix = paste0(sprintf('%04i',iMC),'_')
                for(i in 1:nBR) {
                  fileOut = file.path(
                    target, 
                    paste0(prefix,'qy',sp,'_',i,'.dat.gz')
                  )
                  write.table(
                    cbind(wavlBR,qySample[iMC,,i]),
                    sep=' ', row.names=FALSE, col.names=FALSE,
                    file = gzfile(fileOut)
                  )
                }
              }
            }
          }
        }
      })
    
  })

# Plots ####

shiny::observe({
  req(photoDB())
  reacs = photoDB()$REACTANTS
  tag   = reacs[ photoDB()$CHANNEL == 1 ] 
  nums = 0:length(tag)  
  names(nums) = c("0",tag)
  photoSampleReacsFiltered(nums)
  shiny::updateSelectInput(
    session  = session,
    "photoSampleReaction",
    choices  = nums,
    selected = "0"
  )    
  
})

shiny::observeEvent(
  input$photoSampleMinus,
  {
    iReac = as.numeric(input$photoSampleReaction)
    if(iReac > 0) {
      iReac = iReac - 1
      shiny::updateSelectInput(
        session = session,
        "photoSampleReaction",
        selected = iReac
      )  
    }
  })

shiny::observeEvent(
  input$photoSamplePlus,
  {
    req(photoSampleReacsFiltered())
    reacs = photoSampleReacsFiltered()
    
    iReac = as.numeric(input$photoSampleReaction)
    if(iReac < length(reacs)-1) {
      iReac = iReac + 1
      shiny::updateSelectInput(
        session = session,
        "photoSampleReaction",
        selected = iReac
      )
    }
  })

## XS ####
output$plotSampleXS = renderPlot({
  
  reso = input$photoSampleXSReso
  photoSampleDir = paste0(photoPublic, '_', photoEditOrigVersion())
  source_dir = file.path(photoSampleDir, paste0(reso,'nm'))
  req(dir.exists(source_dir))
  
  req(input$photoSampleReaction)
  iReac = as.numeric(input$photoSampleReaction)
  req(iReac != 0)

  nMC = as.numeric(input$photoSamplePlotSize)  
  
  sp = getSpecies(names(photoSampleReacsFiltered())[iReac+1])[1]
  iMC = 0
  prefix=paste0(sprintf('%04i',iMC),'_')
  file = file.path(source_dir,paste0(prefix,'se',sp,'.dat.gz'))
  x = read.table(
    file,
    col.names=1:10,
    header=FALSE, 
    fill=TRUE)
  wavl = x[,1]
  X = matrix(NA,nrow=nrow(x),ncol=nMC+1)
  X[,1] = x[,2]
  for (iMC in 1:nMC) {
    prefix=paste0(sprintf('%04i',iMC),'_')
    file = file.path(source_dir,paste0(prefix,'se',sp,'.dat.gz'))
    if(!file.exists(file))
      next
    x = read.table(
      file,
      col.names=1:10,
      header=FALSE, 
      fill=TRUE)
    X[,iMC+1] = x[,2]
  }
  
  log = ''
  ylim = c(0, 1.1*max(X[,1]))
  if(input$photoSampleXSLog){
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
    wavl, X,
    type = 'l', 
    lwd  = 3,
    lty  = 1,
    xaxs = 'i',
    xlab = 'Wavelength [nm]',
    xlim = input$photoSampleWLPlotRange,
    yaxs = 'i',
    log  = log,
    ylim = ylim,
    ylab = expression(paste('Cross-section [', cm ^ 2, ']')),
    col  = gPars$cols_tr[5],
    main = paste0(sp,' / ',reso,'nm / ',photoEditOrigVersion())
  )
  grid()
  lines(wavl,X[,1],col=gPars$cols[2])
  box()
},
height = plotHeight)

## BRs ####
output$plotSampleBR = renderPlot({

  req(input$photoSampleReaction)
  iReac = as.numeric(input$photoSampleReaction)
  req(iReac != 0)
  
  reso = input$photoSampleXSReso
  photoSampleDir = paste0(photoPublic, '_', photoEditOrigVersion())
  source_dir = file.path(photoSampleDir, paste0(reso,'nm'))
  req(dir.exists(source_dir))
  
  sp = getSpecies(names(photoSampleReacsFiltered())[iReac+1])[1]
  
  # Search BR files for sp
  prefix  =  paste0(sprintf('%04i',0),'_')
  pattern = paste0(prefix,'qy',sp,'_')
  files   = list.files(path = source_dir, pattern = pattern)
  nBR     = length(files)
  
  # Check nb of available samples
  nMC = as.numeric(input$photoSamplePlotSize)
  pattern = paste0('qy',sp,'_1')
  files   = list.files(path = source_dir, pattern = pattern)
  nMCa = length(files) - 1
  if(nMCa < nMC) {
    id = shiny::showNotification(
      strong(paste0(
        'Number of available samples (',
        nMCa,
        ') smaller than requested !')
      ),
      closeButton = TRUE,
      duration = NULL,
      type = 'warning'
    )
    # Pursue with available samples
    nMC = nMCa
  }

  # Get data
  wavl = read.table(
    file.path(source_dir,paste0(prefix,'qy',sp,'_1.dat.gz'))
  )[,1]

  sampleBR0 = matrix(0, nrow = length(wavl), ncol = nBR)
  prefix=paste0(sprintf('%04i',0),'_')
  for (j in 1:nBR) {
    x = read.table(
      file.path(source_dir,paste0(prefix,'qy',sp,'_',j,'.dat.gz'))
    )
    sampleBR0[,j] = x[,2]
  }
  sampleBR = array(data = 0,dim  = c(nMC,length(wavl),nBR) )
  for(iMC in 1:nMC) {
    prefix=paste0(sprintf('%04i',iMC),'_')
    for (j in 1:nBR) {
      x = read.table(
        file.path(source_dir,paste0(prefix,'qy',sp,'_',j,'.dat.gz'))
      )
      sampleBR[iMC,,j] = x[,2]
    }
  }

  reac = which(
    photoDB()$REACTANTS == names(photoSampleReacsFiltered())[iReac+1]
  )[1]
  channels = reac:(reac+nBR-1)
  
  cols    = rep(gPars$cols,2)
  cols_tr = rep(gPars$cols_tr,2)
  lty     = c(rep(1,length(gPars$cols)),rep(2,length(gPars$cols)))
  
  ylab = 'Branching ratios'
  ylim = c(-0.01, 1.4)
  
  if (input$photoSampleBRDisplay == 0 | nBR == 1) {
    qy0 = sampleBR0
    qy  = sampleBR
    leg = photoDB()$PRODUCTS[channels]
    
    
  } else if (input$photoSampleBRDisplay == 1) {
    ionic = c()
    for (i in 1:nBR)
      ionic[i] = 'E' %in% getSpecies(photoDB()$PRODUCTS[channels[i]])
    qy0 = matrix(data = 0, nrow = length(wavl), ncol = 2)
    qy0[, 1] = rowSums(sampleBR0[, !ionic, drop = FALSE])
    qy0[, 2] = rowSums(sampleBR0[,  ionic, drop = FALSE])
    qy = array(data = 0, dim  = c(nMC, length(wavl), 2))
    qy[, , 1] = rowSums(sampleBR[, , !ionic, drop = FALSE], dims = 2)
    qy[, , 2] = rowSums(sampleBR[, ,  ionic, drop = FALSE], dims = 2)
    leg = c("Neutrals", "Ions")
    
  } else if (input$photoSampleBRDisplay == 2) {
    qy0 = rowSums(sampleBR0)
    qy = array(data = 0, dim  = c(nMC, length(wavl), 1))
    qy[, , 1] = rowSums(sampleBR, dims = 2)
    leg = c("Sum-to-one")
    
  } else if (input$photoSampleBRDisplay == 3) {
    mu  = apply(sampleBR, c(2,3), mean, na.rm = TRUE)
    sig = apply(sampleBR, c(2,3), sd  , na.rm = TRUE)
    mu[mu==0] = 1
    qy0 = sig / mu
    leg = photoDB()$PRODUCTS[channels]
    ylab = 'Relative uncertainty'
    ylim = c(-0.01,1.4*max(qy0, na.rm = TRUE))
    
  } else if (input$photoSampleBRDisplay == 4) {
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
    qy0 = matrix(data = 0, nrow = length(wavl), ncol = 3)
    qy0[,1] = apply(qy[, !ionic, drop = FALSE], 1, nzMean)
    qy0[,2] = apply(qy[,  ionic, drop = FALSE], 1, nzMean)
    qy0[,3] = apply(qy, 1, nzMean)
    leg = c("Neutrals", "Ions", "Global")
    ylab = 'Mean relative uncertainty'
    ylim = c(-0.01,1.4*max(qy0, na.rm = TRUE))
  }
  
  par(mar = c(4, 4, 2, 1),
      mgp = gPars$mgp,
      tcl = gPars$tcl,
      lwd = gPars$lwd,
      pty = 'm',
      cex = 1.5)
  
  matplot(
    wavl, qy0,
    type = 'l', 
    lwd  = 3,
    xaxs = 'i',
    xlab = 'Wavelength [nm]',
    xlim = input$photoSampleWLPlotRange,
    yaxs = 'i',
    ylim = ylim,
    ylab = ylab,
    col  = cols,
    lty  = lty,
    main = ''
  )
  grid()
  
  if (input$photoSampleBRDisplay < 3) {
    for(iMC in 1:nMC) {
      matlines(
        wavl, qy[iMC,,],
        lwd  = 3,
        col  = cols_tr,
        lty  = lty
      )
    }
    matlines(
      wavl, qy0, 
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
height = plotHeight)
