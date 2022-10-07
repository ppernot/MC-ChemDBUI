photoSampleReacsFiltered = shiny::reactiveVal()

# Sample ####
observeEvent(
  input$photoSampleBtn,
  {
    req(photoDB())
    
    nMC = 1 + as.numeric(input$photoSampleSize) # Account for nominal
    
    sortBR = input$photoBRSampleSort
    
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
        fileList = list.files(path = target, full.names = TRUE)
        file.remove(fileList)
      }
    }
    
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
              qy = matrix(0, ncol = np, nrow = length(wl1))
              for (i in 1:np)
                qy[ ,i] = downSample(wl, S[, i+i0]/xs, 
                                     reso = reso)$xs[first:last]
              wavlBR = wl1
              
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
              
              qy  = matrix(0, ncol = length(files), nrow = length(wl1))
              for (i in seq_along(files)) {
                S = read.table(file = file.path(source,files[i]),
                               header = TRUE, check.names = FALSE)
                wl  = S[, 1]
                xs  = S[, 2]
                qy[ ,i] = downSample(wl, xs, reso = reso)$xs[first:last]
              }
              wavlBR = wl1
            }
            # print(paste0('Got BRs for ',sp))
            # print(tail(qy))
            
            ## Define common wavl range for XS and qy
            ## BRs outside of XS range are useless...
            lims   = range(wavlXS)
            selBR  = wavlBR >= lims[1] & wavlBR <= lims[2]
            wavlBR = wavlBR[selBR]
            qy     = qy[selBR,]
            # print(str(wavlBR))
            # print(rowSums(qy))
            
            # Expand BR range to XS if necessary
            lims = range(wavlBR)
            sel  = wavlXS >= lims[1] & wavlXS <= lims[2] 
            if(sum(!sel) != 0 ) {
              qy1 = matrix(0, ncol = ncol(qy), nrow = length(wavlXS))
              qy1[sel,] = qy
              # print(sel)
              if(!sel[1]) {
                # print('Extrap. left')
                # Extrapolate towards short wavl
                i0 = which(sel)[1]
                for(i in 1:(i0-1))
                  qy1[i,] = qy[1,]
              }
              if(!sel[length(sel)]) {
                # print('Extrap. right')
                # Extrapolate towards long wavl
                i0 = length(sel) - which(rev(sel))[1] + 1
                for(i in (i0+1):nrow(qy1))
                  qy1[i,] = qy[nrow(qy),]
              }
              # print(qy1)
              qy = qy1
              wavlBR = wavlXS
            }
            # print(rowSums(qy))
            
            
            # Sample XS ####
            incProgress(1/len, 
                        detail = paste0(' XS ',sp,'/',type,'/',reso,' nm'))
            for (iMC in 0:nMC) {
              prefix = paste0(sprintf('%04i', iMC), '_')
              # Systematic perturbation of XS
              if(iMC == 0) {
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
           
            # Sample BRs ####
            incProgress(1/len, 
                        detail = paste0(' BRs ',sp,'/',type,'/',reso,' nm'))
            if(nBR == 1) {
              # A single channel: no uncertainty in BR
              qySample = array(
                data = 1,
                dim = c(nMC, ncol(qy), nBR)
              )
              
            } else {
              
              qy = qy / rowSums(qy)
              
              # Generate Diri-based samples at each wavelength
              ionic = c()
              for (i in 1:nBR)
                ionic[i] = 'E' %in% getSpecies(photoDB()$PRODUCTS[chans[i]])
              
              if ( (sum(ionic) * sum(!ionic)) == 0) {
                # Diri sampling
                qySample = diriSample(
                  qy,
                  ru = ifelse( sum(ionic) == 0, photoRuBRN, photoRuBRI),
                  nMC = nMC,
                  eps = photoEps,
                  sortBR = sortBR)
              } else {
                # Nested sampling
                qySample = hierSample(
                  qy, 
                  ionic = ionic,
                  ru = c(photoRuBRNI, photoRuBRN, photoRuBRI),
                  nMC = nMC,
                  eps = photoEps,
                  sortBR = sortBR)
              }
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
  
  nMC = as.numeric(input$photoSamplePlotSize)  
  
  sp = getSpecies(names(photoSampleReacsFiltered())[iReac+1])[1]
  
  # Search BR files for sp
  prefix  =  paste0(sprintf('%04i',0),'_')
  pattern = paste0(prefix,'qy',sp,'_')
  files   = list.files(path = source_dir, pattern = pattern)
  nBR     = length(files)

  wavl = read.table(
    file.path(source_dir,paste0(prefix,'qy',sp,'_1.dat.gz'))
  )[,1]

  qySample = array(data = 0,dim  = c(nMC+1,length(wavl),nBR) )
  for(iMC in 0:nMC) {
    prefix=paste0(sprintf('%04i',iMC),'_')
    for (j in 1:nBR) {
      x = read.table(
        file.path(source_dir,paste0(prefix,'qy',sp,'_',j,'.dat.gz'))
      )
      qySample[iMC+1,,j] = x[,2]
    }
  }

  reac = which(
    photoDB()$REACTANTS == names(photoSampleReacsFiltered())[iReac+1]
  )[1]
  
  if (input$photoSampleBRDisplay == 0 | nBR == 1) {
    qy = qySample
    leg = photoDB()$PRODUCTS[reac:(reac+nBR-1)]
    
  } else if (input$photoSampleBRDisplay == 1) {
    ionic = c()
    for (i in 1:nBR)
      ionic[i] = 'E' %in% getSpecies(photoDB()$PRODUCTS[reac + i - 1])
    qy = array(data = 0, dim  = c(nMC + 1, length(wavl), 2))
    qy[, , 1] = rowSums(qySample[, , !ionic, drop = FALSE], dims = 2)
    qy[, , 2] = rowSums(qySample[, ,  ionic, drop = FALSE], dims = 2)
    leg = c("Neutrals","Ions")
    
  } else {
    qy = array(data = 0, dim  = c(nMC + 1, length(wavl), 1))
    qy[, , 1] = rowSums(qySample, dims = 2)
    leg = c("Sum-to-one")
  }
  
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
      wavl, qy[1,,],
      type = 'l', 
      lwd  = 3,
      xaxs = 'i',
      xlab = 'Wavelength [nm]',
      xlim = input$photoSampleWLPlotRange,
      yaxs = 'i',
      ylim = c(-0.01, 1.2),
      ylab = paste('Branching ratios'),
      col  = cols_tr,
      lty  = lty,
      main = ""
    )
    grid()
    for(iMC in 1:nMC) {
      matlines(
        wavl, qy[iMC+1,,],
        lwd  = 3,
        col  = cols_tr,
        lty  = lty
      )
    }
    # Overdraw nominal curves
    matlines(
      wavl, qy[1,,], 
      lwd = 4, 
      col = cols, 
      lty = lty
    )
    
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

# Statistics ####
output$photoStats = renderPrint({
  # req(collateDone())
  # 
  # id = shiny::showNotification(
  #   h4('Computing statistics...'),
  #   closeButton = TRUE,
  #   duration = NULL,
  #   type = 'message'
  # )
  # 
  # fileName = file.path(
  #   paste0(ionsPublic,'_',ionsEditOrigVersion()),
  #   'run_0000.csv')
  # scheme   = as.data.frame(
  #   data.table::fread(file = fileName, header = FALSE, sep = ';')
  # )
  # nbReac   = nrow(scheme)
  # reactants = products = params = type = reacTagFull = list()
  # for (i in 1:nbReac) {
  #   line = replace(scheme[i, ], scheme[i, ] == " ", "")
  #   line = replace(line, is.na(line), "")
  #   reactants[[i]] = line[1:3][line[1:3]!=""]
  #   products[[i]]  = line[4:7][line[4:7]!=""]
  #   type[[i]]      = line[ncol(scheme)]
  #   reacTagFull[[i]] = paste0(
  #     paste(reactants[[i]], collapse = ' + '), 
  #     ' --> ', 
  #     paste(products[[i]], collapse = ' + ')
  #   )
  # }
  # 
  # # Build species list from reactants and products
  # species   = unique(unlist(c(reactants,products)))
  # nbSpecies = length(species)
  # reacts    = unique(unlist(reactants))
  # nbReacts  = length(reacts)
  # prods     = unique(unlist(products))
  # nbProds   = length(prods)
  # 
  # cat('Nb reactions       = ', nbReac   ,'\n')
  # cat('Nb species         = ', nbSpecies,'\n\n')
  # 
  # # Nature of species
  # types = c("neutrals","ions")
  # chems = c("hydrocarbons","N-bearing","O-bearing")
  # heavy = c("C0","C1","C2","C3","C4","C5","C6","Cmore")
  # resu = matrix(0,
  #               nrow = 1 + length(chems),
  #               ncol = 2 + length(types))
  # rownames(resu) = c(chems,'total')
  # colnames(resu) = c(types,'radicals','total')
  # for (type in types)
  #   resu['total',type] = sum(
  #     selectSpecies(species,c(type,chems,heavy)))
  # resu['total','radicals'] = sum(
  #   selectSpecies(species,c(types,'radicals',chems,heavy)))
  # for (chem in chems) {
  #   resu[chem,'total'] = sum(
  #     selectSpecies(species,c(types,chem,heavy)))
  #   resu[chem,'radicals'] = sum(
  #     selectSpecies(species,c(types,'radicals',chem,heavy)))
  # }
  # for (type in types)
  #   for (chem in chems)
  #     resu[chem,type] = sum(
  #       selectSpecies(species,c(type,chem,heavy)))
  # resu['total','total'] = sum(resu[,'total'])
  # cat('Compositions (excl. dummies)\n')
  # cat('----------------------------\n\n')
  # print(resu)
  # cat('\n\n')
  # 
  # # Loss and prod
  # L = R = matrix(0,ncol=nbSpecies,nrow=nbReac)
  # for (m in 1:nbReac) {
  #   reac = unlist(reactants[m])
  #   prod = unlist(products[m] )
  #   for (n in 1:nbSpecies) {
  #     search=species[n]
  #     L[m,n] = length(which( search == reac )) # Loss
  #     R[m,n] = length(which( search == prod )) # Prod
  #   }
  # }
  # selLossless = colSums(L) == 0
  # nbLossless = sum(selLossless)
  # if(nbLossless != 0)
  #   lossless = sort(species[selLossless])
  # selProdless = colSums(R) == 0
  # nbProdless = sum(selProdless)
  # if(nbProdless != 0)
  #   prodless = sort(species[selProdless])
  # 
  # cat('Reactants   : ',nbReacts,'\n')
  # cat('Products    : ',nbProds,'\n\n')
  # cat('Number of lossless species = ',nbLossless)
  # if(nbLossless != 0) 
  #   cat(strwrap(paste0(lossless, collapse = ', '),width=60,prefix='\n'),'\n\n')
  # cat('Number of prodless species = ',nbProdless)
  # if(nbProdless != 0) 
  #   cat(strwrap(paste0(prodless, collapse = ', '),width=60,prefix='\n'),'\n\n')
  # 
  # id0=shiny::removeNotification(id)
})