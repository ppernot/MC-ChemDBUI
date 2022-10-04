observeEvent(
  input$photoSampleBtn,
  {
    req(photoDB())
    
    sampleSize = 1 + as.numeric(input$photoSampleSize) # Account for nominal
    
    # Cross-sections ####
    for (iReac in seq_along(photoDB()$REACTANTS)) {
      
      # Treat first channel only
      if(photoDB()$CHANNEL[iReac] != 1)
        next
      
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
        wl   = xsl$wavelength
        xs   = xsl$photoabsorption
        uF0  = xsl$uncF # Priority on DB parameter
        
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
      }
      
      if(!is.null(uF0))
        uF = as.numeric(uF0)
      else
        uF = photoDefaultuF
      
      for (reso in photoXSResolutions){
        target = file.path(photoPublic,photoEditOrigVersion(),
                           paste0(reso,'nm'))
        
        # Interpolate on regular grid
        xsl  = downSample(wl, xs, reso = reso)
        wl   = xsl$wl
        xs   = xsl$xs
        
        # Remove tailing zeroes
        first = which(xs != 0)[1]
        last  = length(xs) - which(rev(xs) != 0)[1] + 1
        wl   = wl[first:last]
        xs   = xs[first:last]
        
        for (iMC in 1:sampleSize) {
          prefix = paste0(sprintf('%04i', iMC), '_')
          # Systematic perturbation of XS
          if(iMC == 1) {
            Frnd = 1
          } else {
            rnd =  truncnorm::rtruncnorm(1,-3,3,0,1) # Avoid outliers
            # Frnd = rlnorm(1, meanlog = 0, sdlog = log(uF))
            Frnd = exp( log(uF) * rnd )
          }
          write.table(
            cbind(wl, xs * Frnd),
            sep = ' ',
            row.names = FALSE,
            col.names = FALSE,
            file = gzfile(
              file.path(target, paste0(prefix, 'se', sp, '.dat.gz'))
            )
          )
        }
      }
    }
  })

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
