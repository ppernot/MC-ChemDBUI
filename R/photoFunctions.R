# Misc functions ####
# getSpecies = function(
#   file = 'PhotoScheme.dat'
# ) {
#   nbReac = 0
#   reactants = products = params = type = orig = locnum = list()
# 
#   scheme  = read.fwf(file = file, widths = rep(11, 12))
#   scheme  = t(apply(scheme, 1, function(x)
#     gsub(" ", "", x)))
#   for (i in 1:nrow(scheme)) {
#     nbReac = nbReac + 1
#     terms = scheme[i, 1:2]
#     reactants[[nbReac]] = terms[!is.na(terms) &
#                                   terms != "" & terms != "HV"]
#     terms = scheme[i, 3:6]
#     products[[nbReac]]  = terms[!is.na(terms) &
#                                   terms != "" & terms != "HV"]
#     terms = scheme[i, 7:12]
#     params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
#     type[[nbReac]]      = 'photo'
#     locnum[[nbReac]]    = i
#     orig[[nbReac]]      = file
#   }
#   species = levels(as.factor(unlist(reactants)))
#   nbSpecies = length(species)
# 
#   # Build loss matrix
#   L = matrix(0, ncol = nbSpecies, nrow = nbReac)
#   for (m in 1:nbReac) {
#     reac = unlist(reactants[m])
#     for (n in 1:nbSpecies) {
#       search = species[n]
#       L[m, n] = length(which(search == reac)) # Loss
#     }
#   }
#   colnames(L) = species
# 
#   return(
#     list(
#       species = species,
#       L = L,
#       params = params,
#       products = products
#     )
#   )
# }

getXShdf5 = function (
  species,
  source_dir = './Leiden/',
  getxs = TRUE # Set to FALSE to get only uncF
) {
  # Misc. info
  info = read.csv(
    file = file.path(source_dir, 'cross_section_properties.csv'),
    header = TRUE,
    skip = 18,
    stringsAsFactors = FALSE
  )
  
  # Get out if species not available in dataset
  if (!(species %in% trimws(info$species)))
    return(NULL)

  # Species OK: proceed
  rownames(info) = trimws(info$species)
  uncCode =  trimws(info[species, 'cross_sec_unc'])
  uF = switch(
    uncCode,
    'A+' = 1.2,
    'A'  = 1.3,
    'B'  = 2  ,
    'C'  = 10
  )

  xs = list()
  xs[['uncF']] = uF

  if(getxs) {
    # Cross sections
    file = file.path(source_dir, 'cross_sections', 
                     species, paste0( species, '.hdf5') )
    fh5 = hdf5r::H5File$new(file, mode = 'r')

    for (key in c('photoabsorption',
                  'photodissociation',
                  'photoionisation',
                  'wavelength')) {
      xs[[key]] = fh5[[key]]$read()
    }
  }

  return(xs)
}
downSample = function(wl,xs,reso = 1) {
  # Define new grid
  lims = round(range(wl))
  wl1  = seq(min(lims), max(lims), by = reso)
  nl1  = length(wl1)

  # Compute increments on original grid
  dwl = diff(wl)
  ## Remove intervals of null width
  sel = dwl != 0
  dwl = dwl[sel]
  sel = c(sel, FALSE)
  pI = xs[sel] * dwl
  wl = wl[sel]

  # Interpolate from [ir]regular to regular grid
  p1 = rep(NA, length(wl1) - 1)
  i = 0
  for (il in wl1[1:(nl1 - 1)]) {
    i = i + 1
    # Which points within interval
    sel = wl >= il - 0.5 * reso &
      wl <  il + 0.5 * reso
    if (sum(sel) == 0) {
      # No point in interval: increase upper limit
      sel1 = which(wl > il + 0.5 * reso)[1]
      if (is.na(sel1)) {
        # No point within full range: contribution is null
        p1[i] = 0
      } else {
        # Found some points
        if (i == 1) {
          # If first point on grid assign first value
          p1[i] = pI[sel1] / dwl[sel1]
        } else {
          # If not first point, linear interpolation from previous value
          x0 = il - reso
          v0 = p1[i - 1]
          x1 = wl[sel1]
          v1 = pI[sel1] / dwl[sel1]
          p1[i] = v0 + (v1 - v0) / (x1 - x0) * reso
        }
      }
    } else {
      # At least one point in regular interval: sum contributions
      p1[i] = sum(pI[sel]) / sum(dwl[sel])
    }
  }
  # Remove last point
  wl1 = wl1[1:(length(wl1) - 1)]

  return(list(wl = wl1, xs = p1))
}
# nds = function(ns,dist) {
#   command=paste("echo ",ns," '",dist,"' | ./Bin/Rnested.x")
#   # quotes around dist avoid shell interpretation
#   tc=textConnection(system(command,intern=T))
#   liste=scan(tc,quiet=TRUE)
#   close(tc)
#   nleaves=liste[1]
#   nlast=nleaves*ns+1
#   nds=matrix(liste[2:nlast],ncol=nleaves,byrow=T)
#   return(nds)
# }
gamDiri = function(x,r) { # Eq.9 in Plessis2010
  x  = x / sum(x)  # Ensure normalization
  return( 1 / r^2 * (sum(x*(1-x)) / sum(x*sqrt(x*(1-x))))^2 - 1 )
}
gamDiriGM = function(x,r) { # Geom Mean gamma
  x = x / sum(x)
  p = prod( (1 - x - x*r^2) / (x*r^2))
  return( p^(1/length(x)) )
}
threshComp = function(X,r) {
  # Threshold to avoid stray Diri samples
  # Might be problematic in dimensions > 2 and for large r !!!
  amin = 1.1 * r^2 / (1 + r^2) 
  sel  = X < amin
  if(sum(sel) > 0) {
    X[ sel] = amin
    X[!sel] = X[!sel] / sum(X[!sel]) * (1 - sum(sel)*amin)
  }
  return(X)
}
diriSample0 = function(br, ru, nMC, eps, newGam = TRUE, fDirg = TRUE) {
  # Flat sampling by Diri or Dirg
  # 
  qySample = matrix(0, nrow = nMC, ncol = length(br))

  # Count non-zero channels
  br = br / sum(br)
  br[br <= eps] = 0
  br = br / sum(br)
  sel_nz = br != 0

  if( sum(sel_nz) <= 1 ) {
    # 1 channel: no uncertainty
    qySample[,sel_nz] = 1

  } else {
    if(newGam) {
      # Dirg
      bru = br[sel_nz] * ru
      if(fDirg)
        bru = bru * fScaleDirg(sum(sel_nz))
      stringBR = paste0(
        'Dirg(',
        paste0(br[sel_nz], collapse = ','),
        ';',
        paste0(bru, collapse = ','),
        ')')
      
    } else {
      X = br[sel_nz]
      X = threshComp( X, ru)
      gamma = gamDiriGM(X, ru)
      # Dirichlet
      stringBR = paste0(
        'Diri(',
        paste0(X, collapse = ','),
        ';',
        gamma,
        ')')
    }
    # Sample by Nested.x
    qySample[,sel_nz] = nds(nMC, stringBR)
  }

  return(qySample)
}
diriSample = function(qy, ru = 0.1, nMC = 500, eps = 1e-4, 
                      newGam = TRUE, fDirg = TRUE) {

  nc = ncol(qy)
  nw = nrow(qy)
  qySample = array(
    data = 0,
    dim  = c(nMC,nw,nc)
  )

  for (il in 1:nw)
    qySample[ , il, ] = diriSample0(qy[il,], ru, nMC, eps, newGam, fDirg)

  return(qySample)
}
defDir = function(n){
  if(n<=1)
    ''
  else
    paste0('*Diun(',n,')')
}
hierSample  = function(qy, ionic, ru = c(0.1,0.1,0.1), 
                       nMC = 500, eps = 1e-4, 
                       newGam = TRUE, fDirg = TRUE) {
  # Nested sampling when ionic and !ionic channels present
  # *** Treat only non-zero channels ***
  nc = ncol(qy); nw = nrow(qy)
  qySample = array(
    data = 0,
    dim  = c(nMC,nw,nc)
  )
  
  for (il in 1:nw) {
    
    br = qy[il,]
    br[br <= eps] = 0 # Threshold
    br = br / sum(br) # Renormalize

    brNeu = sum(br[!ionic])
    indxN = NULL
    if (brNeu > 0) {
      brN = br[!ionic] / brNeu # Renormalize for neutrals sub-tree
      sel_nzN = brN > 0
      indxN = which(!ionic)[sel_nzN]
      stringBRN = defDir(sum(sel_nzN))
      if (sum(sel_nzN) > 1) {
        if (newGam) {
          r = ru[2]
          brNu = r * brN
          if(fDirg)
            brNu = brNu * fScaleDirg(sum(sel_nzN))
          stringBRN = paste0(
            '*Dirg(',
            paste0(brN[sel_nzN], collapse = ','),
            ';',
            paste0(brNu[sel_nzN], collapse = ','),
            ')')
        } else {
          r = ru[2]
          X = brN[sel_nzN]
          X = threshComp( X, r)
          gamma = gamDiriGM(X, r)
          stringBRN = paste0(
            '*Diri(',
            paste0(X, collapse = ','), ';',
            gamma, ')')
        }
      }
    }
    
    brIon = sum(br[ionic])
    indxI = NULL
    if (brIon > 0) {
      brI  = br[ionic] / brIon # Renormalize for ionic sub-tree
      sel_nzI = brI > 0
      indxI  = which(ionic)[sel_nzI]
      stringBRI = defDir(sum(sel_nzI))
      if (sum(sel_nzI) > 1) {
        if (newGam) {
          r = ru[3]
          brIu = r * brI
          if(fDirg)
            brIu = brIu * fScaleDirg(sum(sel_nzI))
          stringBRI = paste0(
            '*Dirg(',
            paste0(brI[sel_nzI], collapse = ','),
            ';',
            paste0(brIu[sel_nzI], collapse = ','),
            ')')
        } else {
          r = ru[3]
          X = brI[sel_nzI]
          X = threshComp( X, r)
          gamma = gamDiriGM(X, r)
          stringBRI = paste0(
            '*Diri(',
            paste0(X, collapse = ','), ';',
            gamma, ')')
        }
      }
    }
    
    # Collate sampled channels indices
    ret_nz = c()
    if(!is.null(indxN))
      ret_nz = c(ret_nz,indxN)
    if(!is.null(indxI))
      ret_nz = c(ret_nz,indxI)
    
    # Sampling
    if(length(ret_nz) == 1) {
      qySample[, il, ret_nz] = 1.0
      
    } else {
      if(brIon == 0)
        stringBR = substring(stringBRN,2) # Remove leading star
      
      else if(brNeu == 0)
        stringBR = substring(stringBRI,2) # Remove leading star
      
      else {
        r = ru[1] 
        X = c(brNeu,brIon)
        if(newGam) {
          # Interpret r as absolute uncertainty
          bru = r * X
          if(fDirg)
            bru = bru * fScaleDirg(2)
          stringBR = paste0(
            'Dirg(',X[1],stringBRN,',',X[2],stringBRI,';',
            paste0(bru, collapse = ','),')'
          )
        } else {
          X = threshComp( X, r)
          gammaNI = gamDiri(X, r)
          # gammaNI = gamDiriGM(X, r)
          stringBR = paste0(
            'Diri(',X[1],stringBRN,',',X[2],stringBRI,';',
            gammaNI,')'
          )
        } 
        
      }
      qySample[, il, ret_nz] = nds(nMC, stringBR)
    }    
    
  }
  
  return(qySample)
}
arrangeSample = function(S, useRanks = FALSE) {

  # Reorder MC samples of branching ratios for each wavelength 
  # in order to introduce serial correlation.
  # - for each sample at wavelength i, look for
  #   the *closest* sample at wavelength i+1.
  # - two distances possible: L1 based on ranks, and
  #   L2 based on BR values.
  # - sample at i+1 is reordered and used as reference for
  #   next step.      
     
  dd = dim(S)
  nMC = dd[1]; nWl = dd[2]; nBR = dd[3]
  pow = ifelse(useRanks, 1, 2)
  
  for(iWl in 1:(nWl-1)) {
    X = S[ , iWl,]
    Y = S[ , iWl+1,]
    if(useRanks) {
      for(iBR in 1:nBR) {
        X[ ,iBR] = rank(S[ ,iWl  ,iBR])
        Y[ ,iBR] = rank(S[ ,iWl+1,iBR])
      }
    }
    iord = vector("integer",nMC)
    for (iMC in 1:(nMC - 1)) {
      Z = abs(t(Y) - X[iMC,])^pow
      closest = which.min(colSums(Z))
      iord[iMC] = closest
      Y[closest,] = rep(nMC^3, nBR) # Replace by very large value
    }
    iord[nMC]    = setdiff(1:nMC, iord)
    S[ , iWl+1,] = S[iord, iWl+1,]
  }
  
  return(S)
}
gm = function(X) { 
  # Geometric mean
  x = X[X != 0]
  if(length(x) == 0) 
    return(NA)
  else
    return( prod(x)^(1/length(x)) )
} 
fScaleDirg = function(d) {
  # Scale for b params of Dirg(a;b) to recover
  # prescribed mean relative uncert.
  # fd = c(
  #   c( NA, 0.99, 0.93, 0.92, 0.92, 0.93, 0.94, 0.95, 0.96, 0.96),
  #   rep(0.96,50)
  # )
  # fd[d] * ( 1 + 1/d )
  1
}
