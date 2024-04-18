k_kooij = function(p, tempRange, T0=300, T1=300) {
  p[1] * (tempRange / T0)^p[2] * exp(-p[3]/tempRange) *
    p[4] * exp(p[5] * abs(1/tempRange -1/ T1))
}
k_assoc0 = function(p, tempRange, M, T0 = 300) {
  k0M  = k_kooij(p[1:5], tempRange, T0) * M
  return(k0M)
}
k_assocmd = function(p, tempRange, M, T0 = 300) {
  # kAss from Dobrijevic2016, corrected...
  k0M  = k_kooij(p[1:5],   tempRange, T0) * M
  kInf = k_kooij(p[6:10],  tempRange, T0)
  kR   = k_kooij(p[11:15], tempRange, T0)
  Fc   = p[16]
  
  N   = 0.75 - 1.27 * log10(Fc)
  lPr = log10(k0M / kInf)
  fExp = 1 + (lPr / N)^2
  lF1  = log10(Fc) / fExp
  k  = kInf * (kR + k0M * 10^lF1) / (kInf + k0M)
  
  return(k)
}
k_assocvv = function(p, tempRange, M, T0 = 1) {
  # kAss from Vuitton2019
  # Note: kInf and k0 in reverse order from k_assocMD and T0 = 1
  kInf = k_kooij(p[1:5],   tempRange, T0)
  k0M  = k_kooij(p[6:10],  tempRange, T0) * M
  kR   = k_kooij(p[11:15], tempRange, T0)
  Fc   = p[16]
  
  if (max(kR) > min(0.99 * kInf) ) {
    k = kInf
  } else {
    C = -0.4 - 0.67 * log10(Fc)
    N = 0.75 - 1.27 * log10(Fc)
    lPr = log10(k0M / kInf)
    fExp = 1 + ((lPr + C) / (N - 0.14 * (lPr + C)))^2
    lF1 = log10(Fc) / fExp
    kInf1 = kInf - kR
    k = kR + (kInf1 * k0M * 10^lF1) / (kInf1 + k0M)
  }
  
  return(k)
}

k_assoctroe = function(p, tempRange, M, T0 = 300) {
  fTroe = function(k0, kInf, P) {
    # Simplified Troe formula
    Fc = 0.64
    lFc= log10(Fc)
    Pr = k0 * P / kInf
    cE = -0.4 -0.67*lFc
    nE = 0.75 -1.27*lFc
    dE = 0.14
    c1 = log10(Pr) + cE
    fE = 1.0 + (c1/(nE-dE*c1))**2
    return( kInf * (Pr/(1.0+Pr)) * Fc**(1.0/fE) )  
  }
  k0   = k_kooij(p[1:5],  tempRange, T0)
  kInf = k_kooij(p[6:10], tempRange, T0)
  return( fTroe(k0, kInf, M) )
}
