kooij = function(p, tempRange, M, T0=300, T1=300) {
  p[1] * (tempRange / T0)^p[2] * exp(-p[3]/tempRange) *
    p[4] * exp(p[5] * abs(1/tempRange -1/ T1))
}
k3body = function(p, tempRange, M, T0 = 300) {
  # kAss from Dobrijevic2016
  k0   = kooij(p[1:5], tempRange, M, T0)
  kInf = kooij(p[6:10], tempRange, M, T0)
  kr   = kooij(p[11:15], tempRange, M, T0)
  Fc   = p[16]
  Ni   = 1
  
  F1 = exp(log(Fc) / (1 + (log(k0 * M + kInf) / Ni) ^ 2))
  k  = kInf * (k0 * M * F1 + kr) / (k0 * M + kInf)
  
  return(k)
}
kEq18VV = function(p, tempRange, M, T0 = 1) {
  # kAss from Vuitton2019
  # Note: kInf and k0 in reverse order from k3body
  kInf = kooij(p[1:5], tempRange, M, T0)
  k0   = kooij(p[6:10], tempRange, M, T0)
  kR   = kooij(p[11:15], tempRange, M, T0)
  Fc   = p[16]
  
  if (max(kR) > min(0.99 * kInf) ) {
    k = kInf
  } else {
    C = -0.4 - 0.67 * log10(Fc)
    N = 0.75 - 1.27 * log10(Fc)
    lPr = log10(k0 * M / kInf)
    fExp = 1 + ((lPr + C) / (N - 0.14 * (lPr + C))) ^ 2
    lF1 = log10(Fc) / fExp
    kInf1 = kInf - kR
    k = kR + 10 ^ lF1 * kInf1 * k0 * M / (kInf1 + k0 * M)
  }
  
  return(k)
}
# # Troe formula
# Pr   = k0 * M / kInf
# cExp = -0.4 - 0.67 * log10(fc)
# NExp = 0.75 - 1.27 * log10(fc)
# dExp = 0.14
# fExp = 1 + ((log10(Pr) + cExp) / (NExp - dExp * (log10(Pr) + cExp)))^2
# broadF = fc ^ (1 / fExp)
# k = kInf * (Pr / (1 + Pr)) * broadF
