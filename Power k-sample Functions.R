library(maxLik)

generateVector = function(randomGenerator, N) {
  a = randomGenerator(N)
  N = length(a)
  a = (a - mean(a)) / sqrt(var(a) * (N - 1) / N)
  a
}

kcgf = function(tau, xd, f) {
  #This function calculates K(tau), K"(tau) and K"(tau), where
  #t is a 2*(k-1) dimensional matrix, K(t) is a real number.
  #Inputs:t  in K(t).      
  #       xd=the population: a=1,a=2, ..., a=N.
  #       f=vector of ratios of group sizes (n_1,..  .,n_{k-1}) and population size (N), 
  #Outputs:cgf2$cgf:K(tau)
  #        cgf2$dcgf:K'(tau)
  #        cgf2$ddc gf:K"(tau)
  #           
  taum = matrix(tau, length(f), 2)
  nn = length(xd)
  xsd = sqrt((var(xd) * (nn - 1)) / nn)
  xdb = mean(xd)
  x = (xd - xdb) / xsd
  f1 = 1 - sum(f)
  nn1 = rep(1, nn)
  x2 = rbind(nn1, x)
  kt = mean(log(f1 + f %*% exp(taum %*% x2)))
  h = t(t(f * exp(taum %*% x2)) / c(f1 + f %*% exp(taum %*% x2)))
  hx = t(x * t(f * exp(taum %*% x2)) / c(f1 + f %*% exp(taum %*% x2)))
  hhx = rbind(h, hx)
  dk = c(h %*% t(x2) / nn)
  ddk = (cbind(rbind(diag(c(h %*% nn1)), diag(c(h %*% x))),
               rbind(diag(c(h %*% x)), diag(c(hx %*% x))))
         - hhx %*% t(hhx)) / nn
  list(cgf = kt, dcgf = dk, ddcgf = ddk)
}

kLam = function(x1, t0, xd, f) {
  #This function solves  K'(t)=x for a given x=(f,x1).   
  #Inputs: x1: a given k-1-dimension vector.
  #        t0: an initial value for Newton's method of finding the 
  #            root of  k'(t)-x for a given x.
  #     xd: population.
  #     f=n/N
  #Outputs: Lam$lam=\Lambda(x)
  #         Lam$t=t(x), 3-dimension vector, t(x)%*%x-K(t(x))=sup{t.x-K(t): t}
  #         Lam$$kdd K"(t(x))
  #Get t0=t0(x), t0(x)%*%x-K(t0(x))=sup{t%*%x-K(t): t}
  
  x = c(f, x1 * 0.999)
  for (i in 1 : 10) {   
    ktt = kcgf(t0, xd, f)
    
    #print(det(ktt$ddcgf))
    t0 = t0 + solve(ktt$ddcgf) %*% (x - ktt$dcgf)
  }
  ktt = kcgf(t0,xd,f)
  list(lam = -ktt$cgf + x %*% t0, t = t0, kdd = ktt$ddcgf, err = x - ktt$dcgf)
}

F = function(x, N) {
  N * sum(x ^ 2)
}

getMean = function(b, f) {
  k = length(f) + 1
  N = sum(n)
  f1 = c(f, 1 - sum(f))
  n = f * N
  
  means = numeric(k - 1)
  lastIndex = 0
  for (i in 1 : (k - 1)) {
    nextIndex = lastIndex + n[i]
    means[i] = sum(b[(lastIndex + 1) : nextIndex]) / N;
    lastIndex = nextIndex
  }
  
  means
}

PowerForMu = function(mu, distr, f, Z, M) {
  k = length(f) + 1
  N = length(mu)
    
  pFValues = numeric(Z)
  #pLValues = numeric(Z)
  pTValues = numeric(Z)
  
  for (i in 1 : Z) {
    a = generateVector(distr, N)
    a = a + mu
    xsd = sqrt((var(a) * (N - 1)) / N)
    xdb = mean(a)
    a = (a - xdb) / xsd
    Lbar = getMean(a, f)
    Fbar = c(Lbar, -sum(Lbar))
  print(Lbar)
    
    FVal = F(Fbar, N)
    LVal = c(kLam(Lbar, rep(0, 2 * (k - 1)), a, f)$lam)
    
    FValues = NULL
    #LValues = NULL
    
    xxx = NULL; for(j in 1 : M) xxx = c(xxx, getMean(sample(a), f))
    xxx = matrix(xxx, k - 1, M)
    FValues = NULL; for (j in 1 : M) FValues = c(FValues, F(c(xxx[,j], -sum(xxx[,j])), N))
    

    pFValues[i] = mean(FValues > FVal)
    #pLValues[i] = mean(TValues > TVal)
    pTValues[i] = kHu(a, f, LVal, 10)
  }
  
  PowerF = mean(pFValues < .05)
  #PowerL = mean(pLValues < .05)
  PowerT = mean(pTValues < .05)
  
  #list (PowerF = PowerF, PowerL = PowerL, PowerT = PowerT)
  list (PowerF = PowerF, PowerT = PowerT)
}


kHu = function(xd, f, v, M)
{
  #this function calculates SP tail probs for k-sample perm test based on Lam(x)
  #input xd data, f group proportions for k-1, v P(Lam(X)>v),M= no of MC in integral
  #Formula (12) used for delta. Note |V0|^.5 and |V_\tau_0| both equal  sqrt(f1...fk)
  nn = length(xd)
  n = f * nn
  xsd = sqrt((var(xd) * (nn - 1))/nn)
  xdb = mean(xd)
  x = (xd - xdb) / xsd
  
  f1 = 1 - sum(f)
  lf = length(f)
  tailp1 = 1 - pchisq(nn * v ^ 2, lf) #first order approximation
  VV = kcgf(rep(0, 2 * lf), xd, f)$ddcgf[(lf + 1) : (2 * lf),
                                         (lf + 1) : (2 * lf)]
  svdV = svd(VV)
  Vh = svdV$u %*% diag(svdV$d ^ (1 / 2)) %*% t(svdV$v)
  Gu = 0
  sdelus = matrix(0, M, lf + 3)
  
  Mn = 0
  for (i in 1 : M) {                             
    s = rnorm(lf)
    s = s / sqrt(sum(s ^ 2))
    
    s = c(Vh %*% s)
    r = v
    t0 = rep(0, 2 * lf)
    for (j in 1 : 20) {
      x1 = c(f, r * s)
      tem = function(t0) {x1 = c(f, r * s); kcgf(t0, xd, f)$cgf - t0 %*% x1 }
      klt = optim(t0, tem, control = list(maxit = 50), method = "Nelder-Mead")
      t0 = c(klt$par)
      r = r + c((v - sqrt(abs(2 * (kcgf(t0, xd, f)$cgf - t0 %*% x1)))) * 
                  sqrt(abs(2 * (kcgf(t0, xd, f)$cgf - t0 %*% x1))) /
                  sum(s * t0[(lf + 1) : (2 * lf)]))
    }
    klt = optim(t0, tem, control = list(maxit = 2000))
    t0 = c(klt$par)
    delus = f1 * prod(f) * r ^ (lf - 1) /
      (v ^ (lf - 2) * sum(s * klt$par[(lf + 1) : (2 * lf)]) *
         (det(kcgf(t0, xd, f)$ddcgf) ^ .5))  

    Gu = ifelse(is.na(delus), Gu, Gu + delus)
    Mn = ifelse(is.na(delus), Mn, Mn + 1)
    sdelus[i,] = c(delus, s, v - sqrt(-2 * klt$value),
                   sum(s * klt$par[(lf + 1) : (2 * lf)]))
  }
  cn = (nn ^ (lf / 2)) / (2 ^ (lf / 2 - 1) * gamma(lf / 2))      
  tailp = ifelse(M == 0, 0, tailp1 + nn ^ (-1) * cn * v ^ (lf - 2) *
                   exp(-nn * v ^ 2 / 2) * (Gu / Mn - 1)) #LR
    
  tailpBN = ifelse(M == 0, 0, 1 - pchisq(nn * (v - log(Gu / Mn) /
                                                 (nn * v)) ^ 2, lf))
  list(tailpchi2 = tailp1, tailp = c(tailp, tailpBN), sdelus = sdelus)
}

mcanov = function(xd, f, MC = 10000) {
#For data vector xd in groups of proportions f (f_k left out) gives MC for Lam and X^2 
  nn = length(xd)
  n = f * nn
  xsd = sqrt((var(xd) * (nn - 1)) / nn)
  xdb = mean(xd)
  x = (xd - xdb) / xsd
  f1 = 1 - sum(f)
  lf = length(f)
  N = c(cumsum(n), nn)
  xb = rep(0, lf)
  for(i in 1 : lf) xb[i] = sum(x[(N[i] + 1) : N[i + 1]]) / nn
  MCpv = rep(0, MC)
  MCsqpv = rep(0, MC)
  for (k in 1 : MC) {
    xs = sample(x, nn)
    xsb = rep(0, lf)
    for(i in 1 : lf) xsb[i] = sum(xs[(N[i] + 1) : N[i + 1]]) / nn
    MCpv[k] = kLam(xsb, rep(0, 2 * lf), xd, f)$lam
    xsbk = c(xsb, -sum(xsb))
    MCsqpv[k] = nn * xsbk %*% diag(1 / c(f, f1)) %*% xsbk
  }

  list(mcreps = cbind(MCpv, MCsqpv))
}
