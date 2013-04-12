# used for permutations
library(gtools)

# used for maxNR
library(maxLik)

# used for caching results for the same values
library(memoise)

generateRandomMatrix = function(randomGenerator, b, k) {
  A = matrix(randomGenerator(b * k), b, k)
  A = A - apply(A, 1, mean)
  A = A / sqrt(mean(A ^ 2))
}

getb = function(A) dim(A)[1]
getk = function(A) dim(A)[2]

perm = memoise(function(A, i) {
  k = getk(A)
  permutations(k, k - 1, A[i,])
})

expb = memoise(function(A, i, beta) exp(perm(A, i) %*% beta))

kappa = memoise(function(A, beta) {
  b = getb(A)
  k = getk(A)
  s = 0
  for (i in 1 : b) s = s + log(sum(expb(A, i, beta)) / prod(1 : k))
  s / b
})
  
kappad = memoise(function(A, beta) {
  b = getb(A)
  k = getk(A)
  s = 0
  for (i in 1 : b) {
    expb = expb(A, i, beta)
    s = s + t(perm(A, i)) %*%  expb / sum(expb)
  }
  
  s / b
})

kappad2 = memoise(function(A, beta) {
  b = getb(A)
  k = getk(A)
  s = 0
  for (i in 1 : b) {
    perm = perm(A, i)
    expb = expb(A, i, beta)
    s = s - (t(perm) %*% expb %*% t(expb) %*% perm) / sum(expb) ^ 2 +
      (t(perm) %*% diag(c(expb)) %*% perm) / sum(expb)
  }

  s / b
})

orderedMeans = memoise(function(A) {
  b = getb(A)
  k = getk(A)
  ASorted = matrix(numeric(), b, k)
  for(i in 1 : b)
    ASorted[i,] = sort(A[i,])
  
  colMeans(ASorted)
})


isPointInsideDomain = memoise(function(A, x, tolerance = 1) {
  k = getk(A)
  
  y = x
  if (length(y) == k - 1)
    y[length(y) + 1] = -sum(y)
  
  stopifnot(length(y) == k)
  
  for (usedIndices in 1 : (k - 1)) {
    subVectors = combinations(k, usedIndices, y, set = FALSE)
    subVectorSums = apply(subVectors, 1, sum)
    sumMins = sum(orderedMeans(A)[1 : usedIndices])
    sumMaxs = sum(orderedMeans(A)[(k - usedIndices + 1) : k])
    if (!all(subVectorSums > sumMins * tolerance))
      return(FALSE);
    if (!all(subVectorSums < sumMaxs * tolerance))
      return(FALSE);
  }
  return(TRUE)
})


putInsideDomain = memoise(function(A, x, tolerance) {
  if (isPointInsideDomain(A, x, tolerance))
  {
    return (x)
  }
  
  norm = sqrt(sum(x ^ 2))
  
  inside = 0
  outside = norm
  
  # 2^(-15) < 0.0001
  for (i in 1 : 15) {
    r = (inside + outside) / 2
    y = x / norm * r
    if (!isPointInsideDomain(A, y, tolerance))
      outside = r
    else
      inside = r
  }
  
  return (x / norm * inside)
})

Lambda = memoise(function(A, x, tolerance = 0.9999) {
  k = getk(A)
  x = putInsideDomain(A, x, tolerance)
  
  ll = function(beta) { beta %*% x - kappa(A, beta) }
  lll = function(beta) { x - kappad(A, beta)[,1] }
  llll = function(beta) { -kappad2(A, beta) }
  lam = try(maxNR(ll, grad = lll, hess = llll, start = rep(0, k - 1), tol = 1 - tolerance), silent = TRUE)
  
  if ((length(lam) == 1) && (class(lam) == "try-error")) {
    lam = Lambda(A, x * tolerance, tolerance)
  }
  
  lam
})

Sqr = memoise(function(M) { svd(M)$u %*% diag(sqrt(svd(M)$d)) %*% t(svd(M)$v) })

UniformOnSphere = function(n) {
  x = rnorm(n)
  r = sqrt(sum(x ^ 2))
  x / r
}

IntOnSphereMonteCarlo = function(h, d, M) {
  r = 0
  for(i in 1 : M)
    r = r + h(UniformOnSphere(d))
  r / M
}

V0h = memoise(function(A) {
  k = getk(A)
  Sqr(kappad2(A, rep(0, k - 1)))  
})

transformIntoSphere = function(A, r, s) {
  (r * V0h(A) %*% s)[,1]
}

delta = memoise(function(A, u, s, deltaIterations = 5) {
  k = getk(A)
  r = u
  V0h = V0h(A)
  
  for (i in 1 : iterations) {
    tr = transformIntoSphere(A, r, s)
    L = Lambda(A, tr)
    r = r - ((L$max - u ^ 2 / 2) / (L$est %*% V0h %*% s))[1,1]
  }
  
  betah = L$est
  
  (det(kappad2(A, betah)) ^ (- 1 / 2) * r ^ (k - 2) * det(V0h)) / (u ^ (k - 3) * abs(t(s) %*% V0h %*% betah))
})

TailDistrMonteCarlo = function(A, u, M, deltaIterations = 5) {
  k = getk(A)
  b = getb(A)
  
  delta2 = function(s) delta(A, u, s, deltaIterations)
  
  cb = (b ^ ((k - 1) / 2)) / (2 ^ ((k - 1) / 2 - 1) * gamma((k - 1) / 2))
  Gb = IntOnSphereMonteCarlo(delta2, k - 1, M)
  
  LugRice = 1 - pchisq(b * u ^ 2, k - 1) + cb / b * u ^ (k - 3) * exp(-b * u ^ 2 / 2) * (Gb - 1)
  ustar = u - log(Gb) / (b * u)
  NielCox = 1 - pchisq(b * ustar ^ 2, k - 1)
  list(LugRice = LugRice, NielCox = NielCox)
}
