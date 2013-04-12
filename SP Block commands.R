# should be changed if files are placed into a different folder
setwd("c:/Inga/")

source('SP Block functions.R')

# make predictable 'random' results
# 12345 could be changed to different value to produce other 'random' results
# comment out next line if you need truly random value
set.seed(12345)

b = 7
k = 4
A = generateRandomMatrix(rexp, b, k)
iterations = 10
u = (6 : 15) / 10
M = 10

xxx = NULL; for (i in 1 : iterations) xxx = c(xxx, apply(t(apply(A, 1, sample, size = k - 1)), 2, mean))
xxx = matrix(xxx, k - 1, iterations)

Lstat = NULL; for (i in 1 : iterations) Lstat = c(Lstat, Lambda(A, xxx[,i])$max)
Fstat = b * apply(xxx ^ 2, 2, sum)

apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
round(1 - pchisq(b * u ^ 2, k - 1), 4)

tplr = NULL
for (up in u) tplr = c(tplr, TailDistrMonteCarlo(A, up, M)$LugRice)
round(tplr, 4)

tpnc = NULL
for (up in u) tpnc = c(tpnc, TailDistrMonteCarlo(A, up, M)$NielCox)
round(tpnc, 4)