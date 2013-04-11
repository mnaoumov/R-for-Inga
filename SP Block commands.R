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