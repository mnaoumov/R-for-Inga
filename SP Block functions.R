generateRandomMatrix = function(randomGenerator, b, k) {
  A = matrix(randomGenerator(b * k), b, k)
  A = A - apply(A, 1, mean)
  A = A / sqrt(mean(A ^ 2))
}