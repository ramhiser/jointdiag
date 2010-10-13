Q <- c(2,1) * diag(2)

T1 <- matrix(c(3,1,1,4), nrow = 2, byrow = TRUE)
T2 <- matrix(c(1, -1, -1, 4), nrow = 2, byrow = TRUE)
T3 <- matrix(c(3, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)

T.list <- list(T1, T2, T3)

# On scratch paper, I showed that the resulting set of deltas should be the 3x2 matrix
#	12, 4
#	4, 4
#	12, 1
#
# In other words:
actual.delta <- matrix(c(12, 4, 4, 4, 12, 1), ncol = 2, byrow = TRUE)

delta <- compute.deltas(Q, T.list)

cat("Our automated calculation of delta is equal to what we computed on scratch paper?", all(actual.delta == delta), "\n")

# On scratch paper, using weights of 1 (n_i = 1, for all i) I computed T to be the 2 x 2 matrix
actual.T <- matrix(c(39 / 12, 5 / 8, 5 / 8, 19 / 12), ncol = 2, byrow = TRUE)

T <- aaply(seq_len(k), 1, function(i) {
	weights[i] * (delta[i, 1] - delta[i, 2]) / (delta[i, 1] * delta[i, 2]) * T.list[[i]]
})
T <- colSums(T)

cat("Our automated calculation of T is equal to what we computed on scratch paper?", all(actual.T == T), "\n")


Q <- g_algorithm(T.list, verbose = TRUE)$Q