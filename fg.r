# FG-algorithm
# Given in Flury and Gautschi (1986).
# Simultaneously transforms k positive definite matrices of dimension p x p to near diagonal form.

library(plyr)
library(gregmisc) # For combinations()

fg <- function(mat.list, p, alpha = 0, epsilon.f = 0.0001, epsilon.g = 0.0001) {
	# TODO: Ensure that the matrices are each p x p before commencing.
	B <- diag(p)
	f_algorithm(B)
}

f_algorithm <- function(A, B, p) {
	all.pairs <- combinations(p, 2)
	
	aaply(all.pairs, 1, function(pair) {
		l <- pair[1]
		j <- pair[2]
		# TODO: Compute T_i for i, j and l
		# TODO: Perform the G-algorithm on (T_1, ..., T_k) to get Q
		# Hstar <- H %*% Q
	})
}

# T is a list of k matrices.
g_algorithm <- function(T.list, epsilon.g = 0.0001) {
	Q <- diag(2)
	# TODO: Compute delta_{ij}'s
	# TODO: Compute T
	# TODO: Compute the (normalized) eigenvectors of T.
	# Q = (q1, q2)
	#	q1 <- first eigenvector of T
	#	q2 <- second eigenvector of T
	# Compare the 2-norm of Q_{g-1} - Q
}

# A measure of "deviation from diagonality" of a positive definite symmetrix matrix x.
diagonal.dev <- function(x) {
	det(diag(x)) / det(x)
}