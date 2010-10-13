# FG-algorithm
# Given in Flury and Gautschi (1986).
# Simultaneously transforms k positive definite matrices of dimension p x p to near diagonal form.

library(plyr)
library(gregmisc) # For combinations()

# Performs the FG-algorithm as given in the Flury and Gautschi (1986) paper.
#	mat.list is a list of k positive definite symmetric matrices of dimension p x p.
fg <- function(mat.list, initial.B = diag(p), p, epsilon.f = 0.0001, epsilon.g = 0.0001, verbose = FALSE) {
	# TODO: Ensure that the matrices are each p x p before commencing.
	f.alg.out <- f_algorithm(A.list = mat.list, B = initial.B, p = p, verbose = verbose)
	f.alg.out
}

# Performs the F-algorithm as given in the Flury and Gautschi (1986) paper.
# Each iteration of the F-algorithm requires the G-algorithm.
f_algorithm <- function(A.list, B = diag(p), p, epsilon.f = 0.0001, max.iterations = 1000, norm = 2, verbose = FALSE) {
	k <- length(T.list)
	# Step F_0 and F_1: Initial approximation to the orthogonal matrix minimizing the "deviation from diagonality" metric.
	#	By default, the initial approximation is the p x p identity matrix.
	current.B <- B
	previous.phi <- prod(laply(A.list, function(A.i) {
		#diagonal.dev(t(current.B) %*% A.i %*% current.B)
		diagonal.dev(A.i, B)
	}))
	
	is.done <- FALSE
	current.iteration <- 1
	
	# We perform the F-algorithm until convergence or until we reach the maximum number of iterations,
	# which are specified by the user (default is 1000).
	while(!is.done) {
		# We generate a matrix of all pairs (l, j), 1 <= l < j <= p.
		all.pairs <- combinations(p, 2)

		# Step F_2: Repeat steps F_21 to F_24 for all pairs (l, j), 1 <= l < j <= p.
		a_ply(all.pairs, 1, function(pair) {
			l <- pair[1]
			j <- pair[2]
			b.l <- B[, l]
			b.j <- B[, j]
			
			# Step F_21: Compute T_i (i = 1, ..., k) and H.
			#	The resulting T_i are positive definite symmetric.
			H <- cbind(b.l, b.j)	
			T.list <- compute.T.list(A.list, b.l, b.j)
			
			# Step F_22: Perform the G-algorithm on (T_1, ..., T_k) to get an orthogonal 2 x 2 matrix Q.
			#	The resulting T_i are positive definite symmetric.
			Q <- g_algorithm(T.list = T.list, norm = norm)$Q
			
			# Step F_23: Orthogonal rotation of the two columns of H by an angle alpha.
			#	Hstar is matrix of dimension p x 2.
			#	The first column of Hstar is denoted b.l.star.
			#	The second column of Hstar is denoted b.j.star.
			Hstar <- H %*% Q
			
			# Step F_24: In the matrix B, we replace:
			#	the column b.l with b.l.star (the first column of Hstar),
			#	the column b.j with b.j.star (the second column of Hstar).
			B[, l] <<- Hstar[, 1]
			B[, j] <<- Hstar[, 2]
		})
		
		phi <- prod(laply(A.list, function(A.i) {
			#diagonal.dev(t(current.B) %*% A.i %*% current.B)
			diagonal.dev(A.i, B)
		}))
		
		# We calculate the "deviation from diagonality" for the matrix B from the current step and the previous step.
		# Then we look at the difference of the two.
		diff.dev.diag <- previous.phi - phi
		
		if(verbose) {
			cat("B\n")
			print(B)
			cat("\nIteration", current.iteration, "of F-algorithm with current deviation measure", previous.phi, "\n")
			cat("Difference between deviation measures:", diff.dev.diag, "\n\n")
		}
		
		# If the difference of the "deviations from diagonality" is less than epsilon, we have convergence and stop.
		if(diff.dev.diag < epsilon.f) {
			is.done <- TRUE
		}
		# If we have exceeded the maximum number of iterations specified by the user (default is 1000),
		# we stop the algorithm and warn the user.
		else if(current.iteration >= max.iterations) {
			is.done <- TRUE
			warning(paste("F-algorithm did not converge after", max.iterations, "iterations\n"))
		}
		# If we have not yet converged, we proceed to the next step of the algorithm and use
		# the value of B as our updated B_i.  Also, we increment the current iteration.
		else {
			current.B <- B
			previous.phi <- phi
			current.iteration <- current.iteration + 1
		}
	}
	if(verbose) {
		cat("F-algorithm stopped after", current.iteration, "iterations.\n")
	}
	list(B = B, num.iterations = current.iteration)
}

# Performs the G-algorithm as given in the Flury and Gautschi (1986) paper.
# T is a list of k positive definite symmetric 2 x 2 matrices.
# Q is an orthogonal 2 x 2 matrix.
g_algorithm <- function(T.list, Q = diag(2), weights = rep(1, k), epsilon.g = 0.0001, max.iterations = 1000, norm = 2, verbose = FALSE) {
	k <- length(T.list)
	current.Q <- Q
	
	is.done <- FALSE
	current.iteration <- 1
	while(!is.done) {
		# Step G_1: Compute the deltas
		delta <- compute.deltas(Q, T.list)
		
		# Step G_2: Compute T, which is a linear combination of the T_i's involving delta.
		#	See page 172.
		T <- aaply(seq_len(k), 1, function(i) {
			weights[i] * (delta[i, 1] - delta[i, 2]) / (delta[i, 1] * delta[i, 2]) * T.list[[i]]
		})
		T <- colSums(T)
		
		# Step G_3: Compute the (normalized) eigenvectors of T.
		#	The first eigenvector of T becomes the first column of Q.
		#	The second eigenvector of T becomes the second column of Q.
		Q <- eigen(T)$vectors
		
		# We compute the matrix norm of the difference of the Q from the current step and the previous step.
		matrix.norm <- pnorm(current.Q - Q, norm)
		
		if(verbose) {
			cat("Q\n")
			print(Q)
			cat("\nIteration", current.iteration, "of G-algorithm with norm", matrix.norm, "\n")
		}
		
		# If the norm of the difference is less than epsilon, we have convergence and stop.
		if(matrix.norm < epsilon.g) {
			is.done <- TRUE
		}
		# If we have exceeded the maximum number of iterations specified by the user (default is 1000),
		# we stop the algorithm and warn the user.
		else if(current.iteration >= max.iterations) {
			is.done <- TRUE
			warning(paste("G-algorithm did not converge after", max.iterations, "iterations\n"))
		}
		# If we have not yet converged, we proceed to the next step of the algorithm and use
		# the value of Q as our updated Q_i.  Also, we increment the current iteration.
		else {
			current.Q <- Q
			current.iteration <- current.iteration + 1
		}
	}
	if(verbose) {
		cat("G-algorithm stopped after", current.iteration, "iterations.\n")
	}
	list(Q = Q, num.iterations = current.iteration)
}

# Returns a (k x 2) matrix of the delta_{ij}'s as given in Equation (2.3) on p. 172 of paper.
# T.list is a list of k positive definite symmetric matrices of size 2 x 2.
# Q is an orthogonal 2 x 2 matrix.
compute.deltas <- function(Q, T.list) {
	laply(T.list, function(T.i) {
		# i is fixed.
		q1 <- Q[,1]
		q2 <- Q[,2]
		c(t(q1) %*% T.i %*% q1, t(q2) %*% T.i %*% q2)
	})
}

# Returns a list of k positive definite symmetric matrices of dimension 2 x 2.
# See Step F_22 of paper on page 171.
# A.list is a list of k positive definite symmetric matrices of size p x p.
# b.l and b.j are vectors of length p.
compute.T.list <- function(A.list, b.l, b.j) {
	llply(A.list, function(A.i) {
		T.i <- matrix(c(
			t(b.l) %*% A.i %*% b.l,
			t(b.l) %*% A.i %*% b.j,
			t(b.j) %*% A.i %*% b.l,
			t(b.j) %*% A.i %*% b.j
		), nrow = 2, byrow = 2)
		T.i
	})
}

# A measure of "deviation from diagonality" of a positive definite symmetrix matrix x.
diagonal.dev <- function(A.i, B) {
	prod(diag(t(B) %*% A.i %*% B)) / det(A.i)
}

# This computes the p-norm of a matrix, x.
# By default, we use the 2-norm (Euclidean norm).
pnorm <- function(x, p = 2) {
	sum(x^p)^(1/p)
}