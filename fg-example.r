# The two matrices A1 and A2 are given in the Flury and Gautschi (1986) paper.
# We run our implementation of the FG-algorithm to ensure that we obtain the correct results.
A1 <- matrix(
		c(45, 10, 0, 5, 0, 0,
		10, 45, 5, 0, 0, 0,
		0, 5, 45, 10, 0, 0,
		5, 0, 10, 45, 0, 0,
		0, 0, 0, 0, 16.4, -4.8,
		0, 0, 0, 0, -4.8, 13.6),
	nrow = 6, ncol = 6, byrow = TRUE)
	
A2 <- matrix(
		c(27.5, -12.5, -0.5, -4.5, -2.04, 3.72,
		-12.5, 27.5, -4.5, -0.5, 2.04, -3.72,
		-0.5, -4.5, 24.5, -9.5, -3.72, -2.04,
		-4.5, -0.5, -9.5, 24.5, 3.72, 2.04,
		-2.04, 2.04, -3.72, 3.72, 54.76, -4.68,
		3.72, -3.72, -2.04, 2.04, -4.68, 51.24),
	nrow = 6, ncol = 6, byrow = TRUE)
	
A <- list(A1 = A1, A2 = A2)

source("fg.r")
B <- fg(A, p = 6, verbose = TRUE)