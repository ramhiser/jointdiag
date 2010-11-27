jointDiag.cov <- function(df, shrink = FALSE, shrink.param = 0.01) {
	if(shrink == FALSE) {
		shrink.param <- 0
	}
	covs <- daply(df, .(labels), function(class.df) {
		n.k <- nrow(class.df)
		p <- ncol(class.df) - 1
		(1 - shrink.param) * (n.k - 1) * cov(class.df[,-1]) / n.k + shrink.param * diag(p)
	})
	covs <- aperm(covs, perm = c(2,3,1))

	dimnames(covs) <- NULL
	covs
}

jointDiag.transform <- function(df, B) {
	transformed.df <- ddply(df, .(labels), function(class.df) {
		x <- as.matrix(class.df[,-1])
		dimnames(x) <- NULL
		data.frame(x %*% t(B))
	})
	transformed.df
}