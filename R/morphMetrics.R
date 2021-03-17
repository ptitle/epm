pullNeighbor <- function(Row, Names) {
	Order <- order(Row)
	Vals <- Row[Order]
	Labels <- Names[Order]
	NN_dist <- Vals[2]
	NN_id <- Labels[2]
	NN_mean <- mean(Vals[2:length(Vals)])
	NN_sd <- ifelse(length(Labels) > 2, sd(Vals[2:length(Vals)]), 0)
	NN_max <- max(Vals[2:length(Vals)])
	Output <- data.frame(di = NN_dist, NN_id = NN_id, mean_dist = NN_mean, sd_dist = NN_sd, max_dist = NN_max)
	return(Output)
}


genRi <- function(Row, Mins, Maxs, Ntax, Nrep = 20) {
	Sims <- replicate(Nrep, rbind(Row, sapply(1:length(Mins), function(x) runif(Ntax - 1, min = Mins[x], max = Maxs[x]))))
	Dists <- sapply(1:dim(Sims)[3], function(sim) pullNeighbor(as.matrix(dist(Sims[,,sim]))[1,], Names = 1:dim(Sims)[2]))
	Ri <- unlist(Dists["di",])
	return(mean(Ri))
}


nnDist <- function(Mat, dMat = NULL, Nrep = 20) {
	if (all(is.na(Mat))) {
		return(NA)
	}

	if (is.vector(Mat) | nrow(Mat) == 1) {
		# This matrix only has one species, so disparity is 0.
		outMat <- data.frame(pi = 0, ri = 0, di = 0, NN_id = 0, mean_dist = 0, sd_dist = 0, max_dist = 0)
	} else {
		Means <- apply(Mat, 2, mean)
		SDs <- apply(Mat, 2, sd)
		Mins <- Means - sqrt(3) * SDs # eq 6 in Foote 90
		Maxs <- Means + sqrt(3) * SDs # eqn 7 in Foote 90
		Ntax <- nrow(Mat)
		if (is.null(dMat)) {
			dMat <- dist(Mat) # edit this line to use whatever distance function you want. Also edit line in genRi if not using Euclidean distance
		}
		dMat <- as.matrix(dMat)
		outList <- apply(dMat, 1, pullNeighbor, Names = colnames(dMat))
		Ris <- apply(Mat, 1, genRi, Mins = Mins, Maxs = Maxs, Ntax = Ntax, Nrep = Nrep)
		outMat <- do.call(rbind, outList)
		outMat <- data.frame(pi = (outMat$di - Ris) / Ris, ri = Ris, outMat)
	}
	return(outMat)
}