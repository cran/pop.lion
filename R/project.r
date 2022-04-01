#########################################################################
# MONTECARLO SINGLE CORE
#########################################################################

project <- function(
	years,
	runs,
	survival,
	litter_distribution,
	pop_initial,
	conflict_age,
	conflict_mortality,
	hunting_age,
	hunting_mortality,
	hunter_error,
	K_indiv,
	K_pride,
	K_coali,
	K_edged,
	seed,
	details
) {

	if (missing(litter_distribution)) {
		litter_distribution <- c(0.12,0.30,0.35,0.19,0.04)
	}

	if (missing(pop_initial)) {
		pop_initial <- 5
	}

	if (missing(conflict_age)) {
		conflict_age <- array(120, dim=c(2), dimnames=list(c("female", "male")))
	}

	if (missing(conflict_mortality)) {
		conflict_mortality <- array(0, dim=c(12*years+1, 2), dimnames=list(NULL, c("female", "male")))
	}

	if (missing(hunting_age)) {
		hunting_age <- array(120, dim=c(2), dimnames=list(c("female", "male")))
	}

	if (missing(hunting_mortality)) {
		hunting_mortality <- array(0, dim=c(12*years+1, 2), dimnames=list(NULL, c("female", "male")))
	}

	conflict_mortality2 <- array(0, dim=c(12*years, 2))
	j <- 0
	for (sex in 1:2) {
		j <- j+1
		conflict_mortality2[,j] <- conflict_mortality[, sex]
	}

	hunting_mortality2 <- array(0, dim=c(12*years, 2))
	j <- 0
	for (sex in 1:2) {
		j <- j+1
		hunting_mortality2[,j] <- hunting_mortality[, sex]
	}

	if (missing(hunter_error)) {
		hunter_error <- 0
	}

	if (missing(K_edged)) {
		K_edged <- 0
	}

	if (K_edged > K_pride) {
		K_edged = K_pride
	}

	if (missing(seed)) {
		seed <- 1
	}
	set.seed(seed)

	if (missing(details)) {
		details <- F
	}

	months <- (years*(12)+1)
	stats <- 35

	output <- .Call("C_montecarlo",
									as.integer(years),
									as.integer(runs),
									as.double(as.vector(data.matrix(survival))),
									as.double(as.vector(litter_distribution)),
									as.integer(as.vector(data.matrix(conflict_age))),
									as.integer(as.vector(data.matrix(hunting_age))),
									as.double(as.vector(data.matrix(conflict_mortality2))),
									as.double(as.vector(data.matrix(hunting_mortality2))),
									as.integer(hunter_error),
									as.integer(pop_initial),
									as.integer(K_indiv),
									as.integer(K_pride),
									as.integer(K_coali),
									as.integer(K_edged)
	)

	results <- list()
	results$runs <- aperm(`dim<-`(t(output$runs), c(stats, months, runs)), c(2, 1, 3))

	if(details) {
		ni <- matrix(output$individuals, nrow=length(output$individuals)/months, ncol=months, byrow=T)
		results$individuals <- ni[nrow(ni):1,]
	}

	colnames(results$runs) <- c(
		"NINDIV",
		"NPRIDES",
		"NCOALIS",
		"NCOALIS_RESIDENT",
		"NCOALIS_VAGRANT",
		"NPRIDES_RESIDENT",
		"NPRIDES_VAGRANT",
		"COALISIZE_RESIDENT",
		"COALISIZE_VAGRANT",
		"PRIDESIZE_RESIDENT",
		"PRIDESIZE_VAGRANT",
		"NFEMALES",
		"NMALES",
		"TAKEOVERS",
		"LITTERS",
		"AGE",
		"NMALES_ADULT",
		"NFEMALES_ADULT",
		"NFEMALES_DISPERSED",
		"NPRIDES_EDGED",
		"NPRIDES_CORE",
		"PRIDESIZE_CORE",
		"PRIDESIZE_EDGED",
		"NFEMALES_CORE",
		"NCUBS_CORE",
		"NFEMALES_EDGE",
		"NCUBS_EDGE",
		"COALITION_TENURE",
		"NINDIV_CORE",
		"NINDIV_EDGE",
		"FAILED_HUNTS",
		"NPRIDES_PER_COALI",
		"IBI",
		"NMALES_EDGE",
		"NMALES_CORE"
	)

	results$parameters <- list(	years = years,
															runs = runs,
															survival = survival,
															litter_distribution = litter_distribution,
															conflict_age = conflict_age,
															conflict_mortality = conflict_mortality,
															hunting_age = hunting_age,
															hunting_mortality = hunting_mortality,
															hunter_error = hunter_error,
															pop_initial = pop_initial,
															K_indiv = K_indiv,
															K_pride = K_pride,
															K_coali = K_coali,
															K_edged = K_edged,
															seed = seed
	)

	x <- gc(verbose = FALSE)

	return(results)

}
