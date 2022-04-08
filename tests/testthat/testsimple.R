library('pop.lion')

years = 25

survival <-  matrix(1, nrow=180, ncol=2)
survival[1:12, 1:2] <- 0.97^(1/12)
survival[13:24, 1:2] <- 0.98^(1/12)
survival[25:96, 1:2] <- 0.99^(1/12)
survival[97:108, 1:2] <- 0.98^(1/12)
survival[109:120, 1:2] <- 0.96^(1/12)
survival[121:132, 1:2] <- 0.94^(1/12)
survival[133:144, 1:2] <- 0.92^(1/12)
survival[145:156, 1:2] <- 0.90^(1/12)
survival[157:168, 1:2] <- 0.87^(1/12)
survival[169:180, 1:2] <- 0.83^(1/12)

litter_distribution <- c(0.10, 0.30, 0.35, 0.20, 0.05)

conflict_age <- array(4*12, dim=c(2), dimnames=list(c("female", "male")))
conflict_mortality <- array(0, dim=c(12*years, 2), dimnames=list(NULL, c("female", "male")))
conflict_mortality[24:36,] <- 15.2

hunting_age <- array(5*12, dim=c(2), dimnames=list(c("female", "male")))
hunting_mortality <- array(0, dim=c(12*years, 2), dimnames=list(NULL, c("female", "male")))
hunting_mortality[72:84,"male"] <- 10

projection <- project(
	years = years,
	runs = 100,
	survival = survival,
	litter_distribution = litter_distribution,
	pop_initial = 5,
	conflict_age = conflict_age,
	conflict_mortality = conflict_mortality,
	hunting_age = hunting_age,
	hunting_mortality = hunting_mortality,
	hunter_error = 0,
	K_indiv = 400,
	K_pride = 20,
	K_coali = 20,
	K_edged = 10,
	seed = 1,
	details = FALSE
)

oldpar <- par(mfrow = c(1,1))
par(mfrow=c(2,2))
plot_projection(projection, "NINDIV")
plot_projection(projection, "NPRIDES")
plot_projection(projection, "NCOALIS")
plot_projection(projection, "LITTERS")
par(oldpar)

# Population size at the end of the simulation:
apply(projection$runs[,"NINDIV",], 1, mean)[12*years+1]
