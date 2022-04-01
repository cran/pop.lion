\name{project}
\alias{project}
\alias{C_montecarlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lion population projections
}
\description{
Run stochastic lion population projections.
}
\usage{
project(years,
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
     details)
}

\arguments{

  \item{years}{
A number: number of years to simulate the population.
}

  \item{runs}{
A number: number of times (or Monte Carlo runs) to simulate the population.
}

  \item{survival}{
A matrix: average monthly survival for each sex.
}

  \item{litter_distribution}{
A vector: probability distribution of litter sizes (1-5 cubs) in the population.
}

  \item{pop_initial}{
A number: number of prides (and coalitions). A simulation starts with an equal number of prides and coalitions.
}

  \item{conflict_age}{
A vector: the minimum age in months at which lions can be killed by conflict for females and males.
}

  \item{conflict_mortality}{
An array: mortality added at the edge by conflict for every month of the simulation and for females and males. Expressed in percentage, a value of 15.2 will be understood by the model as 15.2 per cent. Values can be double. The array has 12 * years rows.
}

  \item{hunting_age}{
A vector: the minimum age in months at which lions can be killed by trophy hunting for females and males.
}

  \item{hunting_mortality}{
An array: mortality added at the edge by trophy hunting for every month of the simulation and for females and males. Expressed in number of individuals, a value of 15 will be understood by the model as 15 killed every month. A value of 0.5 will be understood as 6 lions killed per year. The array has 12 * years rows.
}

\item{hunter_error}{
A number: hunter error.
}

  \item{K_indiv}{
A number: maximum number of individuals in the population.
}

  \item{K_pride}{
A number: maximum number of prides in the population.
}

  \item{K_coali}{
A number: maximum number of coalitions in the population.
}

  \item{K_edged}{
A number: number of prides in the population that are located at the edge of the reserve and therefore vulnerabe to hunting and poaching.
}

  \item{seed}{
(optional) A number: seed of the random number generator.
}

  \item{details}{
(optional) A boolean: indicate whether individual events are exported. This can generate large simulation objects.
}


}
\details{
Run stochastic lion population projections with an Individual-Based Model (IBM) compiled in C.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{runs}{a 3-dimensional array of numbers of individuals with dimension c(years, statistics, runs)}

  \item{individuals}{a 2-dimensional array of individuals events}

  \item{parameters}{a list of parameters of the projection}
%% ...
}

\examples{

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

# Population size at the end of the simulation:
apply(projection$runs[,"NINDIV",], 1, mean)[12*years+1]

}
