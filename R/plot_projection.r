#########################################################################
# PLOT POPULATION COUNTS
#########################################################################

plot_projection <- function(projection, title) {

	# if (title == "Packs") varname <- "packs"
	# if (title == "Pairs") varname <- "pairs"
	# if (title == "Reproductions") varname <- "repros"
	# if (title == "Population") varname <- "pop_size"
	
	varname <- title

	years <- projection$parameters$years
	runs <- projection$runs

	dmed <- apply(runs[,varname,], 1, mean)
	dlci <- apply(runs[,varname,], 1, quantile, 0.025)
	duci <- apply(runs[,varname,], 1, quantile, 0.975)
	dmonths <- 1:nrow(runs)
	dyears <- c(1, (1:projection$parameters$years)*12 + 1)

	plot(0,0, type="l", xlab="Months", ylab="", xlim=c(1,nrow(runs)), ylim=c(0, max(duci)), main=title)
	polygon(x=c(dmonths, rev(dmonths)), y=c(dlci, rev(duci)), col="gray95", border=NA)
	lines(dmonths, dlci, lty="dotted", col="darkgrey")
	lines(dmonths, duci, lty="dotted", col="darkgrey")
	lines(dmonths, dmed, lty="solid")

}

