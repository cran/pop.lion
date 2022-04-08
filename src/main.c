/* main.c
 *
 * Copyright (C) 2015-2022 Guillaume Chapron.
 * gchapron@carnivoreconservation.org
 * with contributions from Matthew Wijers, Andrew Loveridge and David Macdonald
 *
 * This file is part of 'pop.lion', a R package to simulate lion populations
 *
 * 'pop.lion' is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * 'pop.lion' is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with 'pop.lion'. If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "mc.h"
#include "pop.h"
#include "tools.h"

#define extern
# include "globals.h"
#undef extern

/*******************************************************************************
 Monte Carlo
 *******************************************************************************/

SEXP C_montecarlo (
                   SEXP SEXP_years,
                   SEXP SEXP_runs,

                   SEXP SEXP_surv,
                   SEXP SEXP_litter_distribution,

                   SEXP SEXP_conflict_age,
                   SEXP SEXP_hunting_age,
                   SEXP SEXP_conflict_mortality,
                   SEXP SEXP_hunting_mortality,
                   SEXP SEXP_hunter_error,

                   SEXP SEXP_initial_prides_coalitions,

                   SEXP SEXP_K_individuals,
                   SEXP SEXP_K_prides,
                   SEXP SEXP_K_coalitions,
                   SEXP SEXP_K_edged) {

    int nprot = 0;

    R_number_of_years = INTEGER(SEXP_years)[0];
    R_number_mc_runs = INTEGER(SEXP_runs)[0];

    number_of_months = 12*R_number_of_years;
    number_of_months++;

    R_survival_av = malloc(2 * sizeof(double*));
    for (int i = 0; i < 2; i++) {
        R_survival_av[i] = malloc(LENGTH(SEXP_surv)/2 * sizeof(double));
        for (int j = 0; j < LENGTH(SEXP_surv)/2; j++) {
            R_survival_av[i][j] = REAL(SEXP_surv)[j + i * LENGTH(SEXP_surv)/2];
        }
    }

    R_litter_distribution = malloc(LENGTH(SEXP_litter_distribution) * sizeof(double));
    for (int i = 0; i < LENGTH(SEXP_litter_distribution); i++) {
        R_litter_distribution[i] = REAL(SEXP_litter_distribution)[i];
    }

    R_hunter_error = INTEGER(SEXP_hunter_error)[0];

    R_conflict_age_female = INTEGER(SEXP_conflict_age)[0];
    R_conflict_age_male = INTEGER(SEXP_conflict_age)[1];
    R_hunting_age_female = INTEGER(SEXP_hunting_age)[0];
    R_hunting_age_male = INTEGER(SEXP_hunting_age)[1];

    R_mortality_cols = NUMBER_OF_SEXES;

    R_conflict_mortality = malloc(LENGTH(SEXP_conflict_mortality)/R_mortality_cols * sizeof(double*));
    for (int i = 0; i < LENGTH(SEXP_conflict_mortality)/R_mortality_cols; i++) {
        R_conflict_mortality[i] = malloc(R_mortality_cols * sizeof(double));
        for (int j = 0; j < R_mortality_cols; j++) {
            R_conflict_mortality[i][j] = REAL(SEXP_conflict_mortality)[i + j * LENGTH(SEXP_conflict_mortality)/R_mortality_cols];
        }
    }

    R_hunting_mortality = malloc(LENGTH(SEXP_hunting_mortality)/R_mortality_cols * sizeof(double*));
    for (int i = 0; i < LENGTH(SEXP_hunting_mortality)/R_mortality_cols; i++) {
        R_hunting_mortality[i] = malloc(R_mortality_cols * sizeof(double));
        for (int j = 0; j < R_mortality_cols; j++) {
            R_hunting_mortality[i][j] = REAL(SEXP_hunting_mortality)[i + j * LENGTH(SEXP_hunting_mortality)/R_mortality_cols];
        }
    }

    R_initial_prides_coalitions = INTEGER(SEXP_initial_prides_coalitions)[0];

    R_K_individuals = INTEGER(SEXP_K_individuals)[0];
    R_K_prides = INTEGER(SEXP_K_prides)[0];
    R_K_coalitions = INTEGER(SEXP_K_coalitions)[0];
    R_K_edged = INTEGER(SEXP_K_edged)[0];

    stats = malloc(sizeof(struct statistics));
    mc_allocate_statistics(stats);
    monte_carlo(stats);

    SEXP R_runs;
    PROTECT(R_runs = allocVector(REALSXP, number_of_months * R_number_mc_runs * NUMBER_OF_STATS)); nprot++;
    for (long i = 0; i < R_number_mc_runs; i++) {
        for (long j = 0; j < number_of_months; j++) {
            for (long k = 0; k < NUMBER_OF_STATS; k++) {
                REAL(R_runs)[ i*number_of_months*NUMBER_OF_STATS + j*NUMBER_OF_STATS + k] = stats->runs[i][j][k];
            }
        }
    }

    SEXP R_individuals;
    t_history *current_hy = stats->history_individuals;
    int nrows = 0;
    while (current_hy != NULL) {
        nrows++;
        current_hy = current_hy->next;
    }

    current_hy = stats->history_individuals;
    PROTECT(R_individuals = allocVector(REALSXP, nrows*number_of_months)); nprot++;
    for (long i = 0; i < nrows; i++) {
        for (long j = 0; j < number_of_months; j++) {
            REAL(R_individuals)[i * number_of_months + j] = current_hy->events_individual[j];
        }
        current_hy = current_hy->next;
    }

    char *names[2] = {"runs", "individuals"};
    SEXP list, list_names;
    PROTECT(list_names = allocVector(STRSXP,2)); nprot++;
    for (int i = 0; i < 2; i++) {
        SET_STRING_ELT(list_names, i , mkChar(names[i]));
    }

    PROTECT(list = allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(list, 0, R_runs);
    SET_VECTOR_ELT(list, 1, R_individuals);
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(nprot);
    
    for (int i = 0; i < 2; i++) {
        free(R_survival_av[i]);
    }
    free(R_survival_av);

    free(R_litter_distribution);
    
    for (int i = 0; i < LENGTH(SEXP_conflict_mortality)/R_mortality_cols; i++) {
        free(R_conflict_mortality[i]);
    }
    free(R_conflict_mortality);

    for (int i = 0; i < LENGTH(SEXP_hunting_mortality)/R_mortality_cols; i++) {
        free(R_hunting_mortality[i]);
    }
    free(R_hunting_mortality);
    
    mc_free_statistics(stats);

    return(list);

}
