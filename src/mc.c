/* mc.c
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
#include <stdlib.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "mc.h"
#include "globals.h"

/////////////////////////////////////////////////////////////////////////////////
// Allocate statistics
/////////////////////////////////////////////////////////////////////////////////

void mc_allocate_statistics(struct statistics *stats) {
    
    stats->runs = malloc(R_number_mc_runs * sizeof(double **));
    
    for (long i = 0; i < R_number_mc_runs; i++) {
        stats->runs[i] = malloc(number_of_months * sizeof(double *));
        for (long j = 0; j < number_of_months; j++) {
            stats->runs[i][j] = malloc(NUMBER_OF_STATS * sizeof(double));
            for (long k = 0; k < NUMBER_OF_STATS; k++) {
                stats->runs[i][j][k] = 0.0;
            }
        }
    }
    
    stats->history_individuals = malloc(sizeof(t_history));
    stats->history_individuals = NULL;
    
}

/////////////////////////////////////////////////////////////////////////////////
// Free results
/////////////////////////////////////////////////////////////////////////////////

void mc_free_results(struct statistics *stats) {
    
    for (long i = 0; i < R_number_mc_runs; i++) {
        for (long j = 0; j < number_of_months; j++) {
            free(stats->runs[i][j]);
        }
        free(stats->runs[i]);
    }
    free(stats->runs);
    
    t_history *next_hy;
    while (stats->history_individuals != NULL) {
        next_hy = stats->history_individuals->next;
        free(stats->history_individuals->events_individual);
        free(stats->history_individuals);
        stats->history_individuals = next_hy;
    }
    free(stats->history_individuals);
    
    free(stats);
    
}

/////////////////////////////////////////////////////////////////////////////////
// MONTE CARLO
/////////////////////////////////////////////////////////////////////////////////

void monte_carlo(struct statistics *stats) {
    
    GetRNGstate();
    
    long steps = R_number_mc_runs/50;
    
    Rprintf("\n|");
    
    for (long i = 0; i < R_number_mc_runs; i++) {
        
        t_population *pop = malloc(sizeof(t_population));
        
        set_population_parameters(pop);
        
        set_deterministic_parameters(pop);
        
        create_initial_population(pop);
        
        do_statistics(pop, i, 0, stats);
        
        for (long j = 1; j <= R_number_of_years; j++) {
            
            // set_stochastic_parameters(pop);
            
            cycle_year(pop, i, j, stats);
            
            if (pop->number_indiv == 0) {
                break;
            }
            
        }
        
        collect_events(pop, stats);
        
        if (steps > 0) {
            if (i % steps == 0) {
                Rprintf("*");
            }
        }
        
        free_population(pop);
        free(pop);
        
    }
    
    if (steps > 0) {
        Rprintf("|");
    }
    
    PutRNGstate();
    
}
