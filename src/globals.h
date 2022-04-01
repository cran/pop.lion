/* globals.h
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
 * 'pop.wolf' is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with 'pop.wolf'. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GLOBALS_H
#define GLOBALS_H

#include "tools.h"

extern long R_number_of_years;
extern long R_number_mc_runs;

extern long number_of_months;

extern struct statistics *stats;

extern double **R_survival_av;
extern double *R_litter_distribution;

extern int R_conflict_age_male;
extern int R_conflict_age_female;
extern int R_hunting_age_male;
extern int R_hunting_age_female;
extern int R_mortality_cols;

extern double **R_conflict_mortality;
extern double **R_hunting_mortality;

extern int R_hunter_error;

extern int R_initial_prides_coalitions;
extern int R_K_individuals;
extern int R_K_prides;
extern int R_K_coalitions;
extern int R_K_edged;

#endif
