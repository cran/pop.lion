/* pop.h
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

#ifndef POP_H
#define POP_H

#include "mc.h"
#include "tools.h"

#define MAX_INDIV 10000

#define NUMBER_OF_STAGES 5
#define NUMBER_OF_SEXES 2
#define NUMBER_OF_EVENTS 5
#define NUMBER_OF_STATS 35

#define FEMALE 0
#define MALE 1

#define AGE_OLDEST 15*12
#define IDX_NINDIV 0

#define IDX_NINDIV_CORE 28
#define IDX_NINDIV_EDGE 29

#define IDX_NFEMALES 11
#define IDX_NMALES 12

#define IDX_NFEMALES_CORE 23
#define IDX_NFEMALES_EDGE 25

#define IDX_NMALES_CORE 34
#define IDX_NMALES_EDGE 33

#define IDX_NCUBS_CORE 24
#define IDX_NCUBS_EDGE 26

#define IDX_NFEMALES_ADULT 17
#define IDX_NMALES_ADULT 16

#define IDX_NCOALIS 2
#define IDX_NPRIDES 1

#define IDX_NPRIDES_PER_COALI 31

#define IDX_NCOALIS_RESIDENT 3
#define IDX_NCOALIS_VAGRANT 4

#define IDX_NPRIDES_RESIDENT 5
#define IDX_NPRIDES_VAGRANT 6

#define IDX_COALISIZE_RESIDENT 7
#define IDX_COALISIZE_VAGRANT 8

#define IDX_PRIDESIZE_RESIDENT 9
#define IDX_PRIDESIZE_VAGRANT 10

#define IDX_NPRIDES_CORE 20
#define IDX_NPRIDES_EDGED 19

#define IDX_PRIDESIZE_CORE 21
#define IDX_PRIDESIZE_EDGED 22

#define IDX_LITTERS 14
#define IDX_AGE 15
#define IDX_IBI 32

#define IDX_NFEMALES_DISPERSED 18

#define IDX_TAKEOVERS 13
#define IDX_COALITION_TENURE 27

#define IDX_FAILED_HUNTS 30

typedef struct t_individual t_individual;
typedef struct t_pride t_pride;
typedef struct t_coalition t_coalition;
typedef struct t_population t_population;

struct t_individual {
    int unique;
    int alive;
    int sex;
    int age;
    int stage;
    int dispersed;
    int mated;
    int age_mated;
    int is_mother;
    int age_mother;
    int estimated_age;
    int *events;

    t_individual *mother;
    GPtrArray *litter;

    t_pride *my_pride;
    t_coalition *my_coalition;

    t_individual *previous;
    t_individual *next;
};

struct t_pride {
    int stage;
    int age_resident;
    int is_edged;
    int age_vagrant;

    GPtrArray *all_members;
    t_coalition *the_coalition;

    t_pride *previous;
    t_pride *next;
};

struct t_coalition {
    int stage;
    int age_resident;
    int age_vagrant;
    int possible_fights;

    GPtrArray *all_members;
    GPtrArray *the_prides;

    t_coalition *previous;
    t_coalition *next;
};

struct t_population {

    int number_indiv;
    int number_indiv_history;

    double *live_stats;

    int number_prides;
    int number_prides_settled;
    int number_prides_edged;

    int number_coalitions;
    int number_coalitions_settled;

    t_individual *all_indiv;
    t_pride *all_prides;
    t_coalition *all_coalitions;

    double survival[NUMBER_OF_SEXES][AGE_OLDEST];

    int initial_prides_coalitions;
    int K_individuals;
    int K_prides;
    int K_coalitions;
    int K_edged;

    int total_failed_hunts;
    int month;

};

void create_initial_population(t_population *pop);

void set_population_parameters(t_population *pop);
void set_deterministic_parameters(t_population *pop);
void set_stochastic_parameters(t_population *pop);

void cycle_year(t_population *pop, long the_seed, long the_year, struct statistics *stats);
void do_statistics(t_population *pop, long the_run, long the_year, struct statistics *stats);
void collect_events(t_population *pop, struct statistics *stats);
void free_population(t_population *pop);

#endif
