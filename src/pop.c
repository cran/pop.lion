/* pop.c
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

#include "pop.h"
#include "globals.h"

/////////////////////////////////////////////////////////////////////////////////
#pragma mark MACROS
/////////////////////////////////////////////////////////////////////////////////

#define CUB 0
#define YEARLING 1
#define SUBADULT_23 2
#define SUBADULT_34 3
#define ADULT 4

#define VAGRANT 0
#define RESIDENT 1

#define NOT_EDGED 0
#define EDGED 1

#define NATAL 0
#define DISPERSE 1

#define NOT_MATED 0
#define MATED 1

#define NOT_MOTHER 0
#define IS_MOTHER 1

#define MAX_GROUP_SIZE 50

#define SEX_RATIO 0.5

#define PRIDE_MEAN_SIZE_ADULT_F 4
#define PRIDE_MEAN_AGE_ADULT_F 72
#define CUB_PER_ADULT_F 0.34
#define YEARLING_PER_ADULT_F 0.74
#define SUBADULT23_PER_ADULT_F 0.74
#define SUBADULT34_PER_ADULT_F 0.74

#define COALI_MEAN_SIZE_ADULT_M 2.6
#define COALI_MEAN_AGE_ADULT_M 72
#define COALI_MAX_SIZE 4

#define AGE_1ST_REPRODUCTION 48
#define AGE_CUB_DIE_WITHOUT_MOTHER 12
#define AGE_CUB_DIE_WITHOUT_SURROGATE 24
#define AGE_SURROGATE 48

#define AGE_SURVIVE_INFANTICIDE 32

#define AGE_DISPERSAL_F_BEGIN 28 // Females leave natal pride if pride is larger than MAX_PRIDE_SIZE
#define AGE_DISPERSAL_F_BEGIN_TRIGGER 20 // < females less than 20 months avoid new males but stay with the pride (Hanby and Bygott 1987)
#define AGE_DISPERSAL_F_END 45
#define AGE_DISPERSAL_F_END_TRIGGER 30 // < We added this so that fights will only trigger dispersal of females younger than this age (Hanby and Bygott 1987) < 20 = stay; > 20 = disperse
#define AGE_DISPERSAL_M_BEGIN 28
#define AGE_DISPERSAL_M_END 45

#define AGE_MAX_PRIDE_VAGRANT 6
#define AGE_MAX_COALI_VAGRANT 6

#define MORTALITY_NOTSETTLE 0.0 // < Probability of members in a vagrant pride dying every month if they don't settle (alternative to a fixed time as in AGE_MAX_PRIDE_VAGRANT)
#define MORTALITY_FIGHT_M 0.15 // < Probability of members of a losing coalition dying after a fight
#define MORTALITY_FIGHT_F 0.15 // < Probablity of mothers being killed during infanticide
#define NUMBER_SPACES_SEARCHED 1 // < How many vacant spaces do prides search for each month when trying to settle
#define NUMBER_PRIDES_SEARCHED 1 // < How many resident prides with no coalitions do vagrant coalitions search for when trying to settle
#define INFANTICIDE 1 // < Probability of a cub dying from infanticide - in case not all cubs are killed
#define RESIDENT_ADVANTAGE 1.5 // < The advantage factor for resident coalitions in fights - resident coalitions have an advantage because of the pride
#define MAX_PRIDE_SIZE 10 // < Maximum pride size which triggers dispersal
#define FIGHT_CHANCE 1 // < Probability of a vagrant coalition engaging in a fight each month - in case we want to make fights less frequent
#define MATING_SUCCESS 0.25 // < Mating might not always lead to conception - this allows for a lag when mating is unsuccessful
#define STARVATION_1YR 0.03 // < Probability of cubs (< 1yr) dying from starvation when the population is at carrying capacity
#define STARVATION_2YR 0.01 // < Probability of cubs (> 1yr) dying from starvation when the population is at carrying capacity

#define PARENT 8
#define PREGNANT 7
#define MATE 6
#define UNSETTLED 5
#define SETTLED 4
#define DISPERSED 3
#define BORN 2
#define ALIVE 1
#define DIED -1
#define DIED_FOUGHT -2 // Never used
#define DIED_CONFLICT -3
#define DIED_HUNTING -4
#define DIED_INFANTICIDE_FIGHT -5
#define DIED_INFANTICIDE_FREE -6
#define DIED_BACKGROUND -7
#define KILLED_FIGHT -8
#define KILLED_FIGHT_F -9
#define DIED_STARVATION -10
#define DIED_INFANTICIDE_LATE -11

void probe(t_population *pop) {

    t_individual *current_idv = pop->all_indiv;

    int x = 0;

    while (current_idv != NULL) {

        x = x + current_idv->age;

        current_idv = current_idv->next;

    }
    Rprintf(" %d, %d, %3.2f | ", x, pop->number_indiv, (double)x/pop->number_indiv);
}

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION PROTOTYPES
/////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 Comparing functions
 *******************************************************************************/

gint compare_pride_size (void *p1, void *p2);
gint compare_coalition_size (void *c1, void *c2);
gint compare_nprides_btw_coalitions (void *f1, void *f2);
gint compare_age (void *i1, void *i2);

/*******************************************************************************
 Parameters
 *******************************************************************************/

void set_population_parameters(t_population *pop);
void set_deterministic_parameters(t_population *pop);
void set_stochastic_parameters(t_population *pop);

/*******************************************************************************
 Create structures
 *******************************************************************************/

t_individual* create_individual(t_population *pop, int the_sex, int the_age, int the_stage);
t_pride* create_pride(t_population *pop, int the_status);
t_coalition* create_coalition(t_population *pop, int the_status);

/*******************************************************************************
 Create initial population
 *******************************************************************************/

t_pride* create_initial_pride(t_population *pop);
t_coalition* create_initial_coalition(t_population *pop);
void create_initial_population(t_population *pop);

/*******************************************************************************
 Individual connections
 *******************************************************************************/

void individual_update_events(t_individual *the_idv, long the_month, int event);

void cub_gets_mother(t_individual *the_cub, t_individual *the_mother);
void cubs_lose_mother(t_individual *the_mother);
void mother_loses_cub(t_individual *the_cub);

void individual_joins_pride(t_individual *the_indiv, t_pride *the_pride);
void individual_joins_coalition(t_individual *the_indiv, t_coalition *the_coalition);

void individual_leaves_pride(t_individual *the_indiv, t_pride *the_pride);
void individual_leaves_coalition(t_individual *the_indiv, t_coalition *the_coalition);

/*******************************************************************************
 Group connections
 *******************************************************************************/

t_pride* pride_leaves_pop(t_population *pop, t_pride *current_pride);
t_coalition* coalition_leaves_pop(t_population *pop, t_coalition *current_coali);

/*******************************************************************************
 Individuals
 *******************************************************************************/

void individuals_die(t_population *pop, long the_month);
void individuals_die_inoldprides(t_population *pop, long the_month);
void individuals_die_inoldcoalitions(t_population *pop, long the_month);
void individuals_remove(t_population *pop);
int get_individual_edgedrisk(t_population *pop, t_individual *the_idv);
void individuals_hunting(t_population *pop, long the_month);
void individuals_age(t_population *pop);
void individuals_disperse(t_population *pop, long the_month);

/*******************************************************************************
 Prides and coalitions
 *******************************************************************************/

void prides_remove(t_population *pop);
void coalitions_remove(t_population *pop);
void coalitions_age(t_population *pop);
void prides_age(t_population *pop);

/*******************************************************************************
 Prides
 *******************************************************************************/

void prides_settle(t_population *pop, long the_month);
void prides_reproduce(t_population *pop, long the_month);

/*******************************************************************************
 Coalitions
 *******************************************************************************/

void coalitions_split(t_population *pop);
void coalitions_meet_prides(t_population *pop, long the_month);
void coalitions_merge(t_population *pop);
void coalitions_fight(t_population *pop, long the_month);

/*******************************************************************************
 Model
 *******************************************************************************/

void cycle_year(t_population *pop, long the_seed, long the_year, struct statistics *stats);
void do_statistics(t_population *pop, long the_run, long the_month, struct statistics *stats);

/*******************************************************************************
 Free memory
 *******************************************************************************/

void individual_free(t_individual *idv);
void free_population(t_population *pop);

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION DEFINITIONS
/////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 Comparing functions
 *******************************************************************************/

gint compare_pride_size (void *p1, void *p2) {

    t_pride *pride1 = *(t_pride**) p1;

    t_pride *pride2 = *(t_pride**) p2;

    return (gint) (pride1->all_members->len - pride2->all_members->len);

}

gint compare_coalition_size (void *c1, void *c2) {

    t_coalition *coali1 = *(t_coalition**) c1;

    t_coalition *coali2 = *(t_coalition**) c2;

    return (gint) (coali1->all_members->len - coali2->all_members->len);

}

gint compare_nprides_btw_coalitions (void *f1, void *f2) {

    t_coalition *coali1 = *(t_coalition**) f1;

    t_coalition *coali2 = *(t_coalition**) f2;

    return (gint) (coali1->the_prides->len - coali2->the_prides->len);

}

gint compare_age (void *i1, void *i2) {

    t_individual *indiv1 = *(t_individual**) i1;

    t_individual *indiv2 = *(t_individual**) i2;

    return (gint) (indiv1->age - indiv2->age);

}

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : PARAMETERS
/////////////////////////////////////////////////////////////////////////////////

void set_population_parameters(t_population *pop) {

    pop->initial_prides_coalitions = R_initial_prides_coalitions;
    pop->K_individuals = R_K_individuals;
    pop->K_prides = R_K_prides;
    pop->K_coalitions = R_K_coalitions;
    pop->K_edged = R_K_edged;

}

void set_deterministic_parameters(t_population *pop) {

    for (int i = 0; i < AGE_OLDEST; i++) {

        pop->survival[FEMALE][i] = R_survival_av[FEMALE][i];
        pop->survival[MALE][i] = R_survival_av[MALE][i];

    }

    //pop->litter_size = R_litter_size_av;

}

/* missing */ void set_stochastic_parameters(t_population *pop) {
    /*
     double mu = 0;
     double sigma = 0;
     double randomv = 0;

     // Pup survival (0-6 months)
     mu = R_survival_av_PUP;
     sigma = R_survival_sd_PUP;
     if (sigma == 0) {
     randomv = mu;
     } else {
     randomv = rbeta(beta_shape(mu, sigma), beta_rate(mu, sigma));
     }
     pop->survival[PUP] = pow(randomv, (double)1/6);

     // Subadult survival
     mu = R_survival_av_SUBADULT;
     sigma = R_survival_sd_SUBADULT;
     if (sigma == 0) {
     randomv = mu;
     } else {
     randomv = rbeta(beta_shape(mu, sigma), beta_rate(mu, sigma));
     }
     pop->survival[SUBADULT] = pow(randomv, (double)1/12);

     // Vagrant survival
     mu = R_survival_av_VAGRANT;
     sigma = R_survival_sd_VAGRANT;
     if (sigma == 0) {
     randomv = mu;
     } else {
     randomv = rbeta(beta_shape(mu, sigma), beta_rate(mu, sigma));
     }
     pop->survival[VAGRANT] = pow(randomv, (double)1/12);

     // ADULT survival
     mu = R_survival_av_ALPHA;
     sigma = R_survival_sd_ALPHA;
     if (sigma == 0) {
     randomv = mu;
     } else {
     randomv = rbeta(beta_shape(mu, sigma), beta_rate(mu, sigma));
     }
     pop->survival[ALPHA] = pow(randomv, (double)1/12);

     // Litter size
     mu = R_litter_size_av;
     sigma = R_litter_size_sd;
     pop->litter_size = rgamma(gamma_shape(mu, sigma), 1/gamma_rate(mu, sigma));
     */
}

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : CREATE ENTITIES
/////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 Create individual
 *******************************************************************************/

t_individual* create_individual(t_population *pop, int the_sex, int the_age, int the_stage) {

    t_individual *new_idv = malloc(sizeof(t_individual));

    pop->number_indiv++;
    pop->number_indiv_history++;

    new_idv->unique = pop->number_indiv_history;
    new_idv->alive = 1;
    new_idv->sex = the_sex;
    new_idv->stage = the_stage;
    new_idv->age = imin2(the_age, AGE_OLDEST);
    new_idv->dispersed = NATAL;
    new_idv->mated = NOT_MATED;
    new_idv->age_mated = 0;
    new_idv->is_mother = NOT_MOTHER;
    new_idv->age_mother = 0;
    new_idv->estimated_age = imin2(the_age, AGE_OLDEST);

    new_idv->events = malloc((12*R_number_of_years+1)*sizeof(int));
    for (int i = 0; i < (12*R_number_of_years+1); i++) {
        new_idv->events[i] = 0;
    }

    new_idv->mother = NULL;
    new_idv->my_pride = NULL;
    new_idv->my_coalition = NULL;
    new_idv->litter = g_ptr_array_sized_new(MAX_GROUP_SIZE);

    if (pop->number_indiv == 1) {

        new_idv->previous = NULL;
        new_idv->next = NULL;
        pop->all_indiv = new_idv;

    } else {

        new_idv->previous = NULL;
        new_idv->next = pop->all_indiv;
        new_idv->next->previous = new_idv;
        pop->all_indiv = new_idv;

    }

    return(new_idv);

}

/*******************************************************************************
 Create pride
 *******************************************************************************/

t_pride* create_pride(t_population *pop, int the_status) {

    t_pride *a_pride = malloc(sizeof(t_pride));

    pop->number_prides++;
    if (the_status == RESIDENT) {
        pop->number_prides_settled++;

        // If the number of prides settled and not at the edge but at the core (so the total number of prides settled minus the number of prides settled at the edge) is less than the space available at the core (so total K minus K at the edge), then prides go to the core (and therefore not the edge).
        if (pop->number_prides_settled - pop->number_prides_edged <= pop->K_prides - pop->K_edged ) {
            a_pride->is_edged = NOT_EDGED;
        } else {
            a_pride->is_edged = EDGED;
            pop->number_prides_edged++;
        }
    } else {
        a_pride->is_edged = NOT_EDGED;
    }

    a_pride->stage = the_status;
    a_pride->age_resident = 0;
    a_pride->age_vagrant = 0;

    a_pride->all_members = g_ptr_array_sized_new(MAX_GROUP_SIZE);
    a_pride->the_coalition = NULL;

    if (pop->number_prides == 1) {

        a_pride->previous = NULL;
        a_pride->next = NULL;
        pop->all_prides = a_pride;

    } else {

        a_pride->previous = NULL;
        a_pride->next = pop->all_prides;
        a_pride->next->previous = a_pride;
        pop->all_prides = a_pride;

    }

    return(a_pride);

}

/*******************************************************************************
 Create coalition
 *******************************************************************************/

t_coalition* create_coalition(t_population *pop, int the_status) {

    t_coalition *a_coali = malloc(sizeof(t_coalition));

    pop->number_coalitions++;

    if (the_status == RESIDENT) {
        pop->number_coalitions_settled++;
    }

    a_coali->stage = the_status;
    a_coali->age_resident = 0;
    a_coali->age_vagrant = 0;

    a_coali->all_members = g_ptr_array_sized_new(MAX_GROUP_SIZE);
    a_coali->the_prides = g_ptr_array_sized_new(MAX_GROUP_SIZE);

    if (pop->number_coalitions == 1) {

        a_coali->previous = NULL;
        a_coali->next = NULL;
        pop->all_coalitions = a_coali;

    } else {

        a_coali->previous = NULL;
        a_coali->next = pop->all_coalitions;
        a_coali->next->previous = a_coali;
        pop->all_coalitions = a_coali;

    }

    return(a_coali);

}

/*******************************************************************************
 Create initial population
 *******************************************************************************/

t_pride* create_initial_pride(t_population *pop) {

    t_pride *a_pride = malloc(sizeof(t_pride));

    pop->number_prides++;

    // distribute initial prides between edge and core with equal probablity

    if (rbinom(1,0.5) == 1){

        if (pop->number_prides_settled - pop->number_prides_edged < pop->K_prides - pop->K_edged ) {
            a_pride->is_edged = NOT_EDGED;
        } else {
            a_pride->is_edged = EDGED;
            pop->number_prides_edged++;
        }

    } else {

        if (pop->number_prides_edged < pop->K_edged ) {
            a_pride->is_edged = EDGED;
            pop->number_prides_edged++;
        } else {
            a_pride->is_edged = NOT_EDGED;
        }

    }

    a_pride->stage = RESIDENT;
    pop->number_prides_settled++;

    a_pride->age_resident = 0;
    a_pride->age_vagrant = 0;

    a_pride->all_members = g_ptr_array_sized_new(MAX_GROUP_SIZE);
    a_pride->the_coalition = NULL;

    int size1 = rpois(PRIDE_MEAN_SIZE_ADULT_F);
    int size2;

    for (int i = 0; i < size1; i++) {

        t_individual *idv_adult = create_individual(pop, FEMALE, rpois(PRIDE_MEAN_AGE_ADULT_F), ADULT);
        individual_joins_pride(idv_adult, a_pride);

        size2 = rpois(CUB_PER_ADULT_F);
        for (int j = 0; j < size2; j++) {

            t_individual *idv_cub = create_individual(pop, (rbinom(1, SEX_RATIO) == 1)?MALE:FEMALE, 0, CUB);
            individual_joins_pride(idv_cub, a_pride);
            cub_gets_mother(idv_cub, idv_adult);

        }

        size2 = rpois(YEARLING_PER_ADULT_F);
        for (int j = 0; j < size2; j++) {

            t_individual *idv_yearling = create_individual(pop, (rbinom(1, SEX_RATIO) == 1)?MALE:FEMALE, 1*12, YEARLING);
            individual_joins_pride(idv_yearling, a_pride);
            cub_gets_mother(idv_yearling, idv_adult);

        }

        size2 = rpois(SUBADULT23_PER_ADULT_F);
        for (int j = 0; j < size2; j++) {

            t_individual *idv_sub23 = create_individual(pop, (rbinom(1, SEX_RATIO) == 1)?MALE:FEMALE, 2*12, SUBADULT_23);
            individual_joins_pride(idv_sub23, a_pride);
            cub_gets_mother(idv_sub23, idv_adult);

        }

        size2 = rpois(SUBADULT34_PER_ADULT_F);
        for (int j = 0; j < size2; j++) {

            t_individual *idv_sub34 = create_individual(pop, (rbinom(1, SEX_RATIO) == 1)?MALE:FEMALE, 3*12, SUBADULT_34);
            individual_joins_pride(idv_sub34, a_pride);
            cub_gets_mother(idv_sub34, idv_adult);

        }

    }

    if (pop->number_prides == 1) {

        a_pride->previous = NULL;
        a_pride->next = NULL;
        pop->all_prides = a_pride;

    } else {

        a_pride->previous = NULL;
        a_pride->next = pop->all_prides;
        a_pride->next->previous = a_pride;
        pop->all_prides = a_pride;

    }

    return(a_pride);

}

t_coalition* create_initial_coalition(t_population *pop) {

    t_coalition *a_coali = malloc(sizeof(t_coalition));

    pop->number_coalitions++;

    a_coali->stage = RESIDENT;
    pop->number_coalitions_settled++;

    a_coali->age_resident = 5;
    a_coali->age_vagrant = 0;

    a_coali->all_members = g_ptr_array_sized_new(MAX_GROUP_SIZE);
    a_coali->the_prides = g_ptr_array_sized_new(MAX_GROUP_SIZE);

    int size = 2;

    for (int i = 0; i < size; i++) {

        t_individual *a_indiv = create_individual(pop, MALE, rpois(COALI_MEAN_AGE_ADULT_M), ADULT);
        individual_joins_coalition(a_indiv, a_coali);

    }

    if (pop->number_coalitions == 1) {

        a_coali->previous = NULL;
        a_coali->next = NULL;
        pop->all_coalitions = a_coali;

    } else {

        a_coali->previous = NULL;
        a_coali->next = pop->all_coalitions;
        a_coali->next->previous = a_coali;
        pop->all_coalitions = a_coali;

    }

    return(a_coali);

}

void create_initial_population(t_population *pop) {

    pop->number_indiv = 0;
    pop->number_indiv_history = 0;
    pop->all_indiv = NULL;

    pop->number_prides = 0;
    pop->number_prides_settled = 0;
    pop->number_prides_edged = 0;
    pop->all_prides = NULL;

    pop->number_coalitions = 0;
    pop->number_coalitions_settled = 0;
    pop->all_coalitions = NULL;

    pop->total_failed_hunts = 0;
    pop->month = 0;

    // Add with 1 to 1 pride-coalition

    for (int i = 0; i < pop->initial_prides_coalitions; i++) {

        t_pride *a_pride = create_initial_pride(pop);
        t_coalition *a_coalition = create_initial_coalition(pop);

        a_pride->the_coalition = a_coalition;
        g_ptr_array_add(a_coalition->the_prides, a_pride);

    }

    pop->live_stats = malloc(NUMBER_OF_STATS*sizeof(double));
    for (int i = 0; i < NUMBER_OF_STATS; i++) {
        pop->live_stats[i] = 0;
    }

}

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : CONNECT INDIVIDUALS
/////////////////////////////////////////////////////////////////////////////////

void cub_gets_mother(t_individual *the_cub, t_individual *the_mother) {

    the_cub->mother = the_mother;
    g_ptr_array_add(the_mother->litter, the_cub);

}

void cubs_lose_mother(t_individual *the_mother) {

    if (the_mother->litter->len != 0) {
        for (int i = 0; i < the_mother->litter->len; i++) {
            ((t_individual*)g_ptr_array_index(the_mother->litter, i))->mother = NULL;
        }
        g_ptr_array_empty(the_mother->litter);
    }

}

void mother_loses_cub(t_individual *the_cub) {

    if (the_cub->mother != NULL) {
        g_ptr_array_remove_fast(the_cub->mother->litter, the_cub);
    }
    the_cub->mother = NULL;

}

void individual_joins_pride(t_individual *the_indiv, t_pride *the_pride) {

    the_indiv->my_pride = the_pride;
    the_indiv->my_coalition = NULL;

    g_ptr_array_add(the_pride->all_members, the_indiv);

}

void individual_joins_coalition(t_individual *the_indiv, t_coalition *the_coalition) {

    the_indiv->my_pride = NULL;
    the_indiv->my_coalition = the_coalition;

    g_ptr_array_add(the_coalition->all_members, the_indiv);

}

void individual_leaves_pride(t_individual *the_indiv, t_pride *the_pride) {

    if (the_indiv->my_pride != NULL) {

        if (the_indiv->mother != NULL) {
            g_ptr_array_remove_fast(the_indiv->mother->litter, the_indiv);
        }
        the_indiv->mother = NULL;

        g_ptr_array_remove_fast(the_pride->all_members, the_indiv);
        the_indiv->my_pride = NULL;

    }

}

void individual_leaves_coalition(t_individual *the_indiv, t_coalition *the_coalition) {

    if (the_indiv->my_coalition != NULL) {

        g_ptr_array_remove_fast(the_coalition->all_members, the_indiv);
        the_indiv->my_coalition = NULL;

    }

}

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : CONNECT GROUPS
/////////////////////////////////////////////////////////////////////////////////

void pride_changes_coalition(t_pride *current_pride, t_coalition *new_coali) {

    g_ptr_array_remove_fast(current_pride->the_coalition->the_prides, current_pride);
    current_pride->the_coalition = new_coali;
    g_ptr_array_add(new_coali->the_prides, current_pride);

}

t_pride* pride_leaves_pop(t_population *pop, t_pride *current_pride) {

    t_pride *next_pride = current_pride->next;

    if (current_pride->stage == RESIDENT) {
        pop->number_prides_settled--;
    }

    if (current_pride->is_edged == EDGED
        && current_pride->stage == RESIDENT) { // have to include "resident" otherwise could could create a space when a vagrant pride leaves
        pop->number_prides_edged--;
    }

    // Tell members their pride is gone
    for (int i = 0; i < current_pride->all_members->len; i++) {
        ((t_individual*)g_ptr_array_index(current_pride->all_members, i))->my_pride = NULL;
    }
    g_ptr_array_empty(current_pride->all_members);

    // Tell coalition this pride is gone
    if (current_pride->the_coalition != NULL) {

        if (current_pride->the_coalition->the_prides->len == 1){
            current_pride->the_coalition->stage = VAGRANT;
            pop->number_coalitions_settled--;
        }
        g_ptr_array_remove_fast(current_pride->the_coalition->the_prides, current_pride);

    }

    // Tell pride it no longer has a coalition
    current_pride->the_coalition = NULL;

    // If the coalition is at both the beginning and the end of the list
    if ( (current_pride->previous == NULL) & (current_pride->next == NULL) ) {
        pop->all_prides = NULL;
    }

    // If the coalition is at the beginning of the list
    else if ( (current_pride->previous == NULL) & (current_pride->next != NULL) ) {
        current_pride->next->previous = NULL; // The next coalition it has no previous (it is the beginning))
        pop->all_prides = current_pride->next; // Update the start of the list
    }

    // If the coalition is at the end of the list
    else if ( (current_pride->previous != NULL) & (current_pride->next == NULL) ) {
        current_pride->previous->next = NULL; // Connect from previous to next
    }

    // Otherwise
    else {
        current_pride->next->previous = current_pride->previous; // Connect from next to previous
        current_pride->previous->next = current_pride->next; // Connect from previous to next
    }

    // Free memory
    free(current_pride);

    pop->number_prides--;

    return(next_pride);

}

t_coalition* coalition_leaves_pop(t_population *pop, t_coalition *current_coali) {

    t_coalition *next_coali = current_coali->next;

    if (current_coali->stage == RESIDENT) {
        pop->number_coalitions_settled--;
    }

    // Tell members their coalition is gone
    for (int i = 0; i < current_coali->all_members->len; i++) {
        ((t_individual*)g_ptr_array_index(current_coali->all_members, i))->my_coalition = NULL;
    }
    g_ptr_array_empty(current_coali->all_members);

    // Tell prides their coalition is gone
    for (int i = 0; i < current_coali->the_prides->len; i++) {
        ((t_pride*)g_ptr_array_index(current_coali->the_prides, i))->the_coalition = NULL;
    }

    // Tell coalition it no longer has prides
    g_ptr_array_free(current_coali->the_prides);

    // If the coalition is at both the beginning and the end of the list
    if ( (current_coali->previous == NULL) & (current_coali->next == NULL) ) {
        pop->all_coalitions = NULL;
    }

    // If the coalition is at the beginning of the list
    else if ( (current_coali->previous == NULL) & (current_coali->next != NULL) ) {
        current_coali->next->previous = NULL; // The next coalition it has no previous (it is the beginning))
        pop->all_coalitions = current_coali->next; // Update the start of the list
    }

    // If the coalition is at the end of the list
    else if ( (current_coali->previous != NULL) & (current_coali->next == NULL) ) {
        current_coali->previous->next = NULL; // Connect from previous to next
    }

    // Otherwise
    else {
        current_coali->next->previous = current_coali->previous; // Connect from next to previous
        current_coali->previous->next = current_coali->next; // Connect from previous to next
    }

    // Free memory
    free(current_coali);

    pop->number_coalitions--;

    return(next_coali);

}

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : INDIVIDUAL EVENTS
/////////////////////////////////////////////////////////////////////////////////

void individuals_die(t_population *pop, long the_month) {

    // This only updates 'alive'

    t_individual *current_idv = pop->all_indiv;
    t_individual *surrogate_idv = NULL;
    t_individual *cub_idv = NULL;
    t_pride *current_pride = NULL;

    double the_surv;
    int pride_has_surrogate_mum = 0;
    int rdie = -1;

    int individuals_core = 0;
    int individuals_edge = 0;

    double K_indiv_core;
    double K_indiv_edge;

    K_indiv_core = (((double)pop->K_prides - pop->K_edged) / (double)pop->K_prides) * (double)pop->K_individuals;
    K_indiv_edge = ((double)pop->K_edged / pop->K_prides) * (double)pop->K_individuals;


    while (current_idv != NULL) {
        if (get_individual_edgedrisk(pop, current_idv) == 0) {
            individuals_core++;
        } else if (get_individual_edgedrisk(pop, current_idv) == 1) {
            individuals_edge++;
        }
        current_idv = current_idv->next;
    }

    current_idv = pop->all_indiv;

    while (current_idv != NULL) {

        int cub_die = 0;
        int cub_killed = 0;

        if (current_idv->alive == 1) {

            the_surv = pop->survival[current_idv->sex][current_idv->age];

            rdie = rbinom(1, the_surv);

            current_idv->alive = rdie;

            // Senescence

            if (current_idv->age >= AGE_OLDEST) {

                current_idv->alive = 0;

            }

            // What happens to cubs

            if (current_idv->alive == 0 && current_idv->litter->len > 0) {

                pride_has_surrogate_mum = 0;

                // Cub will die if:
                // 1) younger than 12 months,
                // 2) younger than 24 months and no mother older than 4 year in pride

                current_pride = current_idv->my_pride;

                for (int i = 0; i < current_pride->all_members->len; i++) {

                    surrogate_idv = g_ptr_array_index(current_pride->all_members, i);

                    if ( surrogate_idv->sex == FEMALE && surrogate_idv->age > AGE_SURROGATE ) {

                        pride_has_surrogate_mum = 1;

                        break;

                    }

                }

                // Kill cubs if their mother died

                for (int i = 0; i < current_idv->litter->len; i++) {

                    cub_idv = g_ptr_array_index(current_idv->litter, i);

                    if (cub_idv->age <= AGE_CUB_DIE_WITHOUT_MOTHER) {

                        cub_idv->alive = 0;

                    } else if (cub_idv->age > AGE_CUB_DIE_WITHOUT_MOTHER
                               && cub_idv->age <= AGE_CUB_DIE_WITHOUT_SURROGATE
                               && pride_has_surrogate_mum == 0) {

                        cub_idv->alive = 0;

                    }

                }

            }

            // Cubs die from starvation when population is near carrying capacity - this is a common cause of mortality in the Serengeti which is near capacity.

            // Sarvation in the CORE

            double pop_density;

            pop_density = (double)individuals_core/K_indiv_core;

            if (current_idv->age <= 12
                && get_individual_edgedrisk(pop, current_idv) == 0
                && current_idv->alive == 1) {

                cub_die = rbinom(1,(pop_density*STARVATION_1YR));

                if (cub_die == 1) {
                    current_idv->alive = 0;
                }
            }

            if (current_idv->age > 12
                && current_idv->age <= 24
                && get_individual_edgedrisk(pop, current_idv) == 0
                && current_idv->alive == 1) {

                cub_die = rbinom(1,(pop_density*STARVATION_2YR));

                if (cub_die == 1) {
                    current_idv->alive = 0;
                }
            }

            // Starvation in the EDGE

            pop_density = (double)individuals_edge/K_indiv_edge;

            if (current_idv->age <= 12
                && get_individual_edgedrisk(pop, current_idv) == 1
                && current_idv->alive == 1) {

                cub_die = rbinom(1,(pop_density*STARVATION_1YR));

                if (cub_die == 1) {
                    current_idv->alive = 0;
                }
            }

            if (current_idv->age > 12
                && current_idv->age <= 24
                && get_individual_edgedrisk(pop, current_idv) == 1
                && current_idv->alive == 1) {

                cub_die = rbinom(1,(pop_density*STARVATION_2YR));

                if (cub_die == 1) {
                    current_idv->alive = 0;
                }
            }

            // Infanticide if born into a pride where a member of the coalition is not the father

            if (current_idv->alive == 1
                && current_idv->age == 1
                && current_idv->my_pride->the_coalition != NULL
                && current_idv->my_pride->the_coalition->age_resident < 5) {

                current_idv -> alive = 0;
                cub_killed = 1;

            }

        }

        if (current_idv->alive == 0) {

            if (rdie == 0) {
                individual_update_events(current_idv, the_month, DIED_BACKGROUND);
            } else if (cub_die == 1) {
                individual_update_events(current_idv, the_month, DIED_STARVATION);
            } else if (cub_killed == 1) {
                individual_update_events(current_idv, the_month, DIED_INFANTICIDE_LATE);
            } else {
                individual_update_events(current_idv, the_month, DIED);
            }

        } else {

            individual_update_events(current_idv, the_month, ALIVE);

        }

        current_idv = current_idv->next;

    }

}

void individuals_die_inoldprides(t_population *pop, long the_month) { // removes prides that have been vagrants more than AGE_MAX_PRIDE_VAGRANT

    t_pride *current_pride = pop->all_prides;
    t_individual *current_indiv = NULL;

    while (current_pride != NULL) {

        if ((current_pride->stage == VAGRANT) & (current_pride->age_vagrant == AGE_MAX_PRIDE_VAGRANT)) {

            // Kill members

            for (int j = 0; j < current_pride->all_members->len; ) {

                current_indiv = g_ptr_array_index(current_pride->all_members, j);
                individual_leaves_pride(current_indiv, current_pride);
                current_indiv->alive = 0;
                individual_update_events(current_indiv, the_month, DIED);

            }

        }

        current_pride = current_pride->next;

    }

}

void individuals_die_inoldcoalitions(t_population *pop, long the_month) { // removes coalitions that have been vagrants more than AGE_MAX_COALI_VAGRANT

    t_coalition *current_coali = pop->all_coalitions;
    t_individual *current_indiv = NULL;

    while (current_coali != NULL) {

        if ((current_coali->stage == VAGRANT) & (current_coali->age_vagrant == AGE_MAX_COALI_VAGRANT)) {

            // Kill members

            for (int j = 0; j < current_coali->all_members->len; ) { // iterating and changing a collection at the same time!!!

                current_indiv = g_ptr_array_index(current_coali->all_members, j);
                current_indiv->alive = 0;
                individual_leaves_coalition(current_indiv, current_coali);
                individual_update_events(current_indiv, the_month, DIED);

            }

        }

        current_coali = current_coali->next;

    }

}

void individual_update_events(t_individual *the_idv, long the_month, int event) {

    the_idv->events[the_month] = event;

    if (event < 0) {

        t_history *new_hy = malloc(sizeof(t_history));
        new_hy->events_individual = malloc((12*R_number_of_years+1)*sizeof(int));
        for (int i = 0; i < (12*R_number_of_years+1); i++) {
            new_hy->events_individual[i] = the_idv->events[i];
        }

        new_hy->next = stats->history_individuals;
        stats->history_individuals = new_hy;
    }

}

void individuals_remove(t_population *pop) {

    t_individual *current_idv = pop->all_indiv;
    t_individual *next_idv = NULL;

    while (current_idv != NULL) {

        if (current_idv->alive == 0) {

            next_idv = current_idv->next;

            // CUB LEAVING MOTHER
            mother_loses_cub(current_idv);

            // MOTHER LEAVING CUBS
            cubs_lose_mother(current_idv);

            // LEAVE PRIDE
            individual_leaves_pride(current_idv, current_idv->my_pride);

            // LEAVE COALITION
            individual_leaves_coalition(current_idv, current_idv->my_coalition);

            // LEAVE FROM POPULATION AND TELL POPULATION IT LEFT

            // If the individual is at both the beginning and the end of the list
            if ( (current_idv->previous == NULL) & (current_idv->next == NULL) ) {
                pop->all_indiv = NULL;
            }

            // If the individual is at the beginning of the list
            else if ( (current_idv->previous == NULL) & (current_idv->next != NULL) ) {
                current_idv->next->previous = NULL; // The next individual it has no previous (it is the beginning))
                pop->all_indiv = current_idv->next; // Update the start of the list
            }

            // If the individual is at the end of the list
            else if ( (current_idv->previous != NULL) & (current_idv->next == NULL) ) {
                current_idv->previous->next = NULL; // Connect from previous to next
            }

            // Otherwise

            else {
                current_idv->next->previous = current_idv->previous; // Connect from next to previous
                current_idv->previous->next = current_idv->next; // Connect from previous to next
            }

            // FREE
            individual_free(current_idv);

            // Decrease population size
            pop->number_indiv--;

            current_idv = next_idv;

        } else { // Did not die

            current_idv = current_idv->next;

        }

    }

}

int get_individual_edgedrisk(t_population *pop, t_individual *the_idv) {

    double edge_risk = 0;
    double edge_risk_total = 0;

    // For a female
    if ( the_idv->my_pride != NULL ) {

        // in a vagrant pride
        if ( the_idv->my_pride->stage == VAGRANT) {

            edge_risk = rbinom(1, (double)pop->K_edged/pop->K_prides); // The proportion of the reserve at the edge

        }

        // in a resident pride
        if ( the_idv->my_pride->stage == RESIDENT) {

            if ( the_idv->my_pride->is_edged == EDGED) {
                edge_risk = 1;
            } else {
                edge_risk = 0;
            }

        }

    }

    // For a male
    else if ( the_idv->my_coalition != NULL ) {

        // in a vagrant coalition
        if ( the_idv->my_coalition->stage == VAGRANT) {

            edge_risk = rbinom(1, (double)pop->K_edged/pop->K_prides);
        }

        // in a resident coalition
        if ( the_idv->my_coalition->stage == RESIDENT) {

            // Look at edged status of prides belonging to his coalition
            for (int i = 0; i < the_idv->my_coalition->the_prides->len; i++) {
                edge_risk_total += ((t_pride*)g_ptr_array_index(the_idv->my_coalition->the_prides, i))->is_edged;
            }

            if ((edge_risk_total/the_idv->my_coalition->the_prides->len) > 0.5){
                edge_risk = rbinom(1, edge_risk_total/the_idv->my_coalition->the_prides->len);
            } else {
                edge_risk = 0;
            }

        }

    }

    return edge_risk;

}

void individuals_hunting(t_population *pop, long the_month) {

    // Put individuals to hunt in an array

    GPtrArray *array_hunted_individuals = g_ptr_array_sized_new(pop->number_indiv);
    t_individual *current_idv;

    int pos;
    double kill_conflict;
    double kill_hunt;
    int quota_effective;
    int quota_effective_proportion;

    /////////////////////////////////////////////////////////////////////////////////
    // CONFLICT MALES
    /////////////////////////////////////////////////////////////////////////////////

    kill_conflict = R_conflict_mortality[the_month][MALE];

    if (kill_conflict > 0) {

        current_idv = pop->all_indiv;

        while (current_idv != NULL) {

            if ( (current_idv->sex == MALE)
                & (current_idv->age >= R_conflict_age_male)
                & (current_idv->alive == 1)
                & (get_individual_edgedrisk(pop, current_idv) == 1)) {
                g_ptr_array_add(array_hunted_individuals, current_idv);
            }

            current_idv = current_idv->next;

        }

        quota_effective_proportion = rbinom(array_hunted_individuals->len, kill_conflict/100);

        while (quota_effective_proportion > 0) {

            pos = (int)runif(0, array_hunted_individuals->len-1);
            current_idv = (t_individual*)g_ptr_array_index(array_hunted_individuals, pos);
            current_idv->alive = 0;
            individual_update_events(current_idv, the_month, DIED_CONFLICT);

            // Hunted individuals no longer available
            g_ptr_array_remove_index_fast(array_hunted_individuals, pos);
            quota_effective_proportion--;

        }

    }

    g_ptr_array_empty(array_hunted_individuals);

    /////////////////////////////////////////////////////////////////////////////////
    // CONFLICT FEMALES
    /////////////////////////////////////////////////////////////////////////////////

    kill_conflict = R_conflict_mortality[the_month][FEMALE];

    if (kill_conflict > 0) {

        current_idv = pop->all_indiv;

        while (current_idv != NULL) {

            if ( (current_idv->sex == FEMALE)
                & (current_idv->age >= R_conflict_age_female)
                & (current_idv->alive == 1)
                & (get_individual_edgedrisk(pop, current_idv) == 1)) {
                g_ptr_array_add(array_hunted_individuals, current_idv);
            }

            current_idv = current_idv->next;

        }

        quota_effective_proportion = rbinom(array_hunted_individuals->len, kill_conflict/100);

        while (quota_effective_proportion > 0) {

            pos = (int)runif(0, array_hunted_individuals->len-1);
            current_idv = (t_individual*)g_ptr_array_index(array_hunted_individuals, pos);
            current_idv->alive = 0;
            individual_update_events(current_idv, the_month, DIED_CONFLICT);

            // Hunted individuals no longer available
            g_ptr_array_remove_index_fast(array_hunted_individuals, pos);
            quota_effective_proportion--;

        }

    }

    g_ptr_array_empty(array_hunted_individuals);

    /////////////////////////////////////////////////////////////////////////////////
    // HUNTING MALES
    /////////////////////////////////////////////////////////////////////////////////

    kill_hunt = R_hunting_mortality[the_month][MALE];
    int kill_hunt_int;
    int kill_hunt_total;

    if ( pop->month % 12 == 0 ) {
        pop->total_failed_hunts = 0;
    }

    pop->total_failed_hunts = 0; // to reset every month

    if (kill_hunt > 0) {

        if (R_hunter_error == 0){

            current_idv = pop->all_indiv;

            while (current_idv != NULL) {

                if ( (current_idv->sex == MALE)
                    & (current_idv->age >= R_hunting_age_male)
                    & (current_idv->alive == 1)
                    & (get_individual_edgedrisk(pop, current_idv) == 1)) {
                    g_ptr_array_add(array_hunted_individuals, current_idv);
                }

                current_idv = current_idv->next;

            }

        } else {

            current_idv = pop->all_indiv;

            while (current_idv != NULL) {

                if ( (current_idv->sex == MALE)
                    & (current_idv->estimated_age >= R_hunting_age_male)
                    & (current_idv->alive == 1)
                    & (get_individual_edgedrisk(pop, current_idv) == 1)) {
                    g_ptr_array_add(array_hunted_individuals, current_idv);
                }

                current_idv = current_idv->next;

            }

        }

        kill_hunt_int = (int)kill_hunt;
        kill_hunt_total = kill_hunt_int + rbinom(1, (double)kill_hunt - kill_hunt_int);

        quota_effective = fmin2(kill_hunt_total, array_hunted_individuals->len);

        pop->total_failed_hunts = ((double)(kill_hunt_total-quota_effective)/kill_hunt_total)*100; // counts how many hunts failed due to lack of available animals

        while (quota_effective > 0) {

            pos = (int)runif(0, array_hunted_individuals->len-1);
            current_idv = (t_individual*)g_ptr_array_index(array_hunted_individuals, pos);
            current_idv->alive = 0;
            individual_update_events(current_idv, the_month, DIED_HUNTING);

            // Hunted individuals no longer available
            g_ptr_array_remove_index_fast(array_hunted_individuals, pos);
            quota_effective--;

        }

    }

    g_ptr_array_empty(array_hunted_individuals);

    /////////////////////////////////////////////////////////////////////////////////
    // Hunting ADULT FEMALES
    /////////////////////////////////////////////////////////////////////////////////

    kill_hunt = R_hunting_mortality[the_month][FEMALE];

    if (kill_hunt > 0) {

        current_idv = pop->all_indiv;

        while (current_idv != NULL) {

            if ( (current_idv->sex == FEMALE)
                & (current_idv->age >= R_hunting_age_female)
                & (current_idv->alive == 1)
                & (get_individual_edgedrisk(pop, current_idv) == 1)) {
                g_ptr_array_add(array_hunted_individuals, current_idv);
            }

            current_idv = current_idv->next;

        }

        kill_hunt_int = (int)kill_hunt;
        kill_hunt_total = kill_hunt_int + rbinom(1, (double)kill_hunt - kill_hunt_int);

        quota_effective = fmin2(kill_hunt_total, array_hunted_individuals->len);

        while (quota_effective > 0) {

            pos = (int)runif(0, array_hunted_individuals->len-1);
            current_idv = (t_individual*)g_ptr_array_index(array_hunted_individuals, pos);
            current_idv->alive = 0;
            individual_update_events(current_idv, the_month, DIED_HUNTING);

            // Hunted individuals no longer available
            g_ptr_array_remove_index_fast(array_hunted_individuals, pos);
            quota_effective--;

        }

    }

    // Finishing
    g_ptr_array_free(array_hunted_individuals);
    individuals_remove(pop);
    prides_remove(pop);

}

void individuals_age(t_population *pop) {

    t_individual *current_idv = pop->all_indiv;
    int pos;

    pop->month++;

    while (current_idv != NULL) {

        current_idv->age++;

        // Become independent from mother
        if (current_idv->age == AGE_DISPERSAL_F_END) {//

            mother_loses_cub(current_idv);

        }

        // Counts months after a female has mated so that birth occurs at month 4.
        if (current_idv->mated == MATED){
            current_idv->age_mated++;
        }

        // Counts months after a female has given birth to calculate interbirth interval.
        if (current_idv->is_mother == IS_MOTHER){
            current_idv->age_mother++;
        }

        // hunter estimated ages

        if (current_idv->age <= 35
            && current_idv->age >= 12){

            pos = (int)runif(1, 100);
            if (pos <= 2){
                current_idv->estimated_age = (int)runif(60,83);
            }
            if ((pos > 2) & (pos <= 19)){
                current_idv->estimated_age = (int)runif(36,59);
            }
            if ((pos > 19) & (pos <= 100)){
                current_idv->estimated_age = (int)runif(12,35);
            }
        }

        if (current_idv->age <= 59
            && current_idv->age >= 36){

            pos = (int)runif(1, 100);
            if (pos <= 5){
                current_idv->estimated_age = (int)runif(84,107);
            }
            if ((pos > 5) & (pos <= 23)){
                current_idv->estimated_age = (int)runif(60,83);
            }
            if ((pos > 23) & (pos <= 88)){
                current_idv->estimated_age = (int)runif(36,59);
            }
            if ((pos > 88) & (pos <= 100)){
                current_idv->estimated_age = (int)runif(12,35);
            }
        }

        if (current_idv->age <= 83
            && current_idv->age >= 60){

            pos = (int)runif(1, 100);
            if (pos <= 2){
                current_idv->estimated_age = (int)runif(12,35);
            }
            if ((pos > 2) & (pos <= 16)){
                current_idv->estimated_age = (int)runif(36,59);
            }
            if ((pos > 16) & (pos <= 68)){
                current_idv->estimated_age = (int)runif(60,83);
            }
            if ((pos > 68) & (pos <= 100)){
                current_idv->estimated_age = (int)runif(84,107);
            }
        }

        if (current_idv->age <= 107
            && current_idv->age >= 84){

            pos = (int)runif(1, 100);
            if (pos <= 5){
                current_idv->estimated_age = (int)runif(36,59);
            }
            if ((pos > 5) & (pos <= 33)){
                current_idv->estimated_age = (int)runif(60,83);
            }
            if ((pos > 33) & (pos <= 100)){
                current_idv->estimated_age = (int)runif(84,107);
            }
        }

        current_idv = current_idv->next;

    }

}

void individuals_disperse(t_population *pop, long the_month) {

    t_pride *current_pride = pop->all_prides;
    t_individual *current_idv = NULL;

    int number_adult_F;
    int age_oldest_young_male;

    int sisters_leave;
    int brothers_leave;

    while (current_pride != NULL) {

        number_adult_F = 0;
        age_oldest_young_male = 0;

        sisters_leave = 0;
        brothers_leave = 0;

        // Calculate number of adult F and age of oldest young male

        for (int i = 0; i < current_pride->all_members->len; i++) {

            current_idv = g_ptr_array_index(current_pride->all_members, i);

            if (current_idv->sex == FEMALE && current_idv->age >= AGE_DISPERSAL_F_BEGIN ) {
                number_adult_F++;
            }

            if (current_idv->sex == MALE) {
                age_oldest_young_male = fmax2(current_idv->age, age_oldest_young_male);
            }

        }

        if (number_adult_F > MAX_PRIDE_SIZE){
            sisters_leave = 1;
            brothers_leave = 1;
        }

        if (age_oldest_young_male >= AGE_DISPERSAL_M_END) {
            brothers_leave = 1;
        }

        // Female dispersal

        if (sisters_leave == 1) {

            t_pride *new_pride = create_pride(pop, VAGRANT);

            for (int i = 0; i < current_pride->all_members->len; i++) {

                current_idv = g_ptr_array_index(current_pride->all_members, i);

                if (current_idv->sex == FEMALE
                    && current_idv->age >= AGE_DISPERSAL_F_BEGIN
                    && current_idv->age < AGE_DISPERSAL_F_END) {

                    individual_leaves_pride(current_idv, current_pride);
                    individual_joins_pride(current_idv, new_pride);
                    individual_update_events(current_idv, the_month, DISPERSED);
                    current_idv->dispersed = DISPERSE;
                }

            }

        }

        // Male dispersal

        if (brothers_leave == 1) {

            t_coalition *new_coali = create_coalition(pop, VAGRANT);

            for (int i = 0; i < current_pride->all_members->len; i++) {

                current_idv = g_ptr_array_index(current_pride->all_members, i);

                if (current_idv->sex == MALE
                    && current_idv->age >= AGE_DISPERSAL_M_BEGIN ) {

                    individual_leaves_pride(current_idv, current_pride);
                    individual_joins_coalition(current_idv, new_coali);
                    individual_update_events(current_idv, the_month, DISPERSED);
                    current_idv->dispersed = DISPERSE;
                }

            }

        }

        current_pride = current_pride->next;

    }

}

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : PRIDE AND COALITION EVENTS
/////////////////////////////////////////////////////////////////////////////////

void prides_remove(t_population *pop) {

    t_pride *current_pride = pop->all_prides;

    while (current_pride != NULL) {

        if (current_pride->all_members->len == 0) {

            current_pride = pride_leaves_pop(pop, current_pride);

        }

        else {

            current_pride = current_pride->next;

        }

    }

}

void coalitions_remove(t_population *pop) {

    t_coalition *current_coali = pop->all_coalitions;

    while (current_coali != NULL) {

        if (current_coali->all_members->len == 0) {

            current_coali = coalition_leaves_pop(pop, current_coali);

        }

        else {

            current_coali = current_coali->next;

        }

    }

}

void coalitions_age(t_population *pop) {

    t_coalition *current_coali = pop->all_coalitions;

    while (current_coali != NULL) {

        if (current_coali->stage == RESIDENT) {
            current_coali->age_resident++;
        }
        if (current_coali->stage == VAGRANT) {
            current_coali->age_vagrant++;
        }

        current_coali = current_coali->next;

    }

}

void prides_age(t_population *pop) {

    t_pride *current_pride = pop->all_prides;

    while (current_pride != NULL) {

        if (current_pride->stage == RESIDENT) {
            current_pride->age_resident++;
        }
        if (current_pride->stage == VAGRANT) {
            current_pride->age_vagrant++;
        }

        current_pride = current_pride->next;

    }

}

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : PRIDE EVENTS
/////////////////////////////////////////////////////////////////////////////////

void prides_settle(t_population *pop, long the_month) {

    int route;
    double core_density;
    double edge_density;
    double core_prides;
    double core_capacity;
    //    double prob;

    double space;
    double p;
    int find;

    t_pride *current_pride = pop->all_prides;
    GPtrArray *array_vagrant_prides = g_ptr_array_sized_new(pop->number_prides);
    t_individual *current_idv;

    while (current_pride != NULL){

        if (current_pride->stage == VAGRANT){
            g_ptr_array_add(array_vagrant_prides, current_pride); // free this array at the end
        }

        current_pride = current_pride->next;

    }

    for (int i = 0; i < array_vagrant_prides->len; i++) {

        current_pride = (t_pride*)g_ptr_array_index(array_vagrant_prides, i);

        space = (double)pop->K_prides - pop->number_prides_settled;
        p = (double)space/(double)(pop->K_prides/NUMBER_SPACES_SEARCHED);

        find = rbinom(1,p);

        int killed;

        if ((find == 1) & (pop->number_prides_settled < pop->K_prides)){

            current_pride->stage = RESIDENT;

            core_prides = pop->number_prides_settled - pop->number_prides_edged;
            core_capacity = pop->K_prides - pop->K_edged;
            core_density = (double)core_prides/core_capacity;
            edge_density = (double)pop->number_prides_edged/pop->K_edged;

            if (core_density == edge_density){
                route = rbinom(1,0.5);
            } else {
                if (core_density > edge_density){
                    route = 0;
                } else {
                    route = 1;
                }
            }

            if ((route == 0)
                & (pop->number_prides_edged < pop->K_edged)) {
                current_pride->is_edged = EDGED;
                pop->number_prides_edged++;
            } else {
                if (pop->number_prides_settled - pop->number_prides_edged < pop->K_prides - pop->K_edged) {
                    current_pride->is_edged = NOT_EDGED;
                } else {
                    current_pride->is_edged = EDGED;
                    pop->number_prides_edged++;
                }
            }

            pop->number_prides_settled++;

        } else {
            for (int j = 0; j < current_pride->all_members->len; j++) {
                current_idv = (t_individual*)g_ptr_array_index(current_pride->all_members, j);

                killed = rbinom(1,MORTALITY_NOTSETTLE);

                if (killed == 1){
                    current_idv->alive = 0;
                    individual_leaves_pride(current_idv, current_pride);
                    individual_update_events(current_idv, the_month, DIED);
                }
            }
            if (current_pride->all_members->len == 0) {
                pride_leaves_pop(pop, current_pride);
            }
        }
    }

    individuals_remove(pop);
    g_ptr_array_free(array_vagrant_prides);

}

void prides_reproduce(t_population *pop, long the_month) {

    t_pride *current_pride = pop->all_prides;
    t_individual *a_mother;
    GPtrArray *array_mothers = g_ptr_array_sized_new(MAX_GROUP_SIZE);

    int f;
    int birth_interval = 0;
    int num_mothers = 0;
    int ibi = 0;

    while (current_pride != NULL) {

        g_ptr_array_empty(array_mothers);

        if ( current_pride->the_coalition != NULL ) {

            for (int i = 0; i < current_pride->all_members->len; i++) {

                a_mother = ((t_individual*)g_ptr_array_index(current_pride->all_members, i));

                if (a_mother->age >= AGE_1ST_REPRODUCTION
                    && a_mother->sex == FEMALE
                    && a_mother->litter->len == 0) {

                    g_ptr_array_add(array_mothers, a_mother);

                }

            }

            int conception;

            for (int i = 0; i < array_mothers->len; i++) {

                a_mother = ((t_individual*)g_ptr_array_index(array_mothers, i));

                conception = rbinom(1,MATING_SUCCESS);

                if (a_mother->mated == NOT_MATED
                    && conception == 1){

                    a_mother->mated = MATED;
                    a_mother->age_mated = 0;


                }

                // If mother has mated and time since mating == 4, then mother will give birth.

                if (a_mother->mated == MATED
                    && a_mother->age_mated == 4){

                    int size = 0; // litter size
                    int draws = 1; // number of draws
                    int classes = 5; // number of classes in the multinomial
                    int vals[classes]; // array to store draws
                    rmultinom(draws, R_litter_distribution, classes, vals);
                    for(int j=0; j < classes; j++) {
                        if (vals[j] != 0) {
                            size = j + 1;
                        }
                    }
                    f = size;

                    if (f > 0) pop->live_stats[IDX_LITTERS]++;

                    for (int l = 0; l < f; l++) {

                        t_individual *new_idv = create_individual(pop,
                                                                  (rbinom(1, SEX_RATIO) == 1)?MALE:FEMALE,
                                                                  0,
                                                                  CUB);
                        individual_update_events(new_idv, the_month, BORN);
                        individual_update_events(a_mother, the_month, PARENT);
                        cub_gets_mother(new_idv, a_mother);
                        individual_joins_pride(new_idv, current_pride);

                    }

                    a_mother -> mated = NOT_MATED;

                    if (a_mother->is_mother == IS_MOTHER){
                        birth_interval += a_mother->age_mother;
                        num_mothers++;
                        a_mother->age_mother = 0;
                    } else {
                        a_mother->is_mother = IS_MOTHER;
                    }

                }
            }

        }

        current_pride = current_pride->next;

    }

    g_ptr_array_free(array_mothers);

    if (num_mothers > 0){
        ibi = birth_interval/num_mothers;
    }

    pop->live_stats[IDX_IBI] = ibi;

}



/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : COALITION EVENTS
/////////////////////////////////////////////////////////////////////////////////

void coalitions_split(t_population *pop) {

    t_coalition *current_coali = pop->all_coalitions;
    t_individual *current_idv = NULL;
    int new_size;

    while (current_coali != NULL) {

        if (current_coali->stage == VAGRANT
            && current_coali->all_members->len > COALI_MAX_SIZE ) {

            g_ptr_array_sort (current_coali->all_members, (GCompareFunc) compare_age);

            new_size = (int)runif(1, COALI_MAX_SIZE);

            t_coalition *new_coali = create_coalition(pop, VAGRANT);

            for (int i = 0; i < new_size; i++) {

                current_idv = g_ptr_array_index(current_coali->all_members, i);

                individual_leaves_coalition(current_idv, current_coali);
                individual_joins_coalition(current_idv, new_coali);

            }

        }

        current_coali = current_coali->next;

    }

}

void coalitions_merge(t_population *pop) {

    // Put coalitions to sample in arrays

    GPtrArray *array_indiv_coali1resident = g_ptr_array_sized_new(pop->number_coalitions);
    GPtrArray *array_indiv_coali1vagrant = g_ptr_array_sized_new(pop->number_coalitions);

    t_coalition *current_coali = pop->all_coalitions;

    t_individual *current_idv1;
    t_individual *current_idv2;
    t_coalition *old_coali;

    int pos1, pos2;

    // COLLECT COALITIONS OF 1

    while (current_coali != NULL) {

        if (current_coali->all_members->len == 1) {

            if (current_coali->stage == RESIDENT) {
                g_ptr_array_add(array_indiv_coali1resident, g_ptr_array_index(current_coali->all_members, 0));
            }

            if (current_coali->stage == VAGRANT) {
                g_ptr_array_add(array_indiv_coali1vagrant, g_ptr_array_index(current_coali->all_members, 0));
            }

        }

        current_coali = current_coali->next;

    }

    // VAGRANT JOINS VAGRANT

    int new_coali_vagrant_vagrant = (int)(array_indiv_coali1vagrant->len)/2;

    while (new_coali_vagrant_vagrant > 0) {

        pos1 = (int)runif(0, array_indiv_coali1vagrant->len-1);

        current_idv1 = (t_individual*)g_ptr_array_index(array_indiv_coali1vagrant, pos1);
        g_ptr_array_remove_index_fast(array_indiv_coali1vagrant, pos1);

        pos2 = (int)runif(0, array_indiv_coali1vagrant->len-1);
        current_idv2 = (t_individual*)g_ptr_array_index(array_indiv_coali1vagrant, pos2);
        g_ptr_array_remove_index_fast(array_indiv_coali1vagrant, pos2);

        old_coali = current_idv1->my_coalition;

        individual_leaves_coalition(current_idv1, current_idv1->my_coalition);
        individual_joins_coalition(current_idv1, current_idv2->my_coalition);

        coalition_leaves_pop(pop, old_coali);

        new_coali_vagrant_vagrant--;

    }

    g_ptr_array_free(array_indiv_coali1vagrant);
    g_ptr_array_free(array_indiv_coali1resident);

}

void coalitions_fight(t_population *pop, long the_month) {

    t_coalition *current_coali = pop->all_coalitions;
    t_coalition *challenging_coali;
    t_coalition *challenged_coali;

    t_pride *current_pride;

    t_individual *current_idv;

    t_individual *a_mother;

    int var  = 0;

    //    int will_challenge;
    double age_challenging;
    double age_challenged;

    //    int size_pride;
    double p;

    //    double size_ratio;
    double strength_challenged;
    double strength_challenging;
    //    double threshold;
    double win_ratio;
    double age_challenged_years;
    double age_challenging_years;
    int age_challenged_round;
    int age_challenging_round;

    double temp;

    double strength_ratio;
    double power;

    int coalition_tenure = 0;

    GPtrArray *array_challenging_coali = g_ptr_array_sized_new(pop->K_coalitions); // RISK IF VARIABLE IS 0
    GPtrArray *array_challenged_coali = g_ptr_array_sized_new(pop->K_coalitions); // number_coalitions_settled

    while (current_coali != NULL) {

        // Look for challenging coalitions
        if ((current_coali->stage == VAGRANT)
            & (current_coali->all_members->len > 0)) {
            g_ptr_array_add(array_challenging_coali, current_coali); //emptythisarray
        }

        // Look for challenged coalitions (settled with prides)
        if ((current_coali->stage == RESIDENT)
            && (current_coali->the_prides->len > 0)
            & (current_coali->all_members->len > 0)) {

            if (pop->number_coalitions_settled <= 0) {
            }

            g_ptr_array_add(array_challenged_coali, current_coali); // emptythisarray
        }

        current_coali = current_coali->next;

    }

    int fight;

    if (array_challenging_coali->len > 0 && array_challenged_coali->len > 0) {

        int pos;

        for (int i = 0; i < array_challenging_coali->len; i++) {

            if (array_challenged_coali->len == 0) {
                break;
            }

            fight = rbinom(1,FIGHT_CHANCE);

            if (fight == 1) {

                pos = (int)runif(0, array_challenged_coali->len-1);

                challenging_coali = (t_coalition*)g_ptr_array_index(array_challenging_coali, i);
                challenged_coali = (t_coalition*)g_ptr_array_index(array_challenged_coali, pos);

                // measure the srength of the challenge

                strength_challenged = 0;
                strength_challenging = 0;

                // calculate the total strength of the coalitions
                for (int j = 0; j < challenging_coali->all_members->len; j++) {

                    temp = 0;

                    age_challenging = ((t_individual*)g_ptr_array_index(challenging_coali->all_members, j))->age;
                    age_challenging_years = (double)age_challenging/12;
                    age_challenging_round = (double)age_challenging_years + 0.5;

                    if (age_challenging_round >= 6){
                        temp = (double)((age_challenging_round-6)+1);
                        temp = 101 - pow(temp,2.0);
                        strength_challenging = strength_challenging + temp;
                    } else {
                        temp = (((6-age_challenging_round)*2.25)+1);
                        temp = 101 - pow(temp, 2.0);
                        strength_challenging = strength_challenging + temp;
                    }

                }

                for (int j = 0; j < challenged_coali->all_members->len; j++) {

                    temp = 0;

                    age_challenged = ((t_individual*)g_ptr_array_index(challenged_coali->all_members, j))->age;
                    age_challenged_years = (double)age_challenged/12;
                    age_challenged_round = (double)age_challenged_years + 0.5;

                    if (age_challenged_round >= 6){
                        temp = (double)((age_challenged_round-6)+1);
                        temp = (101 - pow(temp,2.0))*RESIDENT_ADVANTAGE;
                        strength_challenged = strength_challenged + temp;
                    } else {
                        temp = (((6-age_challenged_round)*2.25)+1);
                        temp = (101 - pow(temp, 2.0))*RESIDENT_ADVANTAGE;
                        strength_challenged = strength_challenged + temp;
                    }

                }

                // Success of challenge

                strength_ratio = (double)strength_challenged/strength_challenging;
                power = pow((challenged_coali->all_members->len - challenging_coali->all_members->len),2.0)+2.0;

                win_ratio = pow(strength_ratio, power);

                p = (double)win_ratio/(win_ratio+1);

                int killed;
                int kill;

                if ( rbinom(1, p) == 0 ) { // Outcome of the challenge - chance of pride males losing the challenge

                    var++;
                    int pos_pride;

                    // For all prides from this coalition - coalition only loses one of its prides

                    pos_pride = (int)runif(0, challenged_coali->the_prides->len-1);

                    current_pride = (t_pride*)g_ptr_array_index(challenged_coali->the_prides, pos_pride);


                    pride_changes_coalition(current_pride, challenging_coali);

                    // Infanticide

                    int x;
                    double prob;

                    for (int l = 0; l < current_pride->all_members->len; l++) {
                        current_idv = (t_individual*)g_ptr_array_index(current_pride->all_members, l);

                        x = current_idv->age;
                        prob = (0.0000029539*pow(x,6.0)) - (0.000334047*pow(x,5.0)) + (0.0137274757*pow(x,4.0)) - (0.2391477426*pow(x,3.0)) + (1.4663687561*pow(x,2.0)) - (3.5297032314*x) + 100;
                        prob = prob/100;

                        kill = rbinom(1,prob);

                        if (kill == 1
                            && current_idv->age < AGE_SURVIVE_INFANTICIDE) {
                            current_idv->alive = 0;
                            individual_update_events(current_idv, the_month, DIED_INFANTICIDE_FIGHT);
                        } else {
                            if (current_idv->sex == MALE
                                && current_idv->age < AGE_SURVIVE_INFANTICIDE){
                                current_idv->alive = 0;
                                individual_update_events(current_idv, the_month, DIED_INFANTICIDE_FIGHT);
                            }
                        }

                    }

                    // Mothers have x% chance of being killed in a fight

                    int killed;
                    int adult_females = 0;

                    for (int z = 0; z < current_pride->all_members->len; z++) {

                        a_mother = ((t_individual*)g_ptr_array_index(current_pride->all_members, z));

                        if (a_mother->age >= AGE_DISPERSAL_F_END && a_mother->sex == FEMALE){
                            adult_females++;
                        }

                        if (a_mother->age >= AGE_1ST_REPRODUCTION && a_mother->sex == FEMALE && a_mother->litter->len > 0) { // added female sex?

                            killed = rbinom(1,MORTALITY_FIGHT_F);

                            if(killed == 1){

                                a_mother->alive = 0;
                                individual_update_events(a_mother, the_month, KILLED_FIGHT_F);

                            }

                        }

                    }

                    // FEMALE individuals disperse and create new pride

                    t_pride *new_pride = create_pride(pop, VAGRANT);

                    for (int i = 0; i < current_pride->all_members->len; i++) {

                        current_idv = g_ptr_array_index(current_pride->all_members, i);

                        if ( current_idv->sex == FEMALE
                            && current_idv->age >= AGE_DISPERSAL_F_BEGIN_TRIGGER
                            && current_idv->age < AGE_DISPERSAL_F_END_TRIGGER) { // AGE_DISPERSAL_F_END

                            individual_leaves_pride(current_idv, current_pride);
                            individual_joins_pride(current_idv, new_pride);
                            individual_update_events(current_idv, the_month, DISPERSED);
                            current_idv->dispersed = DISPERSE;
                        }

                    }

                    if (new_pride->all_members->len == 0) {

                        pride_leaves_pop(pop, new_pride);

                    }

                    // MALE individuals disperse and create new coalition

                    t_coalition *new_coali = create_coalition(pop, VAGRANT);

                    for (int i = 0; i < current_pride->all_members->len; i++) {

                        current_idv = g_ptr_array_index(current_pride->all_members, i);

                        if (current_idv->sex == MALE && current_idv->age >= AGE_SURVIVE_INFANTICIDE) { // AGE_SURVIVE_INFANTICIDE

                            individual_leaves_pride(current_idv, current_pride);
                            individual_joins_coalition(current_idv, new_coali);
                            individual_update_events(current_idv, the_month, DISPERSED);
                            current_idv->dispersed = DISPERSE;
                        }

                    }

                    if (new_coali->all_members->len == 0) {

                        coalition_leaves_pop(pop, new_coali);

                    }


                    //} // part of initial for loop

                    // challenged_coali no longer available

                    g_ptr_array_remove_index_fast(array_challenged_coali, pos);

                    coalition_tenure += challenged_coali->age_resident;

                    challenging_coali->stage = RESIDENT;
                    challenging_coali->age_resident = 0;
                    challenging_coali->age_vagrant = 0;
                    pop->number_coalitions_settled++;

                    // Challenged coalition only becomes vagrant if it doesn't hold other prides

                    if(challenged_coali->the_prides->len == 0){
                        challenged_coali->stage = VAGRANT;
                        challenged_coali->age_resident = 0;
                        challenged_coali->age_vagrant = 0;
                        pop->number_coalitions_settled--;
                    }

                    // individuals in defeated coalition have x% chance of dying if there is a physical encounter (determined by range of probabilities)

                    if (p >= 0.15){



                        for (int i = 0; i < challenged_coali->all_members->len; i++) {
                            current_idv = (t_individual*)g_ptr_array_index(challenged_coali->all_members, i);

                            killed = rbinom(1,MORTALITY_FIGHT_M);

                            if (killed == 1){
                                current_idv->alive = 0;
                                individual_leaves_coalition(current_idv, challenged_coali);
                                individual_update_events(current_idv, the_month, KILLED_FIGHT);
                            }
                        }
                        if (challenged_coali->all_members->len == 0) {
                            coalition_leaves_pop(pop, challenged_coali);
                        }
                    }

                } else {

                    if (p <= 0.7){

                        for (int i = 0; i < challenging_coali->all_members->len; i++) {
                            current_idv = (t_individual*)g_ptr_array_index(challenging_coali->all_members, i);

                            killed = rbinom(1,MORTALITY_FIGHT_M);

                            if (killed == 1){
                                current_idv->alive = 0;
                                individual_leaves_coalition(current_idv, challenging_coali);
                                individual_update_events(current_idv, the_month, KILLED_FIGHT);
                            }
                        }
                        if (challenging_coali->all_members->len == 0) {
                            coalition_leaves_pop(pop, challenging_coali);

                        }

                    }

                }

            }

        }

    }

    g_ptr_array_free(array_challenged_coali);
    g_ptr_array_free(array_challenging_coali);

    pop->live_stats[IDX_TAKEOVERS] = var;
    pop->live_stats[IDX_COALITION_TENURE] = coalition_tenure;

    individuals_remove(pop);
}


void coalitions_meet_prides(t_population *pop, long the_month) {

    t_individual *current_idv;
    t_individual *a_mother;

    int pos_pride;

    // Look for prides with no coalition

    GPtrArray *array_pride_no_coali = g_ptr_array_sized_new(pop->number_prides);
    t_pride *current_pride = pop->all_prides;

    while (current_pride != NULL) {

        if (current_pride->the_coalition == NULL && current_pride->stage == RESIDENT) {
            g_ptr_array_add(array_pride_no_coali, current_pride);
        }

        current_pride = current_pride->next;
    }

    g_ptr_array_sort (array_pride_no_coali, (GCompareFunc) compare_pride_size);

    // Look for coalitions with no pride

    GPtrArray *array_coali_no_pride = g_ptr_array_sized_new(pop->number_coalitions);
    t_coalition *current_coali = pop->all_coalitions;

    while (current_coali != NULL) {

        if ((current_coali->the_prides->len == 0 ||
             current_coali->stage == VAGRANT ||
             current_coali->all_members->len > 1) &&
            current_coali->all_members->len > 0) { // prides->len == 0 && current_coali->stage == VAGRANT
            g_ptr_array_add(array_coali_no_pride, current_coali);
        }

        current_coali = current_coali->next;
    }

    g_ptr_array_sort (array_coali_no_pride, (GCompareFunc) compare_nprides_btw_coalitions);

    int kill;
    double space;
    int find;
    double p;
    double p_adjusted;

    for (int i = 0; i < array_coali_no_pride->len;  i++) {

        if (array_pride_no_coali->len == 0){
            break;
        }

        current_coali = (t_coalition*)g_ptr_array_index(array_coali_no_pride, i);

        space = (double)array_pride_no_coali->len;
        p = (double)space/(double)(pop->K_prides/NUMBER_PRIDES_SEARCHED);

        p_adjusted = (double)p/(pow(current_coali->the_prides->len,2.0) + 1);

        if (p_adjusted > 1){
            p_adjusted = 1;
        }

        find = rbinom(1,p_adjusted);

        if ((find == 1)
            & (pop->number_coalitions_settled < pop->K_coalitions)) {

            pos_pride = (int)runif(0, array_pride_no_coali->len-1);
            current_pride = (t_pride*)g_ptr_array_index(array_pride_no_coali, pos_pride);

            // Establish connection
            current_pride->the_coalition = current_coali;
            g_ptr_array_add(current_coali->the_prides, current_pride);

            // Tell coalition it is now resident
            if (current_coali->stage == VAGRANT){

                current_coali->stage = RESIDENT;
                current_coali->age_resident = 0; // affects infanticide!
                current_coali->age_vagrant = 0;
                pop->number_coalitions_settled++; // affects whether more can settle

            }

            // Infanticide

            int x;
            double prob;

            for (int l = 0; l < current_pride->all_members->len; l++) {
                current_idv = (t_individual*)g_ptr_array_index(current_pride->all_members, l);

                x = current_idv->age;
                prob = (0.0000029539*pow(x,6.0)) - (0.000334047*pow(x,5.0)) + (0.0137274757*pow(x,4.0)) - (0.2391477426*pow(x,3.0)) + (1.4663687561*pow(x,2.0)) - (3.5297032314*x) + 100;
                prob = prob/100;

                kill = rbinom(1,prob);

                if (kill == 1
                    && current_idv->age < AGE_SURVIVE_INFANTICIDE) {
                    current_idv->alive = 0;
                    individual_update_events(current_idv, the_month, DIED_INFANTICIDE_FREE);
                } else {
                    if (current_idv->sex == MALE
                        && current_idv->age < AGE_SURVIVE_INFANTICIDE){
                        current_idv->alive = 0;
                        individual_update_events(current_idv, the_month, DIED_INFANTICIDE_FIGHT); // Males < AGE_SURVIVE_INFANTICIDE don't ever survive infanticide
                    }
                }
            }

            // Mothers have x% chances of being killed

            int killed;

            for (int z = 0; z < current_pride->all_members->len; z++) {

                a_mother = ((t_individual*)g_ptr_array_index(current_pride->all_members, z));

                if (a_mother->age >= AGE_1ST_REPRODUCTION && a_mother->sex == FEMALE && a_mother->litter->len > 0) { // added female sex?

                    killed = rbinom(1,MORTALITY_FIGHT_F);

                    if(killed == 1){

                        a_mother->alive = 0;
                        individual_update_events(a_mother, the_month, KILLED_FIGHT_F);

                    }

                }

            }

            // FEMALE individuals disperse and create new pride

            t_pride *new_pride = create_pride(pop, VAGRANT);

            for (int i = 0; i < current_pride->all_members->len; i++) {

                current_idv = g_ptr_array_index(current_pride->all_members, i);

                if ( current_idv->sex == FEMALE
                    && current_idv->age >= AGE_DISPERSAL_F_BEGIN_TRIGGER
                    && current_idv->age < AGE_DISPERSAL_F_END_TRIGGER) {

                    individual_leaves_pride(current_idv, current_pride);
                    individual_joins_pride(current_idv, new_pride);
                    individual_update_events(current_idv, the_month, DISPERSED);
                    current_idv->dispersed = DISPERSE;
                }

            }

            if (new_pride->all_members->len == 0) {

                pride_leaves_pop(pop, new_pride);

            }

            // MALE individuals disperse and create new coalition

            t_coalition *new_coali = create_coalition(pop, VAGRANT);

            for (int i = 0; i < current_pride->all_members->len; i++) {

                current_idv = g_ptr_array_index(current_pride->all_members, i);

                if (current_idv->sex == MALE && current_idv->age >= AGE_SURVIVE_INFANTICIDE ) { //AGE_SURVIVE_INFANTICIDE

                    individual_leaves_pride(current_idv, current_pride);
                    individual_joins_coalition(current_idv, new_coali);
                    individual_update_events(current_idv, the_month, DISPERSED);
                    current_idv->dispersed = DISPERSE;
                }

            }

            if (new_coali->all_members->len == 0) {

                coalition_leaves_pop(pop, new_coali);

            }


            g_ptr_array_remove_index_fast(array_pride_no_coali, pos_pride);

        }

    }

    g_ptr_array_free(array_pride_no_coali);
    g_ptr_array_free(array_coali_no_pride);

    individuals_remove(pop);

}


/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : MODEL
/////////////////////////////////////////////////////////////////////////////////

void cycle_year(t_population *pop, long the_run, long the_year, struct statistics *stats) {

    for (long the_month = 12*(the_year-1) + 1; the_month < 12*(the_year-1) + 1 + 12; the_month++) {

        individuals_age(pop);
        prides_age(pop);
        coalitions_age(pop);

        individuals_die(pop, the_month);
        individuals_die_inoldprides(pop, the_month);
        individuals_die_inoldcoalitions(pop, the_month);

        individuals_hunting(pop, the_month-1);
        // Take for example a 20 year simulation
        // Statistics are returned for 20*12+1 = 241 months, because we want to know the initial population
        // The hunting array has 20*12 = 240 rows, so from 0 to 239
        // Note the mismatch between the hunting array and the simulation loop
        // The first month is noted 0 and is before the loop (just initial population)
        // At month 1 (beginning of the loop), hunting array is accessed as it first row which is 0
        // Thefore, we need month-1
        // We don't need month-1 for the others because that would move back to initial population, month 0
        // The loops goes from 1 to 240 (<241)
        // At the last month of the loop, 240, we access hunting array 239

        individuals_remove(pop);
        prides_remove(pop);
        coalitions_remove(pop);

        individuals_disperse(pop, the_month);
        prides_settle(pop, the_month);

        coalitions_split(pop);
        coalitions_merge(pop);

        coalitions_fight(pop, the_month);
        coalitions_meet_prides(pop, the_month);

        prides_reproduce(pop, the_month);

        do_statistics(pop, the_run, the_month, stats);

    }


}

void do_statistics(t_population *pop, long the_run, long the_month, struct statistics *stats) {

    // Create a stat within pop to be able to track particular events like fights

    double *stat = malloc(NUMBER_OF_STATS*sizeof(double));
    for (int i = 0; i < NUMBER_OF_STATS; i++) {
        stat[i] = 0.0;
    }

    // Population size

    stat[IDX_NINDIV] = pop->number_indiv;

    stat[IDX_NPRIDES] = pop->number_prides;

    stat[IDX_NCOALIS] = pop->number_coalitions;

    stat[IDX_FAILED_HUNTS] = pop->total_failed_hunts;

    // Individuals

    t_individual *current_idv = pop->all_indiv;

    while (current_idv != NULL) {

        // N Females
        if (current_idv->sex == FEMALE) {
            stat[IDX_NFEMALES]++;
        }

        // N Males
        if (current_idv->sex == MALE) {
            stat[IDX_NMALES]++;
        }

        // N  adult Females
        if ((current_idv->sex == FEMALE)
            & (current_idv->age >= AGE_DISPERSAL_F_END)) {
            stat[IDX_NFEMALES_ADULT]++;
        }

        // N adult Males
        if ((current_idv->sex == MALE)
            & (current_idv->age >= AGE_DISPERSAL_M_END)) {
            stat[IDX_NMALES_ADULT]++;
        }

        // N adult males in core prides
        if (current_idv->sex == MALE
            && current_idv->age >= AGE_DISPERSAL_M_END
            && get_individual_edgedrisk(pop, current_idv) == 0) {
            stat[IDX_NMALES_CORE]++;
        }

        // N adult females in core prides
        if (current_idv->sex == FEMALE
            && current_idv->age >= AGE_DISPERSAL_F_END_TRIGGER
            && current_idv->my_pride != NULL
            && current_idv->my_pride->stage == RESIDENT
            && current_idv->my_pride->is_edged == NOT_EDGED){
            stat[IDX_NFEMALES_CORE]++;
        }

        // N cubs in core prides
        if (current_idv->age < AGE_DISPERSAL_F_BEGIN
            && current_idv->my_pride != NULL
            && current_idv->my_pride->stage == RESIDENT
            && current_idv->my_pride->is_edged == NOT_EDGED){
            stat[IDX_NCUBS_CORE]++;
        }

        // N adult males at the edge prides
        if (current_idv->sex == MALE
            && current_idv->age >= AGE_DISPERSAL_M_END
            && get_individual_edgedrisk(pop, current_idv) == 1) {
            stat[IDX_NMALES_EDGE]++;
        }


        //  N adult females in edge prides
        if (current_idv->sex == FEMALE
            && current_idv->age >= AGE_DISPERSAL_F_END_TRIGGER
            && current_idv->my_pride != NULL
            && current_idv->my_pride->stage == RESIDENT
            && current_idv->my_pride->is_edged == EDGED){
            stat[IDX_NFEMALES_EDGE]++;
        }

        // N cubs in edge prides
        if (current_idv->age < AGE_DISPERSAL_F_BEGIN
            && current_idv->my_pride != NULL
            && current_idv->my_pride->stage == RESIDENT
            && current_idv->my_pride->is_edged == EDGED){
            stat[IDX_NCUBS_EDGE]++;
        }

        // N Females dispersed
        if ((current_idv->sex == FEMALE)
            & (current_idv->age >= AGE_DISPERSAL_F_END)
            & (current_idv->dispersed == DISPERSE)) {
            stat[IDX_NFEMALES_DISPERSED]++;
        }

        // N Individuals CORE
        if (get_individual_edgedrisk(pop, current_idv) == 0) {
            stat[IDX_NINDIV_CORE]++;
        }

        // N Individuals EDGE
        if (get_individual_edgedrisk(pop, current_idv) == 1) {
            stat[IDX_NINDIV_EDGE]++;
        }

        // Individual age
        stat[IDX_AGE] += current_idv->age;

        current_idv = current_idv->next;
    }

    // Coalitions

    t_coalition *current_coali = pop->all_coalitions;

    while (current_coali != NULL) {

        if (current_coali->stage == RESIDENT) {
            stat[IDX_NCOALIS_RESIDENT]++;
            stat[IDX_COALISIZE_RESIDENT] += current_coali->all_members->len;
            stat[IDX_NPRIDES_PER_COALI] += current_coali->the_prides->len;
        }
        if (current_coali->stage == VAGRANT) {
            stat[IDX_NCOALIS_VAGRANT]++;
            stat[IDX_COALISIZE_VAGRANT] += current_coali->all_members->len;
        }

        current_coali = current_coali->next;

    }

    // Prides

    t_pride *current_pride = pop->all_prides;

    while (current_pride != NULL) {

        if (current_pride->stage == RESIDENT) {
            stat[IDX_NPRIDES_RESIDENT]++;
            stat[IDX_PRIDESIZE_RESIDENT] += current_pride->all_members->len;
        }
        if (current_pride->stage == VAGRANT) {
            stat[IDX_NPRIDES_VAGRANT]++;
            stat[IDX_PRIDESIZE_VAGRANT] += current_pride->all_members->len;
        }

        // Pride size of edge prides
        if ((current_pride->stage == RESIDENT)
            & (current_pride->is_edged == EDGED)) {
            stat[IDX_PRIDESIZE_EDGED] += current_pride->all_members->len;
        }

        // Pride size of core prides
        if ((current_pride->stage == RESIDENT)
            & (current_pride->is_edged == NOT_EDGED)) {
            stat[IDX_NPRIDES_CORE]++;
            stat[IDX_PRIDESIZE_CORE] += current_pride->all_members->len;
        }

        current_pride = current_pride->next;

    }

    //    Another way:
    stat[IDX_NPRIDES_EDGED] = pop->number_prides_edged;

    // Average for size stats


    if (stat[IDX_NINDIV] > 0) {
        stat[IDX_AGE] = stat[IDX_AGE] / stat[IDX_NINDIV];
    }

    if (stat[IDX_NCOALIS_RESIDENT] > 0) {
        stat[IDX_NPRIDES_PER_COALI] /= stat[IDX_NCOALIS_RESIDENT];
    }

    if (stat[IDX_NCOALIS_RESIDENT] > 0) {
        stat[IDX_COALISIZE_RESIDENT] /= stat[IDX_NCOALIS_RESIDENT];
    }

    if (stat[IDX_NCOALIS_VAGRANT] > 0) {
        stat[IDX_COALISIZE_VAGRANT] /= stat[IDX_NCOALIS_VAGRANT];
    }

    if (stat[IDX_NPRIDES_RESIDENT] > 0) {
        stat[IDX_PRIDESIZE_RESIDENT] /= stat[IDX_NPRIDES_RESIDENT];
    }

    if (stat[IDX_NPRIDES_VAGRANT] > 0) {
        stat[IDX_PRIDESIZE_VAGRANT] /= stat[IDX_NPRIDES_VAGRANT];
    }

    if (stat[IDX_NPRIDES_EDGED] > 0) {
        stat[IDX_PRIDESIZE_EDGED] /= stat[IDX_NPRIDES_EDGED];
    }

    if (stat[IDX_NPRIDES_CORE] > 0) {
        stat[IDX_PRIDESIZE_CORE] /= stat[IDX_NPRIDES_CORE];
    }

    if (stat[IDX_NPRIDES_CORE] > 0) {
        stat[IDX_NFEMALES_CORE] /= stat[IDX_NPRIDES_CORE];
    }

    if (stat[IDX_NPRIDES_CORE] > 0) {
        stat[IDX_NCUBS_CORE] /= stat[IDX_NPRIDES_CORE];
    }

    if (stat[IDX_NPRIDES_EDGED] > 0) {
        stat[IDX_NFEMALES_EDGE] /= stat[IDX_NPRIDES_EDGED];
    }

    if (stat[IDX_NPRIDES_EDGED] > 0) {
        stat[IDX_NCUBS_EDGE] /= stat[IDX_NPRIDES_EDGED];
    }

    if (pop->live_stats[IDX_TAKEOVERS] > 0) {
        pop->live_stats[IDX_COALITION_TENURE] /= pop->live_stats[IDX_TAKEOVERS];
    }

    // Get all

    for (int i = 0; i < NUMBER_OF_STATS; i++) {
        stats->runs[the_run][the_month][i] = stat[i];
    }

    // Get specific stats
    stats->runs[the_run][the_month][IDX_TAKEOVERS] = pop->live_stats[IDX_TAKEOVERS];
    stats->runs[the_run][the_month][IDX_LITTERS] = pop->live_stats[IDX_LITTERS];
    stats->runs[the_run][the_month][IDX_COALITION_TENURE] = pop->live_stats[IDX_COALITION_TENURE];
    stats->runs[the_run][the_month][IDX_IBI] = pop->live_stats[IDX_IBI];

    pop->live_stats[IDX_TAKEOVERS] = 0;
    pop->live_stats[IDX_LITTERS] = 0;
    pop->live_stats[IDX_COALITION_TENURE] = 0;
    pop->live_stats[IDX_IBI] = 0;

    free(stat);

}

void collect_events(t_population *pop, struct statistics *stats) {

    t_individual *current_idv = pop->all_indiv;

    while (current_idv != NULL) {

        t_history *new_hy = malloc(sizeof(t_history));
        new_hy->events_individual = malloc((12*R_number_of_years+1)*sizeof(int));
        for (int i = 0; i < (12*R_number_of_years+1); i++) {
            new_hy->events_individual[i] = current_idv->events[i];
        }

        new_hy->next = stats->history_individuals;
        stats->history_individuals = new_hy;

        current_idv = current_idv->next;
    }

}

/////////////////////////////////////////////////////////////////////////////////
#pragma mark FUNCTION : MEMORY
/////////////////////////////////////////////////////////////////////////////////

void individual_free(t_individual *idv) {

    g_ptr_array_free(idv->litter);
    free(idv->events);

    free(idv);

}

void free_population(t_population *pop) {

    // Free individuals

    t_individual *next_idv;

    while (pop->all_indiv != NULL) {
        next_idv = pop->all_indiv->next;
        individual_free(pop->all_indiv);
        pop->all_indiv = next_idv;
    }

    // Free prides

    t_pride *next_pride;

    while (pop->all_prides != NULL) {
        next_pride = pop->all_prides->next;
        g_ptr_array_free(pop->all_prides->all_members);
        free(pop->all_prides);
        pop->all_prides = next_pride;
    }

    // Free coalitions

    t_coalition *next_coali;

    while (pop->all_coalitions != NULL) {
        next_coali = pop->all_coalitions->next;
        g_ptr_array_free(pop->all_coalitions->all_members);
        free(pop->all_coalitions);
        pop->all_coalitions = next_coali;
    }

    free(pop->live_stats);

}
