/* tools.h
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

#ifndef TOOLS_H
#define TOOLS_H

typedef struct t_history t_history;

struct t_history {
    int *events_individual; // This is passed from the individual. Be careful to not free it when releasing the individual.
    t_history *next;
};

struct statistics {
    double ***runs;
    t_history *history_individuals;
};

double beta_shape(double mu, double sigma);
double beta_rate(double mu, double sigma);
double gamma_shape(double mu, double sigma);
double gamma_rate(double mu, double sigma);

typedef int gint;
typedef void *gpointer;
typedef const void *gconstpointer;
typedef gint (*GCompareFunc) (gconstpointer a, gconstpointer b);

typedef struct _GPtrArray GPtrArray;

struct _GPtrArray {
	gpointer *pdata;
	int len;
	int alloc;
};

#define g_ptr_array_index(array, index_) ((array)->pdata)[index_]

GPtrArray* g_ptr_array_sized_new(int reserved_size);
void g_ptr_array_add (GPtrArray *array, gpointer data);
void g_ptr_array_remove_index_fast(GPtrArray *array, int index);
void g_ptr_array_remove_fast(GPtrArray *array, gpointer data);
void g_ptr_array_empty(GPtrArray *array);
void g_ptr_array_free(GPtrArray *array);
void g_ptr_array_shuffle(GPtrArray *array);
void g_ptr_array_sort(GPtrArray *array, GCompareFunc compare_func);

#endif
