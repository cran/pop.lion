/* tools.c
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

/* This file contains modified functions from the GLIB
 * GLIB - Library of useful routines for C programming
 * Copyright (C) 1995-1997  Peter Mattis, Spencer Kimball and Josh MacDonald
 * Available at ftp://ftp.gtk.org/pub/gtk/.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <Rmath.h>

#include "tools.h"

/*******************************************************************************
 New with size
 *******************************************************************************/

GPtrArray* g_ptr_array_sized_new(int reserved_size) {
	
	GPtrArray *array = malloc(sizeof(GPtrArray));
    
	if (reserved_size > 0) {
		
		array->pdata = malloc(reserved_size * sizeof(gpointer));
		
		for (int i = 0; i < reserved_size; i++) {
			
			array->pdata[i] = NULL;
			
		}
		
	} else {
		
		array->pdata = NULL;
		
	}
    
    array->len = 0;
	
	array->alloc = reserved_size;
	
	return (GPtrArray*) array;
    
}

/*******************************************************************************
 Add
 *******************************************************************************/

void g_ptr_array_add (GPtrArray *array, gpointer data) {
	
	if (array->len == array->alloc) {
		
		array->pdata = realloc(array->pdata, 2*array->alloc*sizeof(gpointer));
		
        for (int i = array->alloc; i < 2*array->alloc; i++) {
			
			array->pdata[i] = NULL;
			
		}
        
        array->alloc = 2 * array->alloc;
		
	}
	
	array->pdata[array->len] = data;
    
    array->len++;
	
}

/*******************************************************************************
 Remove fast at index
 *******************************************************************************/

void g_ptr_array_remove_index_fast(GPtrArray *array, int index) {
		
	if (index == array->len - 1) {
		
		array->pdata[index] = NULL;
		
		array->len--;
		
	}
		
	if (index < array->len - 1) {
		
		array->pdata[index] = array->pdata[array->len-1];
		
		array->pdata[array->len-1] = NULL;
		
		array->len--;
		
	}
	
}

/*******************************************************************************
 Remove fast
 *******************************************************************************/

void g_ptr_array_remove_fast(GPtrArray *array, gpointer data) {
	
	for (int i = 0; i < array->len; i++) {
		
		if (array->pdata[i] == data) {
			
			g_ptr_array_remove_index_fast(array, i);
			
			break;
			
		}
    }
}

/*******************************************************************************
 Empty
 *******************************************************************************/

void g_ptr_array_empty(GPtrArray *array) {
    
    for (int i = 0; i < array->len; i++) {
        
        array->pdata[i] = NULL;
        
    }
    
    array->len = 0;
	
}

/*******************************************************************************
 Free
 *******************************************************************************/

void g_ptr_array_free(GPtrArray *array) {
	
	free(array->pdata);
	
	free(array);
	
}

/*******************************************************************************
 Shuffle array
 *******************************************************************************/

void g_ptr_array_shuffle(GPtrArray *array) {
	
	int i, j;
	gpointer *temp;
	
	for ( i = array->len-1; i>=0; i--) {
		
		j = runif(0, i);
		
		temp = array->pdata[j];
		
		array->pdata[j] = array->pdata[i];
		
		array->pdata[i] = temp;
		
	}
	
}

/*******************************************************************************
 Sort array
 *******************************************************************************/

void g_ptr_array_sort(GPtrArray *array, GCompareFunc compare_func){
    
    qsort(array->pdata,
           array->len,
           sizeof (gpointer),
           compare_func);
    
}

/*******************************************************************************
 Shape for beta distribution
 *******************************************************************************/

double beta_shape(double mu, double sigma) {
	
	return( fmax2(0, ( pow(mu,2) - pow(mu,3) - mu*pow(sigma, 2) ) / pow(sigma,2)) );
	
}

/*******************************************************************************
 Rate for beta distribution
 *******************************************************************************/

double beta_rate(double mu, double sigma) {
	
	return( fmax2(0, ( mu - 2*pow(mu,2) + pow(mu,3) - pow(sigma,2) + mu*pow(sigma, 2) ) / pow(sigma,2)) );
	
}

/*******************************************************************************
 Shape for gamma distribution
 *******************************************************************************/

double gamma_shape(double mu, double sigma) {
	
	double x = 0;
	
	if (sigma > 0) {
		x = pow(mu,2) / pow(sigma,2);
	}
	
	return(x);
	
}

/*******************************************************************************
 Rate for gamma distribution
 *******************************************************************************/

double gamma_rate(double mu, double sigma) {
	
	double x = 0;
	
	if (sigma > 0) {
		x = mu / pow(sigma,2);
	}
	
	return(x);
	
}

