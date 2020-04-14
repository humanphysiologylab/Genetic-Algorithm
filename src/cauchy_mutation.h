//
// Created by andrey on 1/21/20.
//

#ifndef GA_CAUCHY_MUTATION_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <sys/time.h>

#include "random_number_generator.h"


void normalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                               int number_organisms, int number_genes, int number_conductancies,
                     /*out*/   double *genes_output, double *left_border_normalized, double *right_border_normalized,
                     /*other*/ double left_border_value = 0, double right_border_value = 1);

void denormalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                                 int number_organisms, int number_genes, int number_conductancies,
                       /*out*/   double *genes_output,
                       /*other*/ double left_border_value = 0, double right_border_value = 1);

void box_muller_transform(double *uniform_vector, int NUMBER_GENES, long seed);

void cauchy_mutation(double *after_mut, double *after_cross, double *left_border, double *right_border, int NUMBER_ORGANISMS, int NUMBER_GENES);

#define GA_CAUCHY_MUTATION_H

#endif //GA_CAUCHY_MUTATION_H
