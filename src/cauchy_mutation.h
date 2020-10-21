//
// Created by andrey on 1/21/20.
//

#ifndef GA_CAUCHY_MUTATION_H
#define GA_CAUCHY_MUTATION_H

#include <sys/time.h>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include "random_number_generator.h"

void normalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                               int number_organisms, int number_genes, int number_logscale_params,
        /*out*/   double *genes_output, double *left_border_normalized, double *right_border_normalized,
        /*other*/ double left_border_value = 0, double right_border_value = 1);

void denormalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                                 int number_organisms, int number_genes, int number_logscale_params,
        /*out*/   double *genes_output,
        /*other*/ double left_border_value = 0, double right_border_value = 1);

void transform_genes(/*in*/    double *genes_input, double *left_border_input, double *right_border_input,
                               int number_organisms, int number_genes, int number_log10scale_params,
                               const double *left_border_output, const double *right_border_output,
        /*out*/   double *genes_output);

void transform_genes_back(/*in*/    double *genes_input, double *left_border_input, double *right_border_input,
                                    int number_organisms, int number_genes, int number_log10scale_params,
                                    const double *left_border_output, const double *right_border_output,
        /*out*/   double *genes_output);

void cauchy_mutation(double *genes_output, double *genes_input, double *left_border, double *right_border,
                     int NUMBER_ORGANISMS, int NUMBER_GENES);


#endif //GA_CAUCHY_MUTATION_H
