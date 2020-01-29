//
// Created by andrey on 1/21/20.
//

#ifndef GA_CAUCHY_MUTATION_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "random_number_generator.h"

void box_muller_transform(double *uniform_vector, int NUMBER_GENES, long seed);

void cauchy_mutation(double *after_mut, double *after_cross, double *left_border, double *right_border, int NUMBER_ORGANISMS, int NUMBER_GENES);

#define GA_CAUCHY_MUTATION_H

#endif //GA_CAUCHY_MUTATION_H
