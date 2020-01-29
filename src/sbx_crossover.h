//
// Created by andrey on 1/22/20.
//

#ifndef GA_SBX_CROSSOVER_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

void sbx_crossover(double *next_generation, double *after_cross, int *mpool, double *left_border, double *right_border,
                   int NUMBER_ORGANISMS, int NUMBER_GENES);

#include "random_number_generator.h"

#define GA_SBX_CROSSOVER_H

#endif //GA_SBX_CROSSOVER_H
