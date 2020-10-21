//
// Created by andrey on 1/22/20.
//

#ifndef GA_INITIAL_POPULATION_H
#define GA_INITIAL_POPULATION_H

#include "maleckar.h"

void initial_population(double *next_generation, double *left_border, double *right_border, struct State *initial_state,
                        int NUMBER_ORGANISMS, int NUMBER_GENES, int NUMBER_BASELINES, int INIT_FROM_BACKUP_FILE);



#endif //GA_INITIAL_POPULATION_H
