//
// Created by andrey on 1/22/20.
//

#ifndef GA_TOURNAMENT_SELECTION_H

#include <cstdio>
#include <cstdlib>

#include <ctime>
#include "maleckar.h"

void tournament_selection(int *mpool, double *SD, struct State *state_struct, struct State *state_struct_rewrite, int NUMBER_ORGANISMS, int NUMBER_BASELINES);

#define GA_TOURNAMENT_SELECTION_H

#endif //GA_TOURNAMENT_SELECTION_H
