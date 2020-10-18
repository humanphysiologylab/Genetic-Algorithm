//
// Created by andrey on 1/22/20.
//

#ifndef GA_TOURNAMENT_SELECTION_H

#include <cstdio>
#include <cstdlib>

#include <ctime>
#include <vector>
#include <cassert>
#include "maleckar.h"

void tournament_selection(int *mpool, int mutant_number, std::vector<std::pair<double, int>> & sd_n_index, int number_of_ignored_losers);

#define GA_TOURNAMENT_SELECTION_H

#endif //GA_TOURNAMENT_SELECTION_H
