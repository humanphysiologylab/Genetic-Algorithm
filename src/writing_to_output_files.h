//
// Created by andrey on 1/22/20.
//

#ifndef GA_WRITING_TO_OUTPUT_FILES_H
#define GA_WRITING_TO_OUTPUT_FILES_H
#include <vector>
#include <utility>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "maleckar.h"

void writing_to_output_files(int NUMBER_ORGANISMS, int NUMBER_GENES, int NUMBER_BASELINES,
                             const double *genes, struct State *states_elite,
                             const double *AP_control, const double *AP_current,
                             float *best_scaling_factor, float *best_scaling_shift,
                             const std::vector<std::pair<double, int>> & sd_n_index,
                             int index_generation, int period_backup,
                             const int *CL, const int *TIME,
                             FILE *file_dump, FILE *file_ap_best);



#endif //GA_WRITING_TO_OUTPUT_FILES_H
