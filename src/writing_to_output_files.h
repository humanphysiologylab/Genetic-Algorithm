//
// Created by andrey on 1/22/20.
//

#ifndef GA_WRITING_TO_OUTPUT_FILES_H

#include <stdio.h>
#include <stdlib.h>

#include "maleckar.h"

void writing_to_output_files(FILE *best, FILE *avr, FILE *owle, FILE *ctrl_point, FILE *text, FILE *sd, FILE *ap_best,
                             double *SD, int *SD_index, double average, double *best_parameters, int NUMBER_GENES,
                             int NUMBER_ORGANISMS, double *next_generation, int cntr, int recording_frequency,
                             int NUMBER_BASELINES, struct State *elite_state, int *CL, double *AP_current,
                             double *AP_control, float *best_scaling_factor, float *best_scaling_shift, int elite_index,
                             int *TIME, int time_index);

#define GA_WRITING_TO_OUTPUT_FILES_H

#endif //GA_WRITING_TO_OUTPUT_FILES_H
