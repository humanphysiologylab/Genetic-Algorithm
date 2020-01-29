//
// Created by andrey on 1/22/20.
//

#ifndef GA_FITNESS_FUNCTION_H

#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <ctime>

void scalar_multiplication(double* AP_control, double* AP_current, int length, double* XY, double* X1, double* Y1, double* XX, int* ones, int voltage_border);

float SD_calculation(double* AP_control, double* AP_current, float* best_scaling_factor, float* best_scaling_shift, int length);

void create_SD_index(int NUMBER_ORGANISMS, int *SD_index);

void insert_sort (int NUMBER_ORGANISMS, double *SD, int *SD_index);

void fitness_function(double *AP_control, double *AP_current, float *best_scaling_factor, float *best_scaling_shift,
                      int *TIME, double *SD, int *SD_index, int NUMBER_ORGANISMS, int NUMBER_BASELINES, int time_sum);


#define GA_FITNESS_FUNCTION_H

#endif //GA_FITNESS_FUNCTION_H
