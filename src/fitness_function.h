//
// Created by andrey on 1/22/20.
//

#ifndef GA_FITNESS_FUNCTION_H
#define GA_FITNESS_FUNCTION_H

#include <utility>
#include <vector>

void scalar_multiplication(double* AP_control, double* AP_current, int length, double* XY, double* X1, double* Y1, double* XX, int* ones, int voltage_border);

float SD_calculation(double* AP_control, double* AP_current, float* best_scaling_factor, float* best_scaling_shift, int length);



void fitness_function(double *AP_control, double *AP_current, float *best_scaling_factor, float *best_scaling_shift,
                      int *TIME, std::vector<std::pair<double, int>> & sd_n_index, int NUMBER_ORGANISMS, int NUMBER_BASELINES, int time_sum);




#endif //GA_FITNESS_FUNCTION_H
