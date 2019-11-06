/* Fitness function calculation.
Date: 21.05.2019
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void scalar_multiplication(double* AP_control, double* AP_current, int length, double* XY, double* X1, double* Y1, double* XX, int* ones, int voltage_border)
{
    /*Scalar multiplication of two vectors*/
    int c, points_after;
    points_after = 0;
    *XY = 0.;
    *X1 = 0.;
    *Y1 = 0.;
    *XX = 0.;
    *ones = 0.;
    for (c = 0; c<length; c++)
    {
        if ((AP_current[c]>voltage_border)||(points_after == 1)){
            points_after = 1;
            *XY += AP_control[c] * AP_current[c];
            *X1 += AP_control[c] * 1.;
            *Y1 += AP_current[c] * 1.;
            *XX += AP_control[c] * AP_control[c];
            *ones += 1;
        }
    }
}


float calculate_half_height(double *array, int length) {

    float min = 1e10, max = -1e10;

    for (int i = 0; i < length; ++i) {
        if (array[i] < min) {
            min = array[i];
        }
        if (array[i] > max) {
            max = array[i];
        }
    }

    return (min + max) / 2;
}


int calculate_time_shift(double *array_reference, double *array_moving, int length) {

    float half_height_ref = calculate_half_height(array_reference, length);
    float half_height_mov = calculate_half_height(array_moving, length);

    int half_height_ref_index = -1;
    int half_height_mov_index = -1;

    for (int i = 0; i < length; ++i) {
        if ((half_height_ref_index == -1) && (array_reference[i] > half_height_ref)) {
            half_height_ref_index = i;
        }
        if ((half_height_mov_index == -1) && (array_moving[i] > half_height_mov)) {
            half_height_mov_index = i;
        }
        if ((half_height_mov_index != -1) && (half_height_ref_index != -1)) {
            break;
        }
    }

    return half_height_mov_index - half_height_ref_index;
}

float SD_calculation(double *AP_control, double *AP_current, float *best_scaling_factor, float *best_scaling_shift,
                     int *best_time_shift, int length) {

    double Variance;
    double AP_control_scaled, sd;
    double diff_between_potentials;
    int s, ss, points_after;
    double XY, X1, Y1, XX;
    int ones;
    double beta, alpha;

    int voltage_border = -20;          // start point (in mV) for SD calculation
    int minimal_amplitude = 80;         // mV
    int minimal_rest_potential = -70;   // mV

    int maximal_rest_potential = -90;   // mV
    sd = 0;
    ss = 0;

    int time_shift = calculate_time_shift(AP_current, AP_control, length);

    int length_truncated = length;
    double *AP_control_shifted = AP_control;
    double *AP_current_shifted = AP_current;

    if (time_shift < 0) {
        AP_current_shifted = AP_current - time_shift;
        length_truncated = length + time_shift;
    } else {
        AP_control_shifted = AP_control + time_shift;
        length_truncated = length - time_shift;
    }

    scalar_multiplication(AP_control_shifted, AP_current_shifted, length_truncated,
                          &XY, &X1, &Y1, &XX, &ones, voltage_border);

    if (XX != 0) {
        beta = (XY * X1 - Y1 * XX) / (X1 * X1 - ones * XX);
        alpha = (XY - beta * X1) / XX;

        if ((AP_control_shifted[length_truncated - 1] * alpha + beta) > minimal_rest_potential) {
            beta = minimal_rest_potential;
            alpha = (XY - beta * X1) / XX;
        }

        if ((AP_control_shifted[length_truncated - 1] * alpha + beta) < maximal_rest_potential) {
            beta = maximal_rest_potential;
            alpha = (XY - beta * X1) / XX;
        }

        if (alpha < minimal_amplitude){
            alpha = minimal_amplitude;
        }

        points_after = 0;
        for (s = 0; s < length_truncated; s++) {
            if ((AP_current_shifted[s] > voltage_border) || (points_after == 1)) {
                points_after = 1;
                AP_control_scaled = AP_control_shifted[s] * alpha + beta;
                diff_between_potentials = AP_control_scaled - AP_current_shifted[s];
                sd += diff_between_potentials * diff_between_potentials;
                ss += 1;
            }
        }

        Variance = sqrt(sd / (ss));

    } else {
        alpha = 0;
        beta = 0;
        Variance = 9e37;
    }

    *best_scaling_factor = alpha;
    *best_scaling_shift = beta;
    *best_time_shift = time_shift;

    return Variance;
}


void create_SD_index(int NUMBER_ORGANISMS, int *SD_index){
    int ll;
    for(ll=0;ll<NUMBER_ORGANISMS;ll++)
        SD_index[ll]=ll;
}


void insert_sort (int NUMBER_ORGANISMS, double *SD, int *SD_index){
    double v, v1;
    int index, organism;
    for(organism = 1; organism < NUMBER_ORGANISMS; organism++)
    {
        v = SD[organism];
        v1 = SD_index[organism];
        index = organism-1;
        while(index >= 0 && SD[index] > v)
        {
            SD[index+1] = SD[index];
            SD_index[index+1] = SD_index[index];
            index = index - 1;
        }
        SD[index+1] = v;
        SD_index[index+1] = v1;
    }
}


void fitness_function(double *AP_control, double *AP_current, float *best_scaling_factor, float *best_scaling_shift,
                      int *best_time_shift, int *TIME, double *SD, int *SD_index, int NUMBER_ORGANISMS,
                      int NUMBER_BASELINES, int time_sum) {
    int t_current = 0;
    int c, baseline_counter;
    for(c = 0;c < NUMBER_ORGANISMS; c++)
    {
        SD[c]=0;
        for(baseline_counter=0; baseline_counter<NUMBER_BASELINES;baseline_counter++)
        {
            SD[c] += SD_calculation(&AP_control[t_current], &AP_current[t_current + c * (time_sum)],
                                    &best_scaling_factor[baseline_counter + NUMBER_BASELINES * c],
                                    &best_scaling_shift[baseline_counter + NUMBER_BASELINES * c],
                                    &best_time_shift[baseline_counter + NUMBER_BASELINES * c],
                                    TIME[baseline_counter]); //Authomatic baseline scaling is implemented.
            
            t_current += TIME[baseline_counter];
        }
        t_current = 0;
    }
    
    create_SD_index(NUMBER_ORGANISMS, SD_index);
    insert_sort(NUMBER_ORGANISMS, SD, SD_index);

}


