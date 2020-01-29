/* Fitness function calculation.
Date: 21.05.2019
*/

#include "fitness_function.h"

void scalar_multiplication(double *AP_control, double *AP_current, int length, double *XY, double *X1, double *Y1,
                           double *XX, int *ones, int voltage_border) {
    /*Scalar multiplication of two vectors*/
    int c, points_after;
    points_after = 0;
    *XY = 0.;
    *X1 = 0.;
    *Y1 = 0.;
    *XX = 0.;
    *ones = 0.;
    for (c = 0; c < length; c++) {
        if ((AP_current[c] > voltage_border) || (points_after == 1)) {
            points_after = 1;
            *XY += AP_control[c] * AP_current[c];
            *X1 += AP_control[c] * 1.;
            *Y1 += AP_current[c] * 1.;
            *XX += AP_control[c] * AP_control[c];
            *ones += 1;
        }
    }
}

float SD_calculation(double *AP_control, double *AP_current, float *best_scaling_factor, float *best_scaling_shift,
                     int length) {

    double sd = 0;

    for (int i = 0; i < length; i++) {
        double diff = AP_control[i] - AP_current[i];
        sd += diff * diff;
    }

    double rmse = sqrt(sd / length);

    *best_scaling_factor = 1;
    *best_scaling_shift = 0;

    return rmse;
}

void create_SD_index(int NUMBER_ORGANISMS, int *SD_index) {
    int ll;
    for (ll = 0; ll < NUMBER_ORGANISMS; ll++)
        SD_index[ll] = ll;
}

void insert_sort(int NUMBER_ORGANISMS, double *SD, int *SD_index) {
    double v, v1;
    int index, organism;
    for (organism = 1; organism < NUMBER_ORGANISMS; organism++) {
        v = SD[organism];
        v1 = SD_index[organism];
        index = organism - 1;
        while (index >= 0 && SD[index] > v) {
            SD[index + 1] = SD[index];
            SD_index[index + 1] = SD_index[index];
            index = index - 1;
        }
        SD[index + 1] = v;
        SD_index[index + 1] = v1;
    }
}


void fitness_function(double *AP_control, double *AP_current, float *best_scaling_factor, float *best_scaling_shift,
                      int *TIME, double *SD, int *SD_index, int NUMBER_ORGANISMS, int NUMBER_BASELINES, int time_sum) {
    int t_current = 0;
    int c, baseline_counter;
    for (c = 0; c < NUMBER_ORGANISMS; c++) {
        SD[c] = 0;
        for (baseline_counter = 0; baseline_counter < NUMBER_BASELINES; baseline_counter++) {
            SD[c] += SD_calculation(&AP_control[t_current], &AP_current[t_current + c * (time_sum)],
                                    &best_scaling_factor[baseline_counter + NUMBER_BASELINES * c],
                                    &best_scaling_shift[baseline_counter + NUMBER_BASELINES * c],
                                    TIME[baseline_counter]); //Authomatic baseline scaling is implemented.

            if (isnan(SD[c])) {
                printf("SD is NAN!; Set to 1e100\n");
                SD[c] = 1e100;
            }

            t_current += TIME[baseline_counter];
        }
        t_current = 0;
    }

    create_SD_index(NUMBER_ORGANISMS, SD_index);
    insert_sort(NUMBER_ORGANISMS, SD, SD_index);

}


