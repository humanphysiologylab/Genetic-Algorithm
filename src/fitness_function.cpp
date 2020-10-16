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
    double Variance;
    double AP_control_scaled, sd = 0;
    int s, points_after;
    double XY, X1, Y1, XX;
    int ones;
    double beta, alpha;

    int voltage_border = -20;          // start point (in mV) for SD calculation
    int minimal_amplitude = 80;         // mV
    int minimal_rest_potential = -60;   // mV

    int maximal_rest_potential = -90;   // mV


//    scalar_multiplication(AP_control, AP_current, length, &XY, &X1, &Y1, &XX, &ones, voltage_border);

/*    if(XX!=0)
    {
    	beta = (XY * X1 - Y1 * XX)/(X1 * X1 - ones * XX);
   	alpha = (XY - beta * X1)/XX;
//	if (alpha < minimal_amplitude){
//            alpha = minimal_amplitude;
//            beta = (Y1 - alpha * X1)/ones;
//    	}


    	if ((AP_control[0] * alpha + beta) > minimal_rest_potential){
 		beta = minimal_rest_potential;
        	alpha = (XY - beta * X1)/XX;
    	}
	if((AP_control[0] * alpha + beta) < maximal_rest_potential)
    	{
        	beta = maximal_rest_potential;
        	alpha = (XY - beta * X1)/XX;
    	}
	if (alpha < minimal_amplitude) alpha = minimal_amplitude;*/

    points_after = 0;
    for (s = 0; s < length; s++) {
//       		if ((AP_current[s] > voltage_border)||(points_after == 1)){
//        		points_after = 1;
//          		AP_control_scaled = AP_control[s] * alpha + beta;
        const double diff_between_potentials = AP_control[s] - AP_current[s];
        if (AP_current[s] == 0.0) { //?? this never equals True probably
            printf("AP_current[%d] is off\n", s);
            getc(stdin);
        }
        sd += diff_between_potentials * diff_between_potentials;
        
        //       	}
    }

    Variance = sqrt(sd / length);
/*    }
    else
    {
	alpha=0;
	beta=0;
	Variance=9e37;

    }*/



    *best_scaling_factor = 1;//alpha;
    *best_scaling_shift = 0; //beta;

    return Variance;
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


