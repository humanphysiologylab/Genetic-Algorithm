/* Fitness function calculation.
Date: 21.05.2019
*/
#include <math.h>
#include <cstdio> 
#include "fitness_function.h"


void scalar_multiplication(double *AP_control, double *AP_current, int length, double *XY, double *X1, double *Y1,
                           double *XX, int *ones, int voltage_border)
{
    /*Scalar multiplication of two vectors*/
    *XY = 0.;
    *X1 = 0.;
    *Y1 = 0.;
    *XX = 0.;
    *ones = 0.;
    for (int c = 0, points_after = 0; c < length; c++) {
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
                     int length, const double *weight)
{
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
        sd += weight[s] * diff_between_potentials * diff_between_potentials;
        
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



void fitness_function(double *AP_control, double *AP_current, float *best_scaling_factor, float *best_scaling_shift,
                      int *TIME, std::vector<std::pair<double, int>> & sd_n_index, int NUMBER_ORGANISMS, int NUMBER_BASELINES, int time_sum, const double *weight)
{
    for (int c = 0; c < NUMBER_ORGANISMS; c++) {
        double res = 0;
        for (long t_current = 0, baseline_counter = 0; baseline_counter < NUMBER_BASELINES; baseline_counter++) {
            res += SD_calculation(&AP_control[t_current], &AP_current[t_current + c * (time_sum)],
                                    &best_scaling_factor[baseline_counter + NUMBER_BASELINES * c],
                                    &best_scaling_shift[baseline_counter + NUMBER_BASELINES * c],
                                    TIME[baseline_counter], &weight[t_current]); //Authomatic baseline scaling is implemented.
            t_current += TIME[baseline_counter];
        }
        if (isnan(res)) {
            printf("SD is NAN!; Set to 1e100\n");
            res = 1e100;
        }
        //save error to sd_n_index
        sd_n_index[c].first = res;
        sd_n_index[c].second = c;
    }
}
