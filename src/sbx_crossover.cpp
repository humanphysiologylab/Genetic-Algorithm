#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define CROSSRATE 0.9 //probability of crossover

extern double ran2(long *);

void sbx_crossover(double *next_generation, double *after_cross, int *mpool, double *left_border, double *right_border, int NUMBER_ORGANISMS, int NUMBER_GENES){
    int etaC = 10; //The order of the polynomial for the SBX crossover
    int i, ii, jj;
    int num_1, num_2;
    double sw_prob;
    double cr_prob;
    double p1,p2,c1,c2,random;
    double beta, alpha, beta_q;
    long seed, seed_negative;
    
    seed = (long)time(NULL);
    seed_negative = -seed;
    ran2(&seed_negative);

    for (i=0;i<NUMBER_ORGANISMS/2;i++)
    {
        num_1 = mpool[2*i];
        num_2 = mpool[2*i+1];
        cr_prob = ran2(&seed);
        
        if (cr_prob <= CROSSRATE)
        {
            for (ii = 0; ii < NUMBER_GENES; ii++)
            {
                sw_prob = ran2(&seed);
                
                if (next_generation[num_1*NUMBER_GENES+ii] > right_border[ii] || next_generation[num_1*NUMBER_GENES+ii] < left_border[ii]){
                    next_generation[num_1*NUMBER_GENES+ii] = left_border[ii] + (right_border[ii]-left_border[ii])*ran2(&seed);
                }
                if (next_generation[num_2*NUMBER_GENES+ii] > right_border[ii] || next_generation[num_2*NUMBER_GENES+ii] < left_border[ii]){
                    next_generation[num_2*NUMBER_GENES+ii] = left_border[ii] + (right_border[ii]-left_border[ii])*ran2(&seed);
                }
                
                if (sw_prob < 0.5) // "NO"
                {
                    after_cross[2*i*NUMBER_GENES+ii] = next_generation[num_1*NUMBER_GENES+ii];
                    after_cross[(2*i+1)*NUMBER_GENES+ii] = next_generation[num_2*NUMBER_GENES+ii];
                }
                
                else //"YES"
                {
                    if (next_generation[num_1*NUMBER_GENES+ii] < next_generation[num_2*NUMBER_GENES+ii]) {
                        p1 = next_generation[num_1*NUMBER_GENES+ii];
                        p2 = next_generation[num_2*NUMBER_GENES+ii];
                    }
                    else {
                        p1 = next_generation[num_2*NUMBER_GENES+ii];
                        p2 = next_generation[num_1*NUMBER_GENES+ii];
                    }
                    // Calculate beta
                    if (right_border[ii] - p2 < p1-left_border[ii]){
                        beta = (p2-p1)/(2.0 * right_border[ii]-p2-p1);//
                    }
                    else if((p1 == 0.0) && (p2 == 0.0)){
                        after_cross[2*i*NUMBER_GENES+ii] = p1;
                        after_cross[(2*i+1)*NUMBER_GENES+ii] = p2;
                    }
                    else  {
                        beta  = (p2-p1)/(p2+p1-2.0*left_border[ii]);
                    }
                    random = ran2(&seed);
                    
                    alpha = (2.0 - pow(beta,(double)(etaC+1)))*random;
                    
                    if (alpha<=1.0){
                        beta_q = pow(alpha,1.0/((double)(etaC+1)));
                    }
                    else {
                        beta_q = pow(1.0/(2.0-alpha),1.0/((double)(etaC+1)));//
                    }
                    c1 = 0.5*((p1+p2) + beta_q*(p1-p2));
                    c2 = 0.5*((p1+p2) + beta_q*(p2-p1));
                    
                    after_cross[2*i*NUMBER_GENES+ii] = c1;
                    after_cross[(2*i+1)*NUMBER_GENES+ii] = c2;
                }
                
            }
        }
        else
        {
            for(jj = 0; jj < NUMBER_GENES; jj++)
            {
                after_cross[2*i*NUMBER_GENES+jj] = next_generation[num_1*NUMBER_GENES+jj];
                after_cross[(2*i+1)*NUMBER_GENES+jj] = next_generation[num_2*NUMBER_GENES+jj];
            }
        }
    }
}
