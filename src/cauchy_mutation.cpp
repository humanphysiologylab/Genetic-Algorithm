#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592653589793238462643
#define MUTRATE 0.9 //probability of mutation

extern double ran2(long *);


void box_muller_transform(double *uniform_vector, int NUMBER_GENES, long seed){
    double u, v, s, z1, sq, rand_ud;
    double mu, epsilon;
    double vector_len = 0.;
    int item, zi;
    
    mu = 0.;
    epsilon = 1.;
    
    for (zi = 0; zi < NUMBER_GENES; zi++){
        do{
            rand_ud = ran2(&seed);
            if (rand_ud > 0.5){
                u = ran2(&seed);}
            else{
                u = -ran2(&seed);}
            
            rand_ud = ran2(&seed);
            if (rand_ud > 0.5){
                v = ran2(&seed);}
            else{
                v = -ran2(&seed);}
            
            s = sqrt(u*u+v*v);
        } while  ((s == 0) || (s>=1));
        
        sq = sqrt(-2 * log(s)/s);
        
        z1 = u * sq;
        
        uniform_vector[zi] = mu + z1*epsilon;
        
    }
    
    for (item = 0; item < NUMBER_GENES; item++){
        vector_len += uniform_vector[item]*uniform_vector[item];
    }
    vector_len = sqrt(vector_len);
    
    for (item = 0; item < NUMBER_GENES; item++){
        uniform_vector[item] = uniform_vector[item]/vector_len;
    }
}


void cauchy_mutation(double *after_mut, double *after_cross, double *left_border, double *right_border, int NUMBER_ORGANISMS, int NUMBER_GENES){
    double mut_prob;
    double rnd;
    double delta1, delta2, delta_q;
    double sigma = 0.18;
    int mu = 0;
    
    int i = 0;
    int j = 0;
    double r_min = 0;
    double r_max = 0;
    long seed, seed_negative;
    
    double shift, sq_left, sq_right, x_max, x_min;
    double uniform_vector[NUMBER_GENES];
    double right_end[NUMBER_GENES], left_end[NUMBER_GENES];
    
    seed = (long)time(NULL);
    seed_negative = -seed;
    ran2(&seed_negative);
    
    for (i=0; i < NUMBER_ORGANISMS;i++)
    {
        mut_prob = ran2(&seed);
        
        if (mut_prob <= 1-MUTRATE){
            for (j=0; j<NUMBER_GENES; j++){
                after_mut[i*NUMBER_GENES+j] = after_cross[i*NUMBER_GENES+j];
            }
        }
        else {
            int iter = 0;
            double distance_min1 = 1000000.;
            double distance_min2 = 1000000.;
            double distance = 0;
            
            box_muller_transform(uniform_vector, NUMBER_GENES, seed);
            
            /* right border */
            for(iter = 0;iter<NUMBER_GENES;iter++){
                if (uniform_vector[iter]>0){
                    distance = (right_border[iter]-after_cross[i*NUMBER_GENES+iter])/uniform_vector[iter];
                }
                else{
                    distance = (left_border[iter]-after_cross[i*NUMBER_GENES+iter])/uniform_vector[iter];
                }
                if (distance<distance_min1){
                    distance_min1 = distance;
                }
            }
            /* left border*/
            for(iter = 0;iter<NUMBER_GENES;iter++){
                if (-uniform_vector[iter]>0){
                    distance = (right_border[iter]-after_cross[i*NUMBER_GENES+iter])/(-uniform_vector[iter]);
                }
                else{
                    distance = (left_border[iter]-after_cross[i*NUMBER_GENES+iter])/(-uniform_vector[iter]);
                }
                if (distance<distance_min2){
                    distance_min2 = distance;
                }
            }
            /* calculate 2 edge points: */
            for(iter = 0;iter<NUMBER_GENES;iter++){
                right_end[iter] = after_cross[i*NUMBER_GENES+iter] + distance_min1 * uniform_vector[iter];
                left_end[iter] = after_cross[i*NUMBER_GENES+iter] + distance_min2 * (-uniform_vector[iter]);
            }
            
            /* length of 2 vectors:*/
            sq_left = 0;
            sq_right = 0;
            for(iter = 0;iter<NUMBER_GENES;iter++){
                sq_right += (right_end[iter]-after_cross[i*NUMBER_GENES+iter])*(right_end[iter]-after_cross[i*NUMBER_GENES+iter]);
                sq_left += (left_end[iter]-after_cross[i*NUMBER_GENES+iter])*(left_end[iter]-after_cross[i*NUMBER_GENES+iter]);
            }
            x_max = sqrt(sq_right);
            x_min = - sqrt(sq_left);
            
            delta1 = -x_min;
            delta2 = x_max;
            
            r_min = (atan(-delta1/sigma))/pi + 0.5;
            r_max = (atan(delta2/sigma))/pi + 0.5;
            
            rnd = r_min+(r_max-r_min)*ran2(&seed);
            delta_q = mu + sigma*tan(pi*(rnd-1./2));
            
            shift = delta_q;
            
            for(iter = 0;iter<NUMBER_GENES;iter++){
                after_mut[i*NUMBER_GENES+iter] = after_cross[i*NUMBER_GENES+iter] + shift*uniform_vector[iter];
            }
        }
        
    }
}



