#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ord_model/consts.h"

extern double ran2(long *);

void initial_population(double *next_generation, double *left_border, double *right_border, struct State *initial_state, int NUMBER_ORGANISMS, int NUMBER_GENES, int NUMBER_BASELINES, int READ_POPULATION){
    
    long seed, seed_negative;
    int j = 0;
    if (READ_POPULATION == 1){
        FILE *autosave_file;
        autosave_file = fopen("./ga_output/Autosave.txt", "r");
        if(!autosave_file) {printf("No autosave file!\n"); exit(-1);}
    
        for (j = 0; j < NUMBER_GENES*NUMBER_ORGANISMS; j++)
        {
            fscanf(autosave_file,"%lf\n",&next_generation[j]);
        }
        fclose(autosave_file);
    }
    else{
        int organism;
        seed = (long)time(NULL);
        seed_negative = -seed;
        ran2(&seed_negative);
        double rand_num = 0.;
        for (organism = 0; organism < NUMBER_ORGANISMS; organism++){
            /* Conductivities random generation */
            for (j = 0; j < NUMBER_GENES - 2 * NUMBER_BASELINES; j++)
            {
                rand_num = ran2(&seed);
                if (rand_num > 0.5){
                    next_generation[organism*NUMBER_GENES+j] = 1.0+(right_border[j]-1.0)*ran2(&seed);}
                else{
                    next_generation[organism*NUMBER_GENES+j] = 1.0+(left_border[j]-1.0)*ran2(&seed);}
            }
            /* Random generation of [Na+]i and [Ca++]nsr concentrations */
            for (j = 0; j < NUMBER_BASELINES; j++)
            {
                next_generation[(organism+1)*NUMBER_GENES-2*NUMBER_BASELINES+j] = initial_state[j].nai;
                next_generation[(organism+1)*NUMBER_GENES-NUMBER_BASELINES+j] = initial_state[j].cansr;
            }
        }
    }

}
