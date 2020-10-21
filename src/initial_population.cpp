#include <cstdio>
#include <iostream>
#include <cmath>
#include "initial_population.h"
#include "random_number_generator.h"

void initial_population(double *next_generation, double *left_border, double *right_border, struct State *initial_state,
                        int NUMBER_ORGANISMS, int NUMBER_GENES, int NUMBER_BASELINES, int INIT_FROM_BACKUP_FILE) {

    long seed, seed_negative;
    int j = 0;
    if (INIT_FROM_BACKUP_FILE == 1) {
        FILE *autosave_file;
        autosave_file = fopen("./ga_output/backup.txt", "r");
        if (!autosave_file) {
            printf("No autosave file!\n");
            exit(-1);
        }

        for (j = 0; j < NUMBER_GENES * NUMBER_ORGANISMS; j++) {
            fscanf(autosave_file, "%lf\n", &next_generation[j]);
        }
        fclose(autosave_file);
    } else {
        int organism;
        seed = (long) time(NULL);
        seed_negative = -seed;
        ran2(&seed_negative);
        double rand_num = 0.;//?????
        for (organism = 0; organism < NUMBER_ORGANISMS; organism++) {

            /* Conductivities random generation */
            for (j = 0; j < NUMBER_GENES - 3 * NUMBER_BASELINES; j++) {

                if (left_border[j] >= right_border[j]) {
                    std::cout << "Error in input file!\n";
                    std::cout << "Range for parameter № " << j << " is invalid: [" << left_border[j] << ", " << right_border[j] << "]\n";
                    exit(-1);
                }

                if ((left_border[j] < 1) && (1 < right_border[j])) {
                    if (ran2(&seed) > 0.5) {
                        next_generation[organism * NUMBER_GENES + j] = pow(right_border[j], ran2(&seed));
                    } else {
                        next_generation[organism * NUMBER_GENES + j] = pow(left_border[j], 1 - ran2(&seed));
                    }
                } else {
                    double r = ran2(&seed);
                    next_generation[organism * NUMBER_GENES + j] = pow(right_border[j], r) * pow(left_border[j], 1 - r);
                }
            }

            /* Concentrations */
            for (j = 0; j < NUMBER_BASELINES; j++) {
                next_generation[(organism + 1) * NUMBER_GENES - 3 * NUMBER_BASELINES + j] = initial_state[j].Na_i;
                next_generation[(organism + 1) * NUMBER_GENES - 2 * NUMBER_BASELINES + j] = initial_state[j].Ca_rel;
                next_generation[(organism + 1) * NUMBER_GENES - 1 * NUMBER_BASELINES + j] = initial_state[j].K_i;
            }
        }
    }

}
