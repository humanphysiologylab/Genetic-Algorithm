#include "tournament_selection.h"

void tournament_selection(int *mpool, double *SD, struct State *state_struct, struct State *state_struct_rewrite,
                          int NUMBER_ORGANISMS, int NUMBER_BASELINES) {

    int i, j, tt, index_1, index_2;
    int c = 0;
    double copy_1[NUMBER_ORGANISMS];
    double copy_2[NUMBER_ORGANISMS];

    for (tt = 0; tt < NUMBER_ORGANISMS; tt++) {
        copy_1[tt] = SD[tt];
        copy_2[tt] = SD[tt];
    }

    srand(time(NULL));
    while (c < NUMBER_ORGANISMS / 2) {
        //first half of the mating pool
        index_1 = rand() % NUMBER_ORGANISMS;
        index_2 = rand() % NUMBER_ORGANISMS;
        if ((index_1 != index_2) && (copy_1[index_1] != 0) && (copy_1[index_2] != 0)) {
            if (copy_1[index_1] <= copy_1[index_2]) {
                mpool[2 * c] = index_1;
                for (i = 0; i < NUMBER_BASELINES; i++)
                    state_struct_rewrite[i + 2 * c * NUMBER_BASELINES] = state_struct[i + index_1 * NUMBER_BASELINES];
                //printf("Winner in the first group: %d\t\n", mpool[2*c]);
            } else {
                mpool[2 * c] = index_2;
                for (i = 0; i < NUMBER_BASELINES; i++)
                    state_struct_rewrite[i + 2 * c * NUMBER_BASELINES] = state_struct[i + index_2 * NUMBER_BASELINES];
                //printf("Winner in the first group: %d\t\n", mpool[2*c]);
            }
            copy_1[index_1] = 0;
            copy_1[index_2] = 0;
            c += 1;
        }
    }
    c = 0;
    while (c < NUMBER_ORGANISMS / 2) {
        //second half of the mating pool
        index_1 = rand() % NUMBER_ORGANISMS;
        index_2 = rand() % NUMBER_ORGANISMS;
        if ((index_1 != index_2) && (copy_2[index_1] != 0) && (copy_2[index_2] != 0)) {
            if (copy_2[index_1] <= copy_2[index_2]) {
                mpool[2 * c + 1] = index_1;
                for (i = 0; i < NUMBER_BASELINES; i++)
                    state_struct_rewrite[i + (2 * c + 1) * NUMBER_BASELINES] = state_struct[i + index_1 *
                                                                                                NUMBER_BASELINES];
            } else {
                mpool[2 * c + 1] = index_2;
                for (i = 0; i < NUMBER_BASELINES; i++)
                    state_struct_rewrite[i + (2 * c + 1) * NUMBER_BASELINES] = state_struct[i + index_1 *
                                                                                                NUMBER_BASELINES];
            }
            copy_2[index_1] = 0;
            copy_2[index_2] = 0;
            c += 1;
        }
    }

    for (i = 0; i < NUMBER_ORGANISMS; i++)
        for (j = 0; j < NUMBER_BASELINES; j++)
            state_struct[j + i * NUMBER_BASELINES] = state_struct_rewrite[j + i * NUMBER_BASELINES];

}
