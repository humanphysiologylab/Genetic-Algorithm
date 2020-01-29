#ifndef GA_H
#define GA_H

#define GlobalSetupSize 7

FILE *owle, *best, *avr, *ap_best, *text, *sd, *ctrl_point;
time_t start1, end;

double*SD, *next_generation, *after_mut, *after_cross;
float *best_scaling_factor, *best_scaling_shift;
int *TIME;

struct GlobalSetup{
    int number_organisms;
    int number_genes;
    int generations;
    int number_baselines;
    int elites;
    int autosave;
    int recording_frequency;
};


long i, j, h, baseline_counter, time_sum, c, t_current;

#endif
