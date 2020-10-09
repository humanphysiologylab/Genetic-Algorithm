#ifndef GA_H
#define GA_H

#define GlobalSetupSize 7

FILE *f, *owle, *best, *avr, *test, *ap_best, *text, *sd, *ctrl_point;

double *AP_control, *AP_current, *Na_conc, *SD, *next_generation, *after_mut, *after_cross;
double best_organism_ap, scaling_factor, scaling_shift;
float *best_scaling_factor, *best_scaling_shift;

struct GlobalSetup{
    int number_organisms;
    int number_genes;
    int generations;
    int number_baselines;
    int elites;
    int autosave;
    int recording_frequency;
};


long i, j, h, baseline_counter, time_sum, c, time_counter, t_current;

#endif
