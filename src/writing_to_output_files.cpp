#include "writing_to_output_files.h"

void writing_to_output_files(FILE *best, FILE *avr, FILE *owle, FILE *ctrl_point, FILE *text, FILE *sd, FILE *ap_best,
                             double *SD, int *SD_index, double average, double *best_parameters, int NUMBER_GENES,
                             int NUMBER_ORGANISMS, double *next_generation, int cntr, int recording_frequency,
                             int NUMBER_BASELINES, struct State *elite_state, int *CL, double *AP_current,
                             double *AP_control, float *best_scaling_factor, float *best_scaling_shift, int elite_index,
                             int *TIME, int time_index) {

    int i, j, t, baseline_counter;
    double best_organism_ap;

    // 1. SD of best organism in population
    fprintf(best, "%f\n", SD[0]);
    fflush(best);
    // 2. Average SD in generation
    fprintf(avr, "%f\t\n", average);

    // 3. Parameters of the best organism
    for (i = 0; i < NUMBER_GENES; i++) {
        fprintf(owle, "%g\t", best_parameters[i]);
    }
    fprintf(owle, "\n");

    // 4. Parameters of each organism
    for (i = 0; i < NUMBER_ORGANISMS; i++) {
        for (j = 0; j < NUMBER_GENES; j++) {
            fprintf(text, "%g\t", next_generation[i * NUMBER_GENES + j]);
        }
        fprintf(text, "%d\n", cntr);
    }

    // 5. SD values of each organism in population
    for (i = 0; i < NUMBER_ORGANISMS; i++) fprintf(sd, "%lf\t %d\t %d\n", SD[i], SD_index[i], cntr);

    // 6. Autosave
    if ((cntr % recording_frequency == 0) && (cntr != 0)) {
        ctrl_point = fopen("./ga_output/Autosave.txt", "w");
        if (!ctrl_point) {
            printf("No autosave file!\n");
            exit(-1);
        }

        for (i = 0; i < NUMBER_ORGANISMS; i++) {
            for (j = 0; j < NUMBER_GENES; j++) {
                fprintf(ctrl_point, "%g\t", next_generation[i * NUMBER_GENES + j]);
            }
            fprintf(ctrl_point, "\n");
        }
        fclose(ctrl_point);
        printf("\nThe last autosave recording was made at %d generation\n", cntr);
    }

    // 7. AP of the best organism
    // Best organism is in first column, rescaled baseline is in the second.
    t = 0;
    double scaling_factor, scaling_shift;
    for (baseline_counter = 0; baseline_counter < NUMBER_BASELINES; baseline_counter++) {
        for (i = 0; i < TIME[baseline_counter]; i++) {
            scaling_factor = best_scaling_factor[baseline_counter + NUMBER_BASELINES * elite_index];
            scaling_shift = best_scaling_shift[baseline_counter + NUMBER_BASELINES * elite_index];
            best_organism_ap = AP_current[time_index + t + i];
            fprintf(ap_best, "%f\t%f\n", best_organism_ap, AP_control[i + t] * scaling_factor + scaling_shift);
        }
        t += TIME[baseline_counter];
    }
    fflush(ap_best);
    // 8. Elite organisms states
    const char *path = "./states/State_elite_";
    const char *type = ".dat";
    char cl[512], full_path[512];

    for (i = 0; i < NUMBER_BASELINES; i++) {
        sprintf(cl, "%d", CL[i]);
        snprintf(full_path, sizeof full_path, "%s%s%s", path, cl, type);

        FILE *all_state = fopen(full_path, "wb");
        fwrite(&elite_state[i], sizeof(struct State), 1, all_state);
        fclose(all_state);
    }


}
