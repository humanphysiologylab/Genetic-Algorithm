#include "writing_to_output_files.h"

void writing_to_output_files(int NUMBER_ORGANISMS, int NUMBER_GENES, int NUMBER_BASELINES,
                             const double *genes, struct State *states_elite,
                             const double *AP_control, const double *AP_current,
                             float *best_scaling_factor, float *best_scaling_shift,
                             const std::vector<std::pair<double, int>> &sd_n_index,
                             int index_generation, int period_backup,
                             const int *CL, const int *TIME,
                             FILE *file_dump, FILE *file_ap_best) {

    // Parameters dump
    for (int i = 0; i < NUMBER_ORGANISMS; i++) {
        const double sd = sd_n_index[i].first;
        const int i_org = sd_n_index[i].second;
        fwrite(genes + i_org * NUMBER_GENES, sizeof(double), NUMBER_GENES, file_dump);
        fwrite(&sd, sizeof(double), 1, file_dump);
    }
    fflush(file_dump);


    // Backup
    if ((index_generation % period_backup == 0) && (index_generation != 0)) {
        FILE *file_backup;
        file_backup = fopen("./ga_output/backup.txt", "w");
        if (!file_backup) {
            printf("No autosave file!\n");
            exit(-1);
        }

        for (int i = 0; i < NUMBER_ORGANISMS; i++) {
            for (int j = 0; j < NUMBER_GENES; j++) {
                fprintf(file_backup, "%g\t", genes[i * NUMBER_GENES + j]);
            }
            fprintf(file_backup, "\n");
        }

        fclose(file_backup);
        printf("\nThe last autosave recording was made at %d generation\n", index_generation);
    }


    // AP of the best organism, all generations -> binary, last generation -> txt
    // Best organism is the first column, rescaled baseline is the second.
    std::ofstream ap_best_last;
    ap_best_last.open("./ga_output/ap_last.txt");

    int time_sum = 0;
    for (int i = 0; i < NUMBER_BASELINES; ++i) {
        time_sum += TIME[i];
    }
    const int index_elite = sd_n_index[0].second;
    const int t_start = index_elite * time_sum;

    int t = 0;

    for (int i_baseline = 0; i_baseline < NUMBER_BASELINES; i_baseline++) {
        for (int i_time = 0; i_time < TIME[i_baseline]; i_time++) {

            const double scaling_factor = best_scaling_factor[i_baseline + NUMBER_BASELINES * index_elite];
            const double scaling_shift = best_scaling_shift[i_baseline + NUMBER_BASELINES * index_elite];
            const double ap_best_organism = AP_current[t_start + t + i_time];
            const double ap_input_rescaled = AP_control[i_time + t] * scaling_factor + scaling_shift;

            fwrite(&ap_best_organism, sizeof(double), 1, file_ap_best);
            fwrite(&ap_input_rescaled, sizeof(double), 1, file_ap_best);

            ap_best_last << ap_best_organism << " " << AP_control[i_time + t] * scaling_factor + scaling_shift << "\n";

        }
        t += TIME[i_baseline];
    }
    fflush(file_ap_best);
    ap_best_last.close();

    // 8. Elite organisms states
    const char *path = "./states/State_elite_";
    const char *type = ".dat";
    char cl[512], full_path[512];

    for (int i = 0; i < NUMBER_BASELINES; i++) {
        sprintf(cl, "%d", CL[i]);
        snprintf(full_path, sizeof full_path, "%s%s%s", path, cl, type);

        FILE *file_state = fopen(full_path, "wb");

        int sizeof_state = sizeof(struct State);
        double a[sizeof_state];
        state2array(&states_elite[i], a);

        fwrite(a, sizeof(double), sizeof_state, file_state);
        fclose(file_state);
    }
}
