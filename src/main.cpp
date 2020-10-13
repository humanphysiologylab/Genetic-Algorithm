// Main GA module
//
// Genetic Algorithm implementation
//
// Authors:
//     Dmitrii Smirnov <dmitrii.smirnov@phystech.edu>
//     Roman Syunyaev <roman.syunyaev@gmail.com>
//     Alexander Timofeev <richardstallman42@gmail.com>
//     Laboratory of Human Physiology, Moscow Institute of Physics and Technology, 2020
//
// License:
//     Redistribution and use of the code with or without modification
//     are permitted for any scientific, research and educational purposes.
//
//     Redistribution and use of the code with or without modification
//     are prohibited for any commercial purposes.


#include <mpi.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <iostream>
#include <cassert>
#include <filesystem>

#include "cauchy_mutation.h"
#include "fitness_function.h"
#include "initial_population.h"
//#include "ord_model/consts.h"
#include "sbx_crossover.h"
#include "tournament_selection.h"
#include "writing_to_output_files.h"
//#include "ord_model/atrium.h"
#include "maleckar.h"

void scanf_baseline(int j0, int j1, FILE *ff, double *AP_control) {
	//?? what is a baseline
	assert(j0 >= 0);
    for (int j = j0; j < j1; j++) {
        fscanf(ff, "%lf\n", &AP_control[j]);
    }
}

void normalize_baseline(int j0, int j1, double *AP_control) {
    // max (min) is the time of baseline maximum (minimum)
    assert(j0 >= 0);
    if (j0 >= j1)
		return;

    double max = AP_control[j0], min = AP_control[j0];
    for (int i = j0; i < j1; i++) {
        if (min > AP_control[i]) min = AP_control[i];
        if (max < AP_control[i]) max = AP_control[i];
    }
    for (int i = j0; i < j1; i++)
		AP_control[i] = (AP_control[i] - min) / (max - min);
}

int time_array(long &time_sum, FILE *f) {
	//?? dangerous function, we add counter to time_sum, but return it at the same time... why?
    int counter = 0;
	int c;
	do {
		c = fgetc(f);
		if (c == '\n')
			counter++;
	} while (c != EOF);

    time_sum += counter;
    return counter;
}

static char *read_line(char *pcBuf, int iMaxSize, FILE *fStream) {
    //tbh weird function.
    //it returns pcBuf containing next non-zero line from fStream
    //which is not a comment starting with #.
    //seems like p always points at pcBuf.
    //if eof then pcBuf stores \0.
    char *p;

    do {
        p = NULL;
        *pcBuf = '\0';

        if (fgets(pcBuf, iMaxSize, fStream) == NULL) {
            printf("Error! Unable to read line!\n");
            exit(-1);
        }

        if (feof(fStream)) break;
        p = strtok(pcBuf, "\n\r");
    } while ((p == NULL) || (*p == '#'));

    return p;
}





struct GlobalSetup{
    int number_organisms;
    int number_genes;
    int generations;
    int number_baselines;
    int elites;
    int autosave;
    int recording_frequency;
};
const int GlobalSetupItemsNumber = 7;

int main(int argc, char *argv[]) {
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* create an MPI type for struct GlobalSetup */
    int          blocklengths[GlobalSetupItemsNumber] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[GlobalSetupItemsNumber] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Datatype GlobalSetupMPI;
    MPI_Aint     displacements[GlobalSetupItemsNumber];

    displacements[0] = offsetof(GlobalSetup, number_organisms);
    displacements[1] = offsetof(GlobalSetup, number_genes);
    displacements[2] = offsetof(GlobalSetup, generations);
    displacements[3] = offsetof(GlobalSetup, number_baselines);
    displacements[4] = offsetof(GlobalSetup, elites);
    displacements[5] = offsetof(GlobalSetup, autosave);
    displacements[6] = offsetof(GlobalSetup, recording_frequency);

    MPI_Type_create_struct(GlobalSetupItemsNumber, blocklengths, displacements, types, &GlobalSetupMPI);
    MPI_Type_commit(&GlobalSetupMPI);
    /* done */

    /* Create MPI structure for initial state */
    MPI_Datatype StateVectorMPI;
    MPI_Type_contiguous(STATE_ARRAY_SIZE, MPI_DOUBLE, &StateVectorMPI);
    MPI_Type_commit(&StateVectorMPI);
    /* done */





    long t_current;

    long time_sum;
    long baseline_counter;


double *AP_control, *AP_current, *Na_conc, *SD, *next_generation, *after_mut, *after_cross;
double best_organism_ap, scaling_factor, scaling_shift;
float *best_scaling_factor, *best_scaling_shift;



    //open files and make sure it was successful
    int init_status = 0;

    FILE *f, *owle, *best, *avr, *test, *ap_best, *text, *sd, *ctrl_point;

    if (rank == 0) {
		std::filesystem::create_directory("ga_output");
		owle = fopen("./ga_output/low_err.txt", "w");
		best = fopen("./ga_output/best.txt", "w");
		avr = fopen("./ga_output/average.txt", "w");
		ap_best = fopen("./ga_output/AP_best.txt", "w");
		text = fopen("./ga_output/output1.txt", "w");
		sd = fopen("./ga_output/SD.txt", "w");

		if (!(owle && best && avr && ap_best && text && sd)) {
			printf("Cannot open files for output\n");
			init_status = -1;
		}

		if (argc != 2) {
			printf("Error! GA input file required! \n");
			init_status = -1;
		}
	}

    MPI_Bcast(&init_status, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (init_status != 0) {
		MPI_Finalize();
		return 0;
	}
	//initialization complete


    //timer
    double start1, end;
    start1 = MPI_Wtime();






    //??
    MPI_Request reqs[6];
    MPI_Status stats[12];



    double *IA = 0; //Stimulation current, A/F
    double *right_border = 0, *left_border = 0; //min and max possible values of parameters
    int *CL = 0; //Cycle lengths
    
    //1 - include ISO model
    //0 - without ISO model
    int *ISO = 0;
    
    char baseline_file_names[256][256], statedat_file_names[256][256];

    GlobalSetup gs;

    if (rank == 0) {
        const int ciBufSize = 1024;
        char caBuf[ciBufSize];
        FILE *fInput, *fOutput;
        char *InputFile = argv[1];

        fInput = fopen(InputFile, "r");
        if (!fInput) {
            printf("Cannot open GA input file!\n");
            return -1;
        }

        gs.number_organisms = atoi(read_line(caBuf, ciBufSize, fInput));
        gs.number_genes = atoi(read_line(caBuf, ciBufSize, fInput));
        gs.generations = atoi(read_line(caBuf, ciBufSize, fInput));

        //Parameters ranges
        if ((left_border = (double *) malloc(sizeof(double) * gs.number_genes)) == NULL) {
            puts("The 'left_border' array isn't created!");
            exit(-1);
        }
        if ((right_border = (double *) malloc(sizeof(double) * gs.number_genes)) == NULL) {
            puts("The 'right_border' array isn't created!");
            exit(-1);
        }

        for (int i = 0; i < gs.number_genes; i++) {
            char *token, *ptoken;
            token = read_line(caBuf, ciBufSize, fInput);

            if ((ptoken = strtok(token, " \t\n\r")) == NULL) {
                fclose(fInput);
                printf("Error in the input file!\n");
                exit(-1);
            }
            left_border[i] = atof(ptoken);

            if ((ptoken = strtok(NULL, " \t\n\r")) == NULL) {
                fclose(fInput);
                printf("Error in the input file!\n");
                exit(-1);
            }
            right_border[i] = atof(ptoken);
        }

        gs.number_baselines = atoi(read_line(caBuf, ciBufSize, fInput));

        if ((CL = (int *) malloc(sizeof(int) * gs.number_baselines)) == NULL) {
            puts("The 'CL' array isn't created!");
            exit(-1);
        }
        if ((IA = (double *) malloc(sizeof(double) * gs.number_baselines)) == NULL) {
            puts("The 'IA' array isn't created!");
            exit(-1);
        }
        if ((ISO = (int *) malloc(sizeof(int) * gs.number_baselines)) == NULL) {
            puts("The 'ISO' array isn't created!");
            exit(-1);
        }

        for (int i = 0; i < gs.number_baselines; i++)
            CL[i] = atoi(read_line(caBuf, ciBufSize, fInput));
            
        for (int i = 0; i < gs.number_baselines; i++)
            strcpy(baseline_file_names[i], read_line(caBuf, ciBufSize, fInput));

        //read initial states filenames
        for (int i = 0; i < gs.number_baselines; i++)
            strcpy(statedat_file_names[i], read_line(caBuf, ciBufSize, fInput));

        for (int i = 0; i < gs.number_baselines; i++)
            IA[i] = atof(read_line(caBuf, ciBufSize, fInput));

        for (int i = 0; i < gs.number_baselines; i++)
            ISO[i] = atoi(read_line(caBuf, ciBufSize, fInput));

        gs.elites = atoi(read_line(caBuf, ciBufSize, fInput));
        gs.autosave = atoi(read_line(caBuf, ciBufSize, fInput));
        gs.recording_frequency = atoi(read_line(caBuf, ciBufSize, fInput));

        printf("Number of organisms: %d\n", gs.number_organisms);
        printf("Number of optimized parameters: %d\n", gs.number_genes);
        printf("Number of generations: %d\n", gs.generations);
        printf("Number of input baselines: %d\n", gs.number_baselines);
        printf("Number of elite organisms: %d\n", gs.elites);
        printf("Number of CPUs: %d\n", size);


        MPI_Bcast(&gs, 1, GlobalSetupMPI, 0, MPI_COMM_WORLD);
        //end master
    } else {
		//slave receives params from master thats all

        MPI_Bcast(&gs, 1, GlobalSetupMPI, 0, MPI_COMM_WORLD);

        if ((CL = (int *) malloc(sizeof(int) * gs.number_baselines)) == NULL) {
            puts("The 'CL' array isn't created!");
            exit(-1);
        }
        if ((IA = (double *) malloc(sizeof(double) * gs.number_baselines)) == NULL) {
            puts("The 'IA' array isn't created!");
            exit(-1);
        }
        if ((ISO = (int *) malloc(sizeof(int) * gs.number_baselines)) == NULL) {
            puts("The 'ISO' array isn't created!");
            exit(-1);
        }
    }//end slave

    MPI_Bcast(IA, gs.number_baselines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(CL, gs.number_baselines, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ISO, gs.number_baselines, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(baseline_file_names, gs.number_baselines * 256, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(statedat_file_names, gs.number_baselines * 256, MPI_CHAR, 0, MPI_COMM_WORLD);





    int *TIME; //?? array of durations of each AP record?
    
    if ((TIME = (int *) malloc(sizeof(int) * gs.number_baselines)) == NULL) {
        puts("The 'TIME' array isn't created!");
        exit(-1);
    }



    if (rank == 0) {
        if ((SD = (double *) malloc(sizeof(double) * gs.number_organisms)) == NULL) {
            puts("The 'SD' array isn't created!");
            exit(-1);
        }
        if ((best_scaling_factor = (float *) malloc(sizeof(float) * gs.number_organisms * gs.number_baselines)) ==
            NULL) {
            puts("The 'best_scaling_factor' array isn't created!");
            exit(-1);
        }
        if ((best_scaling_shift = (float *) malloc(sizeof(float) * gs.number_organisms * gs.number_baselines)) ==
            NULL) {
            puts("The 'best_scaling_shift' array isn't created!");
            exit(-1);
        }
        if ((next_generation = (double *) malloc(sizeof(double) * gs.number_organisms * gs.number_genes)) == NULL) {
            puts("The 'next_generation' array isn't created!");
            exit(-1);
        }
        if ((after_cross = (double *) malloc(sizeof(double) * gs.number_organisms * gs.number_genes)) == NULL) {
            puts("The 'after_cross' array isn't created!");
            exit(-1);
        }
        if ((after_mut = (double *) malloc(sizeof(double) * gs.number_organisms * gs.number_genes)) == NULL) {
            puts("The 'after mut' array isn't created!");
            exit(-1);
        }
    } else {//duplicating...
        if ((SD = (double *) malloc(sizeof(double) * gs.number_organisms / size)) == NULL) {
            puts("The 'SD' array isn't created!");
            exit(-1);
        }
        if ((best_scaling_factor = (float *) malloc(
                sizeof(float) * gs.number_organisms / size * gs.number_baselines)) == NULL) {
            puts("The 'best_scaling_factor' array isn't created!");
            exit(-1);
        }
        if ((best_scaling_shift = (float *) malloc(sizeof(float) * gs.number_organisms / size * gs.number_baselines)) ==
            NULL) {
            puts("The 'best_scaling_shift' array isn't created!");
            exit(-1);
        }
        if ((next_generation = (double *) malloc(sizeof(double) * gs.number_organisms / size * gs.number_genes)) ==
            NULL) {
            puts("The 'next_generation' array isn't created!");
            exit(-1);
        }
        if ((after_cross = (double *) malloc(sizeof(double) * gs.number_organisms / size * gs.number_genes)) == NULL) {
            puts("The 'after_cross' array isn't created!");
            exit(-1);
        }
        if ((after_mut = (double *) malloc(sizeof(double) * gs.number_organisms / size * gs.number_genes)) == NULL) {
            puts("The 'after mut' array isn't created!");
            exit(-1);
        }

    }





    if (rank == 0) {
        time_sum = 0;
        for (baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
            FILE *filename = fopen(baseline_file_names[baseline_counter], "r");
            if (filename == NULL) {
                printf("Cannot open baseline file: %s\n", baseline_file_names[baseline_counter]);
                return -1;
            }
            TIME[baseline_counter] = time_array(time_sum, filename);
            fclose(filename);
        }
    }
    MPI_Bcast(&time_sum, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(TIME, gs.number_baselines, MPI_INT, 0, MPI_COMM_WORLD);




    if ((AP_control = (double *) malloc(sizeof(double) * (time_sum))) == NULL) {
        puts("The 'AP_control' array isn't created!");
        exit(-1);
    }
    if ((AP_current = (double *) malloc(sizeof(double) * gs.number_organisms * (time_sum))) == NULL) {
        puts("The 'AP_current' array isn't created!");
        exit(-1);
    }







    /*Read initial states for each BASELINE from the files in the directory.*/
    State initial_state[gs.number_baselines];
    State *state_struct, *elite_state, *state_struct_rewrite;

    if (rank == 0) {
        state_struct = (struct State *) malloc(sizeof(struct State) * gs.number_organisms * (gs.number_baselines));
        state_struct_rewrite = (struct State *) malloc(
                sizeof(struct State) * gs.number_organisms * (gs.number_baselines));
    } else {
        state_struct = (struct State *) malloc(
                sizeof(struct State) * gs.number_organisms / size * (gs.number_baselines));
        state_struct_rewrite = (struct State *) malloc(
                sizeof(struct State) * gs.number_organisms / size * (gs.number_baselines));
    }

    elite_state = (struct State *) malloc(sizeof(struct State) * gs.elites * (gs.number_baselines));


   
    if (rank == 0) {
        for (int i = 0; i < gs.number_baselines; i++) {
            FILE *fin = fopen(statedat_file_names[i], "r");
            if (!fin) {
                printf("Cannot open IC file: %s\n", statedat_file_names[i]);
                return -1;
            }
            fread(&initial_state[i], sizeof(struct State), 1, fin); //?? NOT SAFE!
            fclose(fin);
        }
    }
    MPI_Bcast(initial_state, gs.number_baselines, StateVectorMPI, 0, MPI_COMM_WORLD);
    
    


    if (rank == 0) { //?? why only root got initial state in state_struct?
        for (int j = 0; j < gs.number_organisms; j++)
            for (int i = 0; i < gs.number_baselines; i++)
                state_struct[i + j *
                                 gs.number_baselines] = initial_state[i];
            //use same structure for every organism initially
    }






    //read baseline APs
    if (rank == 0) {
        for (t_current = 0, baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
            FILE *file = fopen(baseline_file_names[baseline_counter], "r");
            if (!file) {
                printf("Cannot open baseline file: %s\n", baseline_file_names[baseline_counter]);
                exit(-1);
            }
            scanf_baseline(0, TIME[baseline_counter], file, &AP_control[t_current]);
    //        normalize_baseline(0, TIME[baseline_counter], &AP_control[t_current]); //rescale baseline from 0 to 1. Max is the time of maximum for baseline.
            fclose(file);
            t_current += TIME[baseline_counter];
        }
    }
    MPI_Bcast(AP_control, time_sum, MPI_DOUBLE, 0, MPI_COMM_WORLD);




    if (rank == 0) //?? initial population only at root?
        initial_population(next_generation, left_border, right_border, initial_state, gs.number_organisms,
                           gs.number_genes, gs.number_baselines, gs.autosave);





    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            MPI_Isend(&next_generation[gs.number_genes * gs.number_organisms / size * (i)],
                      gs.number_organisms * gs.number_genes / size, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &reqs[0]);
            MPI_Isend(&state_struct[gs.number_baselines * gs.number_organisms / size * (i)],
                      gs.number_organisms * gs.number_baselines / size, StateVectorMPI, i, 1, MPI_COMM_WORLD, &reqs[1]);
        }

        for (int i = 0; i < gs.number_organisms / size; i++) {
            for (t_current = 0, baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
                action_potential(&state_struct[baseline_counter + i * gs.number_baselines],
                                 &next_generation[i * gs.number_genes], &AP_current[t_current + i * (time_sum)],
                                 CL[baseline_counter], IA[baseline_counter], TIME[baseline_counter],
                                 ISO[baseline_counter], baseline_counter, gs.number_baselines, gs.number_genes);
                t_current += TIME[baseline_counter];
            }
        }
        for (int i = 1; i < size; i++) {
            MPI_Recv(&AP_current[i * (time_sum) * gs.number_organisms / size], (time_sum) * gs.number_organisms / size,
                     MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &stats[1]);
            MPI_Recv(&state_struct[i * gs.number_baselines * gs.number_organisms / size],
                     gs.number_baselines * gs.number_organisms / size, StateVectorMPI, i, 4, MPI_COMM_WORLD, &stats[0]);
            MPI_Recv(&next_generation[i * gs.number_genes * gs.number_organisms / size],
                     gs.number_genes * gs.number_organisms / size, MPI_DOUBLE, i, 8, MPI_COMM_WORLD, &stats[3]);
        }

    } else {
//	printf("line 303, rank %d\n", rank);
        fflush(stdout);
        MPI_Recv(next_generation, gs.number_organisms * gs.number_genes / size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD,
                 &stats[0]);
        MPI_Recv(state_struct, gs.number_organisms * gs.number_baselines / size, StateVectorMPI, 0, 1, MPI_COMM_WORLD,
                 &stats[1]);
        for (int i = 0; i < gs.number_organisms / size; i++) {
            for (t_current = 0, baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
                action_potential(&state_struct[baseline_counter + i * gs.number_baselines],
                                 &next_generation[i * gs.number_genes], &AP_current[t_current + i * (time_sum)],
                                 CL[baseline_counter], IA[baseline_counter], TIME[baseline_counter],
                                 ISO[baseline_counter], baseline_counter, gs.number_baselines, gs.number_genes);
                t_current += TIME[baseline_counter];
            }
        }
        MPI_Send(AP_current, (time_sum) * gs.number_organisms / size, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        MPI_Send(state_struct, gs.number_baselines * gs.number_organisms / size, StateVectorMPI, 0, 4, MPI_COMM_WORLD);
        MPI_Send(next_generation, gs.number_genes * gs.number_organisms / size, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
    }





    //start GA?
    for (int cntr = 0; cntr < gs.generations; cntr++) {
        if (rank == 0) {
            printf("\nGeneration = %d\n", cntr);

            int SD_index[gs.number_organisms];

            fitness_function(AP_control, AP_current, best_scaling_factor, best_scaling_shift, TIME, SD, SD_index,
                             gs.number_organisms, gs.number_baselines, time_sum);

            /*Save final state of best organism to state.dat to get closer to steady state in the next run.*/
            double elite_array[gs.elites];
            int elite_index_array[gs.elites];
            int cc = 0;
            int c = 0;

            while (c < gs.elites) {
                if ((cc == 0) || ((cc != 0) && (SD[cc] != SD[cc - 1]))) {
                    elite_array[c] = SD[cc];
                    elite_index_array[c] = SD_index[cc];
                    for (int i = 0; i < gs.number_baselines; i++)
                        elite_state[c * gs.number_baselines + i] = state_struct[
                                elite_index_array[c] * gs.number_baselines + i];
                    c += 1;
                }
                cc += 1;
                if (cc == gs.number_organisms - 1) {
                    while (c < gs.elites) {
                        elite_array[c] = SD[cc];
                        elite_index_array[c] = SD_index[cc];
                        for (int i = 0; i < gs.number_baselines; i++)
                            elite_state[c * gs.number_baselines + i] = state_struct[
                                    elite_index_array[c] * gs.number_baselines + i];
                        c += 1;
                        cc -= 1;
                    }
                }
            }

            double *elite_organisms;
            elite_organisms = (double *) malloc(sizeof(double) * gs.elites * gs.number_genes);
            for (int i = 0; i < gs.elites; i++)
                for (int j = 0; j < gs.number_genes; j++)
                    elite_organisms[i * gs.number_genes + j] = next_generation[elite_index_array[i] * gs.number_genes +
                                                                               j];

            double average = 0;
            for (int c = 0; c < gs.number_organisms; c++) average += SD[c];
            average = average / gs.number_organisms;


            fflush(stdout);
            printf("Best SD: %.4f mV\n", SD[0]);
            printf("Average SD: %.4f mV\n", average);
            printf("The fittest organism:\n");
            for (int i = 0; i < gs.number_genes - 3 * gs.number_baselines; i++) {
                printf("%g ", next_generation[elite_index_array[0] * gs.number_genes + i]);
            }
            printf("\n");
            printf("Na_i: ");
            for (int i = gs.number_genes - 3 * gs.number_baselines; i < gs.number_genes - 2 * gs.number_baselines; i++) {
                printf("%g ", next_generation[elite_index_array[0] * gs.number_genes + i]);
            }
            printf("\n");
            printf("Ca_sr: ");
            for (int i = gs.number_genes - 2 * gs.number_baselines; i < gs.number_genes - 1 * gs.number_baselines; i++) {
                printf("%g ", next_generation[elite_index_array[0] * gs.number_genes + i]);
            }
            printf("\n");
            printf("K_i: ");
            for (int i = gs.number_genes - 1 * gs.number_baselines; i < gs.number_genes - 0 * gs.number_baselines; i++) {
                printf("%g ", next_generation[elite_index_array[0] * gs.number_genes + i]);
            }
            printf("\n");

            printf("\n\n");
            fflush(stdout);

            //scaling_factor = best_scaling_factor[baseline_counter+gs.number_baselines*elite_index_array[0]];
            //scaling_shift = best_scaling_shift[baseline_counter+gs.number_baselines*elite_index_array[0]];

            writing_to_output_files(best, avr, owle, ctrl_point, text, sd, ap_best, SD, SD_index, average,
                                    &next_generation[elite_index_array[0] * gs.number_genes], gs.number_genes,
                                    gs.number_organisms, next_generation, cntr, gs.recording_frequency,
                                    gs.number_baselines, elite_state, CL, AP_current, AP_control, best_scaling_factor,
                                    best_scaling_shift, elite_index_array[0], TIME, elite_index_array[0] * (time_sum));

            /*Genetic Operators*/
            int mpool[gs.number_organisms];
            tournament_selection(mpool, SD, state_struct, state_struct_rewrite, gs.number_organisms,
                                 gs.number_baselines);
            sbx_crossover(next_generation, after_cross, mpool, left_border, right_border, gs.number_organisms,
                          gs.number_genes);



            double after_cross_normalized[gs.number_organisms * gs.number_genes];
            double left_border_normalized[gs.number_genes], right_border_normalized[gs.number_genes];
            int number_conductancies = 15;

            normalize_genes(/*in*/ after_cross, left_border, right_border,
                                   gs.number_organisms, gs.number_genes, number_conductancies,
                            /*out*/ after_cross_normalized, left_border_normalized, right_border_normalized);

            double after_mut_normalized[gs.number_organisms * gs.number_genes];
            cauchy_mutation(after_mut_normalized, after_cross_normalized,
                            left_border_normalized, right_border_normalized,
                            gs.number_organisms, gs.number_genes);

            denormalize_genes(/*in*/ after_mut_normalized, left_border, right_border,
                                   gs.number_organisms, gs.number_genes, number_conductancies,
                              /*out*/ after_mut);


            /*SLOWEST: AP calculations*/
            t_current = 0;
            for (int i = 1; i < size; i++) {
                MPI_Isend(&after_mut[gs.number_organisms * gs.number_genes / size * (i)],
                          gs.number_organisms * gs.number_genes / size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &reqs[0]);
                MPI_Isend(&state_struct[gs.number_baselines * gs.number_organisms / size * (i)],
                          gs.number_organisms * gs.number_baselines / size, StateVectorMPI, i, 2, MPI_COMM_WORLD,
                          &reqs[1]);
            }
            for (int i = 0; i < gs.number_organisms / size; i++) {
                for (baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
                    action_potential(&state_struct[baseline_counter + i * gs.number_baselines],
                                     &after_mut[i * gs.number_genes], &AP_current[t_current + i * (time_sum)],
                                     CL[baseline_counter], IA[baseline_counter], TIME[baseline_counter],
                                     ISO[baseline_counter], baseline_counter, gs.number_baselines, gs.number_genes);
                    t_current += TIME[baseline_counter];

                }
                t_current = 0;
            }

            for (int i = 1; i < size; i++) {
                MPI_Recv(&AP_current[i * (time_sum) * gs.number_organisms / size],
                         (time_sum) * gs.number_organisms / size, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &stats[1]);
                MPI_Recv(&state_struct[i * gs.number_baselines * gs.number_organisms / size],
                         gs.number_baselines * gs.number_organisms / size, StateVectorMPI, i, 3, MPI_COMM_WORLD,
                         &stats[0]);
                MPI_Recv(&after_mut[i * gs.number_genes * gs.number_organisms / size],
                         gs.number_genes * gs.number_organisms / size, MPI_DOUBLE, i, 8, MPI_COMM_WORLD, &stats[3]);
            }

            fitness_function(AP_control, AP_current, best_scaling_factor, best_scaling_shift, TIME, SD, SD_index,
                             gs.number_organisms, gs.number_baselines, time_sum);

            fflush(stdout);

            /*Elite AP recalculations to get closer to steady state. Note that we suppose number of cores to be more then the number of Elite organisms in this implementation.
             * 1 Elite per core.*/
            for (int h = 1; h < gs.elites; h++) {
                MPI_Isend(&elite_organisms[h * gs.number_genes], gs.number_genes, MPI_DOUBLE, h, 4, MPI_COMM_WORLD,
                          &reqs[0]);
                MPI_Isend(&elite_state[gs.number_baselines * h], gs.number_baselines, StateVectorMPI, h, 6,
                          MPI_COMM_WORLD, &reqs[1]);
            }

            t_current = 0;
            for (baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
                action_potential(&elite_state[baseline_counter], &elite_organisms[0],
                                 &AP_current[SD_index[gs.number_organisms - gs.elites] * (time_sum) + t_current],
                                 CL[baseline_counter], IA[baseline_counter], TIME[baseline_counter],
                                 ISO[baseline_counter], baseline_counter, gs.number_baselines, gs.number_genes);
                t_current += TIME[baseline_counter];
            }

            for (int h = 1; h < gs.elites; h++) {
                MPI_Recv(&AP_current[SD_index[gs.number_organisms - gs.elites + h] * (time_sum)], time_sum, MPI_DOUBLE,
                         h, 5, MPI_COMM_WORLD, &stats[1]);
                MPI_Recv(&elite_state[h * gs.number_baselines], gs.number_baselines, StateVectorMPI, h, 7,
                         MPI_COMM_WORLD, &stats[2]);
                MPI_Recv(&elite_organisms[h * gs.number_genes], gs.number_genes, MPI_DOUBLE, h, 9, MPI_COMM_WORLD,
                         &stats[4]);
            }

            /*Replace worst organisms to elite*/
            int num = 0;
            for (int i = gs.number_organisms - gs.elites; i < gs.number_organisms; i++) {
                if (SD[i] > elite_array[num]) {
                    for (int j = 0; j < gs.number_genes; j++)
                        after_mut[SD_index[i] * gs.number_genes + j] = elite_organisms[num * gs.number_genes + j];
                    for (int j = 0; j < gs.number_baselines; j++)
                        state_struct[SD_index[i] * gs.number_baselines + j] = elite_state[num * gs.number_baselines +
                                                                                          j];
                    num += 1;
                }
            }

            for (int i = 0; i < gs.number_organisms; i++)
                for (int j = 0; j < gs.number_genes; j++)
                    next_generation[i * gs.number_genes + j] = after_mut[i * gs.number_genes + j];

            free(elite_organisms);
        } else {
//    char hostname[256];
//    gethostname(hostname, sizeof(hostname));
//    printf("PID %d on %s ready for attach\n", getpid(), hostname);
            //   fflush(stdout);
//    while (0 == i)
//        sleep(5);

//	printf("line 459, rank %d\n", rank);
            fflush(stdout);

            MPI_Recv(after_mut, gs.number_organisms * gs.number_genes / size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
                     &stats[0]);

//	printf("line 473, rank %d\n", rank);
            MPI_Recv(state_struct, gs.number_organisms * gs.number_baselines / size, StateVectorMPI, 0, 2,
                     MPI_COMM_WORLD, &stats[1]);

//	printf("line 476, rank %d\n", rank);
            t_current = 0;
            for (int i = 0; i < gs.number_organisms / size; i++) {
                for (baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
                    action_potential(&state_struct[baseline_counter + i * gs.number_baselines],
                                     &after_mut[i * gs.number_genes], &AP_current[t_current + i * (time_sum)],
                                     CL[baseline_counter], IA[baseline_counter], TIME[baseline_counter],
                                     ISO[baseline_counter], baseline_counter, gs.number_baselines, gs.number_genes);
                    t_current += TIME[baseline_counter];
                }
                t_current = 0;
            }

            MPI_Send(AP_current, (time_sum) * gs.number_organisms / size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            MPI_Send(state_struct, gs.number_baselines * gs.number_organisms / size, StateVectorMPI, 0, 3,
                     MPI_COMM_WORLD);
            MPI_Send(after_mut, gs.number_genes * gs.number_organisms / size, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);

            if (rank < gs.elites) {

//	printf("line 481, rank %d\n", rank);
                fflush(stdout);
                MPI_Recv(after_mut, gs.number_genes, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &stats[0]);
                MPI_Recv(elite_state, gs.number_baselines, StateVectorMPI, 0, 6, MPI_COMM_WORLD, &stats[1]);
                t_current = 0;
                for (baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
                    action_potential(&elite_state[baseline_counter], &after_mut[0], &AP_current[t_current],
                                     CL[baseline_counter], IA[baseline_counter], TIME[baseline_counter],
                                     ISO[baseline_counter], baseline_counter, gs.number_baselines, gs.number_genes);
                    t_current += TIME[baseline_counter];
                }


                MPI_Send(AP_current, time_sum, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
                MPI_Send(elite_state, gs.number_baselines, StateVectorMPI, 0, 7, MPI_COMM_WORLD);
                MPI_Send(after_mut, gs.number_genes, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
            }
        }
    }









    //finalize
    if (rank == 0) {
        fclose(text);
        fclose(owle);
        fclose(avr);
        fclose(best);
        fclose(ap_best);
        fclose(sd);
    }

    end = MPI_Wtime();
    if (rank == 0)
		printf("Time: %f sec\n", end - start1);

    MPI_Type_free(&GlobalSetupMPI);
    MPI_Type_free(&StateVectorMPI);

    free(left_border);
    free(right_border);
    free(CL);
    free(IA);
    free(ISO);
    free(SD);
    free(TIME);
    free(best_scaling_factor);
    free(best_scaling_shift);
    free(next_generation);
    free(after_cross);
    free(after_mut);
    free(AP_control);
    free(AP_current);
    free(state_struct);
    free(state_struct_rewrite);
    free(elite_state);

    
    MPI_Finalize();
    return 0;
}
