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

#include <omp.h>
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
    srand(time(NULL));
    
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





   
    long time_sum;
    


    double *AP_control, *AP_current, *SD, *next_generation, *after_mut, *after_cross;//, *Na_conc;
    double  scaling_factor, scaling_shift;
    float *best_scaling_factor, *best_scaling_shift;



    //open files and make sure it was successful
    int init_status = 0;

    FILE *f = 0, *owle = 0, *best = 0, *avr = 0, *test = 0, *ap_best = 0, *text = 0, *sd = 0, *ctrl_point = 0;

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
    const double start_time = MPI_Wtime();




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
        printf("Number of MPI nodes: %d\n", size);
        printf("Number of cores at root: %d\n", omp_get_num_threads());

        if ((gs.number_organisms - gs.elites) % 2 != 0) {
            printf("Even number of mutants required!\nAbort");
            return -1;
        }
        
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

    //we will broadcast the content of these file from the root
    //MPI_Bcast(baseline_file_names, gs.number_baselines * 256, MPI_CHAR, 0, MPI_COMM_WORLD);
    //MPI_Bcast(statedat_file_names, gs.number_baselines * 256, MPI_CHAR, 0, MPI_COMM_WORLD);





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
        for (int baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
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
    State *state_struct, *state_struct_rewrite;

    


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
    

    //read baseline APs
    if (rank == 0) {
        for (long t_current = 0, baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
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




    //very first step of GA: random gen of initial parameters (or read from a file)
    if (rank == 0) {
        for (int j = 0; j < gs.number_organisms; j++)
            for (int i = 0; i < gs.number_baselines; i++)
                state_struct[i + j * gs.number_baselines] = initial_state[i];
            //use same structure for every organism initially

        initial_population(next_generation, left_border, right_border, initial_state, gs.number_organisms,
                           gs.number_genes, gs.number_baselines, gs.autosave);
    }
    
    //main GA cycle
    for (int cntr = 0; cntr < gs.generations; cntr++) {

        MPI_Scatter(next_generation, gs.number_organisms * gs.number_genes / size,
            MPI_DOUBLE, next_generation, gs.number_organisms * gs.number_genes / size,
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
        MPI_Scatter(state_struct, gs.number_baselines * gs.number_organisms / size,
            StateVectorMPI, state_struct, gs.number_baselines * gs.number_organisms / size,
            StateVectorMPI, 0, MPI_COMM_WORLD);

        #pragma omp parallel for
        for (int i = 0; i < gs.number_organisms / size; i++) {
            for (long t_current = 0, baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
                action_potential(&state_struct[baseline_counter + i * gs.number_baselines],
                                 &next_generation[i * gs.number_genes], &AP_current[t_current + i * (time_sum)],
                                 CL[baseline_counter], IA[baseline_counter], TIME[baseline_counter],
                                 ISO[baseline_counter], baseline_counter, gs.number_baselines, gs.number_genes);
                t_current += TIME[baseline_counter];
            }
        }

        MPI_Gather(AP_current, (time_sum) * gs.number_organisms / size, MPI_DOUBLE,
            AP_current, (time_sum) * gs.number_organisms / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        MPI_Gather(state_struct, gs.number_baselines * gs.number_organisms / size, StateVectorMPI,
            state_struct, gs.number_baselines * gs.number_organisms / size, StateVectorMPI, 0, MPI_COMM_WORLD);
    
        MPI_Gather(next_generation, gs.number_genes * gs.number_organisms / size, MPI_DOUBLE, next_generation,
            gs.number_genes * gs.number_organisms / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);



        if (rank == 0) {
            printf("\nGeneration %d\n", cntr);
            
            State *elite_state = state_struct;
            State *mutant_state = state_struct + gs.elites;
            const int mutant_number = gs.number_organisms - gs.elites;


            //store pairs of (error, index in state_struct) sorted by error in increasing order
            //thus, first elements are for elite organisms
            std::vector<std::pair<double, int>> sd_n_index(gs.number_organisms);
            

            //TODO fitness_function should be evaluated at each node
            //And there is not need to send AP to the root
            fitness_function(AP_control, AP_current, best_scaling_factor, best_scaling_shift, TIME, sd_n_index,
                             gs.number_organisms, gs.number_baselines, time_sum);

            assert(sd_n_index[0].first < sd_n_index.back().first);
            
            //save elites in elite buffer and later copy it to main array
            State buf_elite_state[gs.elites * gs.number_baselines];
            double * buf_elite_genes = (double *) malloc(sizeof(double) * gs.elites * gs.number_genes);
            for (int i = 0; i < gs.elites; i++) {
                const int elite_index = sd_n_index[i].second;
                for (int j = 0; j < gs.number_baselines; j++)
                    buf_elite_state[i * gs.number_baselines + j] = 
                        state_struct[elite_index * gs.number_baselines + j];
                        
                for (int j = 0; j < gs.number_genes; j++)
                    buf_elite_genes[i * gs.number_genes + j] =
                        next_generation[elite_index * gs.number_genes + j];
            }
            

            //print some output to see that the algo converges
            double average_error = 0;
            for (int i = 0; i < gs.number_organisms; i++)
                average_error += sd_n_index[i].first;
            average_error /= gs.number_organisms;
            
            fflush(stdout);
            printf("Best SD: %.4f mV\n", sd_n_index[0].first);
            printf("Worst SD: %.4f mV\n", sd_n_index.back().first);
            printf("Average SD: %.4f mV\n", average_error);
            
            

            
            /*turn it off for a while
            //Save final state of best organism to state.dat to get closer to steady state in the next run.

            
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
*/


            
            //create mating pool
            
            /*Genetic Operators for mutants*/
            State buf_mutant_state[mutant_number * gs.number_baselines];
            double * buf_mutant_genes = (double *) malloc(sizeof(double) * mutant_number * gs.number_genes);
            
            int mpool[mutant_number]; //mpool: mutant_index -> next_generation_index
            //so let's find it
            tournament_selection(mpool, mutant_number, sd_n_index, gs.elites);

            //we also need to shuffle mpool because it is sorted!
            for (int i = 0; i < mutant_number; i++) {
                int j = i + rand() % (mutant_number - i);
                std::swap(mpool[i], mpool[j]);
            }
            
            //only now copy states from next_generation to buf_mutant_state according to the shuffled mpool!
            for (int i = 0; i < mutant_number; i++) {            
                for (int j = 0; j < gs.number_baselines; j++)
                    buf_mutant_state[i * gs.number_baselines + j] = state_struct[mpool[i] * gs.number_baselines + j];

            }

            //ATTENTION!!!! after_cross is a buffer of size gs.number_organisms!!!! But we use it only for mutants
            sbx_crossover(next_generation, after_cross, mpool, left_border, right_border, mutant_number,
                          gs.number_genes);

            //no need of mpool anymore
            
            

            
            double mutants_after_cross_normalized[mutant_number * gs.number_genes];

            double left_border_normalized[gs.number_genes], right_border_normalized[gs.number_genes];

            const int number_conductancies = 15;

            normalize_genes(/*in*/ after_cross, left_border, right_border,
                                   mutant_number, gs.number_genes, number_conductancies,
                            /*out*/ mutants_after_cross_normalized, left_border_normalized, right_border_normalized);

            double after_mut_normalized[mutant_number * gs.number_genes];
            cauchy_mutation(mutants_after_cross_normalized, after_mut_normalized,
                            left_border_normalized, right_border_normalized,
                            mutant_number, gs.number_genes);

            denormalize_genes(/*in*/ after_mut_normalized, left_border, right_border,
                                   mutant_number, gs.number_genes, number_conductancies,
                              /*out*/ buf_mutant_genes);
                              
            /*Genetic Operators for mutants DONE*/

            //now copy elite to main arrays
            for (int i = 0; i < gs.elites; i++) {
                for (int j = 0; j < gs.number_baselines; j++)
                        state_struct[i * gs.number_baselines + j] =
                            buf_elite_state[i * gs.number_baselines + j];
                        
                for (int j = 0; j < gs.number_genes; j++)
                        next_generation[i * gs.number_genes + j] = 
                            buf_elite_genes[i * gs.number_genes + j];
            }
            //and copy mutants to main arrays
            for (int i = 0; i < mutant_number; i++) {
                const int index = i + gs.elites;
                for (int j = 0; j < gs.number_baselines; j++)
                        state_struct[index * gs.number_baselines + j] =
                            buf_mutant_state[i * gs.number_baselines + j];
                        
                for (int j = 0; j < gs.number_genes; j++)
                        next_generation[index * gs.number_genes + j] = 
                            buf_mutant_genes[i * gs.number_genes + j];
            }

            free(buf_elite_genes);
            free(buf_mutant_genes);
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

    const double end_time = MPI_Wtime();
    if (rank == 0)
		printf("Time: %f sec\n", end_time - start_time);

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
    //free(elite_state);

    
    MPI_Finalize();
    return 0;
}
