// Main GA module
//
// Genetic Algorithm implementation
//
// Authors:
//     Dmitrii Smirnov <dmitrii.smirnov@phystech.edu>
//     Roman Syunyaev <roman.syunyaev@gmail.com>
//     Alexander Timofeev <richardstallman42@gmail.com>
//     Andrey Pikunov <pikunov@phystech.edu>
//     Laboratory of Human Physiology, Moscow Institute of Physics and Technology, 2020
//
// License:
//     Redistribution and use of the code with or without modification
//     are permitted for any scientific, research and educational purposes.
//
//     Redistribution and use of the code with or without modification
//     are prohibited for any commercial purposes.


#include <mpi.h>
#include <iostream>
#include "pcg_random.hpp"

#include "basic_population.h"
#include "genetic_algorithm.h"
#include "test_functions.h"
#include "polynomial_mutation.h"
#include "tournament_selection.h"
#include "sbx_crossover.h"
#include "cauchy_mutation.h"

#include "maleckar_model.h"
#include "cellml_ode_solver.h"
#include "optimization_problem.h"
#include <json.hpp>

#ifdef HIDE_CODE
class Population
{
    void mutation()
    {
/* Genes transformation */
            const double gamma = 1.;
            // Suppose we have X distributed according to cauchy distribution with given gamma and zero mu
            // if we want Q_90 of X to be equal some fixed value x
            // we should set gamma equal to 0.15 * x

            const int number_conc_types = 3;
            const int genes_without_concentrations = gs.number_genes - number_conc_types * gs.number_baselines;
            double array_gamma[gs.number_genes];

            for (int i = 0; i < genes_without_concentrations; ++i) {
                array_gamma[i] = 0.015; // Q_90 = 0.1 [log10_units], multipliers
            }
            for (int i = 0; i < gs.number_baselines; ++i) {
                array_gamma[genes_without_concentrations + i + 0 * gs.number_baselines] = 0.15; // Q_90 = 1 mM, Na_i
                array_gamma[genes_without_concentrations + i + 1 * gs.number_baselines] = 0.015; // Q_90 = 0.1 mM, Ca_rel
                array_gamma[genes_without_concentrations + i + 2 * gs.number_baselines] = 0.15; // Q_90 = 1 mM, K_i
            }

            const double scaler_dimensional = 1 / sqrt(gs.number_genes);
            // is approximately mean projection length of the unit vector onto some direction in `gs.number_genes`-dimensional space





            double left_border_transformed[gs.number_genes];
            double right_border_transformed[gs.number_genes];
            for (int i_genes = 0; i_genes < gs.number_genes; ++i_genes) {
                left_border_transformed[i_genes] = 0;
                if (i_genes < genes_without_concentrations) {
                    right_border_transformed[i_genes] =
                            log10(right_border[i_genes] / left_border[i_genes]) * gamma /
                            (array_gamma[i_genes] / scaler_dimensional);
                } else {
                    right_border_transformed[i_genes] =
                            (right_border[i_genes] - left_border[i_genes]) * gamma /
                            (array_gamma[i_genes] / scaler_dimensional);
                }
            }





            transform_genes(/*in*/ after_cross, left_border, right_border,
                                   gs.number_mutants, gs.number_genes, genes_without_concentrations,
                                   left_border_transformed, right_border_transformed,
                    /*out*/ genes_mutant_after_cross_transformed);

            cauchy_mutation(genes_mutant_after_mut_transformed, genes_mutant_after_cross_transformed,
                            left_border_transformed, right_border_transformed,
                            gs.number_mutants, gs.number_genes, gamma);

            transform_genes_back(/*in*/ genes_mutant_after_mut_transformed, left_border_transformed,
                                        right_border_transformed,
                                        gs.number_mutants, gs.number_genes, genes_without_concentrations,
                                        left_border, right_border,
                    /*out*/ buf_mutant_genes);

    }
};


void normalize_baseline(int j0, int j1, double *AP_control)
{
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




void old_code(int argc, char *argv[])
{
    long time_sum;
    double *AP_control, *AP_current, *SD, *next_generation;//, *Na_conc;
    double scaling_factor, scaling_shift;
    float *best_scaling_factor, *best_scaling_shift;




    //open files and make sure it was successful
    int init_status = 0;

    FILE *file_ap_best = 0, *file_dump = 0;

    if (rank == 0) {
        std::filesystem::create_directory("ga_output");
        file_ap_best = fopen("./ga_output/ap.bin", "wb");
        file_dump = fopen("./ga_output/dump.bin", "wb");

        if (!(file_ap_best && file_dump)) {
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





    double *IA = 0; //Stimulation current, A/F
    double *right_border = 0, *left_border = 0; //min and max possible values of parameters
    int *CL = 0; //Cycle lengths

    //1 - include ISO model
    //0 - without ISO model
    int *ISO = 0;

    GlobalSetup gs;

    char baseline_file_names[256][256], statedat_file_names[256][256];

    if (rank == 0) {
        const int ciBufSize = 1024;
        char caBuf[ciBufSize];
        FILE *fInput = fopen(argv[1], "r");
        if (!fInput) {
            printf("Cannot open GA input file!\n");
            return -1;
        }

        gs.number_organisms = atoi(read_line(caBuf, ciBufSize, fInput));
        gs.number_genes = atoi(read_line(caBuf, ciBufSize, fInput));
        gs.number_generations = atoi(read_line(caBuf, ciBufSize, fInput));

        //Parameters ranges
        try {
            left_border = new double [gs.number_genes];
            right_border = new double [gs.number_genes];
        } catch (const std::bad_alloc & e) {
            std::cerr << "Allocation failed: " << e.what() << '\n';
            exit(-1);
        }


        for (int i = 0; i < gs.number_genes; i++) {
            char *token, *ptoken;
            token = read_line(caBuf, ciBufSize, fInput);

            if ((ptoken = strtok(token, " \t\n\r")) == NULL) {
                fclose(fInput);
                std::cerr << "Error in the input file!" << std::endl;
                exit(-1);
            }
            left_border[i] = atof(ptoken);

            if ((ptoken = strtok(NULL, " \t\n\r")) == NULL) {
                fclose(fInput);
                std::cerr << "Error in the input file!\n" << std::endl;
                exit(-1);
            }
            right_border[i] = atof(ptoken);
        }

        gs.number_baselines = atoi(read_line(caBuf, ciBufSize, fInput));

        try {
            CL = new int [gs.number_baselines];
            IA = new double [gs.number_baselines];
            ISO = new int [gs.number_baselines];
        } catch (const std::bad_alloc & e) {
            std::cerr << "Allocation failed: " << e.what() << '\n';
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

        gs.number_elites = atoi(read_line(caBuf, ciBufSize, fInput));
        gs.INIT_FROM_BACKUP_FILE = atoi(read_line(caBuf, ciBufSize, fInput));
        gs.period_backup = atoi(read_line(caBuf, ciBufSize, fInput));


        gs.number_mutants = gs.number_organisms - gs.number_elites;
        if (gs.number_mutants % 2 != 0) {
            printf("Even number of mutants required!\nWe will try to add one more\n");
            gs.number_mutants++;
            gs.number_organisms++;
        }
        if (gs.number_organisms % size != 0) {
            printf("Number of nodes should divide number of organisms completely\n");
            printf("We will add more mutants and elites for that\n");
            int additional_organisms = size - gs.number_organisms % size;
            gs.number_organisms += additional_organisms;
            gs.number_mutants += additional_organisms;
            if (gs.number_mutants % 2 != 0) {
                gs.number_mutants--;
                gs.number_elites++;
            }
        }

        printf("Number of organisms: %d\n", gs.number_organisms);
        printf("Number of optimized parameters: %d\n", gs.number_genes);
        printf("Number of generations: %d\n", gs.number_generations);
        printf("Number of input baselines: %d\n", gs.number_baselines);
        printf("Number of elite organisms: %d\n", gs.number_elites);
        printf("Number of MPI nodes: %d\n", size);
        printf("Number of cores at root: %d\n", omp_get_max_threads());

        MPI_Bcast(&gs, 1, GlobalSetupMPI, 0, MPI_COMM_WORLD);
        //end master
    } else {
        //slave receives params from master thats all

        MPI_Bcast(&gs, 1, GlobalSetupMPI, 0, MPI_COMM_WORLD);

        try {
            CL = new int [gs.number_baselines];
            IA = new double [gs.number_baselines];
            ISO = new int [gs.number_baselines];
        } catch (const std::bad_alloc & e) {
            std::cerr << "Allocation failed: " << e.what() << '\n';
            exit(-1);
        }
    }//end slave

    MPI_Bcast(IA, gs.number_baselines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(CL, gs.number_baselines, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ISO, gs.number_baselines, MPI_INT, 0, MPI_COMM_WORLD);

    //we will broadcast the content of these file from the root
    //MPI_Bcast(baseline_file_names, gs.number_baselines * 256, MPI_CHAR, 0, MPI_COMM_WORLD);
    //MPI_Bcast(statedat_file_names, gs.number_baselines * 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    State *initial_state, *state_struct, *state_struct_rewrite;

    State * buf_elite_state = 0;

    double * buf_elite_genes = 0;

    State * buf_mutant_state = 0;

    double *buf_mutant_genes = 0;

    double * genes_mutant_after_cross_transformed = 0;
    double * genes_mutant_after_mut_transformed = 0;

    try {
        genes_mutant_after_mut_transformed = new double [gs.number_mutants * gs.number_genes];
        genes_mutant_after_cross_transformed = new double [gs.number_mutants * gs.number_genes];
        buf_mutant_state = new State [gs.number_mutants * gs.number_baselines];
        buf_mutant_genes = new double [gs.number_mutants * gs.number_genes];
        buf_elite_state = new State [gs.number_elites * gs.number_baselines];
        buf_elite_genes = new double [gs.number_elites * gs.number_genes];
        TIME = new int [gs.number_baselines];
        initial_state = new State [gs.number_baselines];
        if (rank == 0) {
            state_struct = new State [gs.number_organisms * gs.number_baselines];
            state_struct_rewrite = new State [gs.number_organisms * gs.number_baselines];

            SD = new double [gs.number_organisms];
            best_scaling_factor = new float [gs.number_organisms * gs.number_baselines];
            best_scaling_shift = new float [gs.number_organisms * gs.number_baselines];
            next_generation = new double [gs.number_organisms * gs.number_genes];
            after_cross = new double [gs.number_organisms * gs.number_genes];
        } else {
            state_struct = new State [gs.number_organisms * gs.number_baselines / size];
            state_struct_rewrite = new State [gs.number_organisms * gs.number_baselines / size];

            SD = new double [gs.number_organisms / size];
            best_scaling_factor = new float [gs.number_organisms * gs.number_baselines / size];
            best_scaling_shift = new float [gs.number_organisms * gs.number_baselines / size];
            next_generation = new double [gs.number_organisms * gs.number_genes / size];
            after_cross = new double [gs.number_organisms * gs.number_genes / size];
        }
    } catch (const std::bad_alloc & e) {
        std::cerr << "Allocation failed: " << e.what() << '\n';
        exit(-1);
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

    try {
        AP_control = new double [time_sum];
        AP_current = new double [time_sum * gs.number_organisms];
    } catch (const std::bad_alloc & e) {
        std::cerr << "Allocation failed: " << e.what() << '\n';
        exit(-1);
    }

    if (rank == 0) {
        for (int i = 0; i < gs.number_baselines; i++) {
            FILE *fin = fopen(statedat_file_names[i], "r");
            if (!fin) {
                printf("Cannot open IC file: %s\n", statedat_file_names[i]);
                return -1;
            }

            double a[STATE_ARRAY_SIZE];
            fread(a, sizeof(double), STATE_ARRAY_SIZE, fin);
            array2state(a, &initial_state[i]);
            fclose(fin);
        }
    }
    MPI_Bcast(initial_state, gs.number_baselines, StateVectorMPI, 0, MPI_COMM_WORLD);


    //read baseline APs
    if (rank == 0) {
        for (long t_current = 0, baseline_counter = 0; baseline_counter < gs.number_baselines; baseline_counter++) {
            FILE *file = fopen(baseline_file_names[baseline_counter], "r");
            if (!file) {
                std::cerr << std::string("Cannot open baseline file:") + baseline_file_names[baseline_counter] << std::endl;
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
                           gs.number_genes, gs.number_baselines, gs.INIT_FROM_BACKUP_FILE);
    }
}

#endif




void main_gen_algo(const char *configFilename)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    pcg_extras::seed_seq_from<std::random_device> seed_source;

    std::ifstream configFile;
    configFile.open(configFilename);
    if (!configFile.is_open() && mpi_rank == 0) {
        std::cerr << "Cannot open main config file" << std::endl;
        throw;
    }
    json config;
    configFile >> config;
    configFile.close();

    MaleckarModel model;
    ODESolver solver;
    MinimizeAPbaselines obj;
    ODEoptimization problem(model, solver, obj);

    double time_read_config = MPI_Wtime();
    problem.read_config(configFilename);
    time_read_config = MPI_Wtime() - time_read_config;

    BasicPopulation popMal(problem,
                    config["n_elites"].get<unsigned>(),
                    config["n_organisms"].get<unsigned>());

    double time_population_init = MPI_Wtime();
    popMal.init(pcg64(seed_source));
    time_population_init = MPI_Wtime() - time_population_init;

    if (mpi_rank == 0) {
        std::cout << "time_read_config, s: " << time_read_config << std::endl;
        std::cout << "time_population_init, s: " << time_population_init << std::endl;
    }

    genetic_algorithm(popMal,
            TournamentSelectionFast(pcg64(seed_source)),
            SBXcrossover(pcg64(seed_source), config["crossrate"].get<double>(), config["eta_crossover"].get<int>()),
            PolynomialMutation<pcg64, pcg_extras::seed_seq_from<std::random_device>>
                                            (seed_source,
                                            config["mutrate"].get<double>(),
                                            config["eta_mutation"].get<int>()),
            config["n_generations"].get<unsigned>());
    
    /*
        genetic_algorithm(popMal,
            TournamentSelectionFast(pcg64(seed_source)),
            NoCrossover(),
            NoMutation(),
            100);
    */

    if (mpi_rank == 0) {
        using Results = decltype(problem)::Results;
        using BaselineResult = decltype(problem)::BaselineResult;
        Results results = problem.get_relative_results();
        BaselineResult res = results[0];
        std::cout << "Printing relative to default results" << std::endl;
        for (const auto & cit: res.constantsResult)
            std::cout << cit.first << " " << cit.second << std::endl;
    }
    /*
    if (mpi_rank == 0) {  
        //now try to simply solve ode model
        ODESolver solver;
        //try malecar on one node
        MaleckarModel model;
        std::vector<double> state(model.state_size());
        double * constants = new double [model.constants_size()];
        model.initConsts(constants);
        model.set_constants(constants);
        model.initState(state.data());
        std::vector<double> orig_state = state;
        int is_correct;
        double t0 = 0, start_record = 1, tout = start_record + 1;
        std::vector<double> ap(1000);
        
        double st = MPI_Wtime();
        solver.solve(model, state, is_correct, t0, start_record, tout, ap);
        st = MPI_Wtime() - st;
        std::cout << "TIME: " << st << std::endl;
    
        FILE *state_final = fopen("state_dump", "w");
        fwrite(state.data(), sizeof(double), state.size(), state_final);
        fclose(state_final);
        
        double err = 0;
        std::cout << std::scientific;
        for (int i = 0; i < model.state_size(); i++) {
            std::cout << i << ": " << (state[i] - orig_state[i]) / orig_state[i] * 100 << "%\n";
            err += pow(state[i] - orig_state[i], 2);
        }
        std::cout << "\nState error: " << sqrt(err) << std::endl;

        
        FILE *ap_file = fopen("ap", "w");
        fwrite(ap.data(), sizeof(double), ap.size(), ap_file);
        fclose(ap_file);
        if (is_correct)
            std::cout << "Good" << std::endl;
        else
            std::cout << "Incorrect integration" << std::endl;
        delete [] constants;
        
    }*/
}


void test_function_example()
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0)
        std::cout << "Usage: ./ga config.json" << std::endl <<
        "Problem with a test function will be solved now" << std::endl;
        
    pcg_extras::seed_seq_from<std::random_device> seed_source;
        
    using FuncToOptimize = RosenbrockFunction<10>;
    FuncToOptimize func;
    FuncOptimization<FuncToOptimize, MinimizeFunc> optim(func);

    BasicPopulation pop(optim, 100, 10000);

    pop.init(pcg64(seed_source));

    genetic_algorithm(pop,
                    TournamentSelectionFast(pcg64(seed_source)),
                    SBXcrossover(pcg64(seed_source)),
                    PolynomialMutation<pcg64, pcg_extras::seed_seq_from<std::random_device>>(seed_source),
                    1000);
    if (rank == 0) {
        auto res = optim.get_result();
        std::cout << "Parameter error: " << func.solution_error(res) << std::endl;
    }
}


int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    if (argc == 2) {
        main_gen_algo(argv[1]);
    } else {
        test_function_example();
    }

    MPI_Finalize();
    return 0;
}
