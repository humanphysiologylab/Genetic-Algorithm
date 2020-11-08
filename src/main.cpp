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
#include "fitness_functors.h"
#include "polynomial_mutation.h"
#include "tournament_selection.h"
#include "sbx_crossover.h"
#include "cauchy_mutation.h"

#include "maleckar_model.h"
#include "cellml_ode_solver.h"


#ifdef HIDE_CODE

class FitnessFunctor
{
    std::vector<double> AP_control;
    std::vector<int> baseline_periods;
public:
    std::vector<int> get_baseline_periods() const
    {
        return baseline_periods;
    }
    int get_sum_of_periods() const
    {
        return AP_control.size();
    }
    double operator()(const double * AP) const
    {
        double res = 0;
        for (int i = 0; i < AP_control.size(); i++)
            res += std::pow(AP[i] - AP_control[i], 2);
        return res;
    }
    FitnessFunctor(std::vector<std::string> & baseline_file_names)
    {
        for (auto & filename: baseline_file_names) {
            FILE *file = fopen(filename.c_str(), "r");
            if (file == NULL) {
                std::cerr << "Cannot open baseline file: " << filename << std::endl;
                throw;
            }
            
            int period;
            double buf;
            for (period = 0;; period++) {
                if (fscanf(file, "%lf\n", &buf) == EOF)
                    break;
                
                AP_control.push_back(buf);
            }
            assert(period > 0);
            baseline_periods.push_back(period);
            fclose(file);
        }
    }
};

class Population
{
    
public:
    FitnessFunctor fitnessFunctor;

    int number_organisms;
    int number_genes;
    int number_baselines;
    int number_elites;
    int INIT_FROM_BACKUP_FILE;
    int period_backup;
    int number_mutants;
    
    int rank;
    int size;
    
    MPI_Datatype GlobalSetupMPI;
    MPI_Datatype StateVectorMPI;
    
    std::vector<double> left_border, right_border;
    
    Population()
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    
        /* create an MPI type for struct GlobalSetup */
        {
            int blocklengths[GlobalSetupItemsNumber];
            std::fill_n(blocklengths, GlobalSetupItemsNumber, 1);
            MPI_Datatype types[GlobalSetupItemsNumber];
            std::fill_n(types, GlobalSetupItemsNumber, MPI_INT);

            MPI_Aint displacements[GlobalSetupItemsNumber];

            displacements[0] = offsetof(GlobalSetup, number_organisms);
            displacements[1] = offsetof(GlobalSetup, number_genes);
            displacements[2] = offsetof(GlobalSetup, number_generations);
            displacements[3] = offsetof(GlobalSetup, number_baselines);
            displacements[4] = offsetof(GlobalSetup, number_elites);
            displacements[5] = offsetof(GlobalSetup, INIT_FROM_BACKUP_FILE);
            displacements[6] = offsetof(GlobalSetup, period_backup);
            displacements[7] = offsetof(GlobalSetup, number_mutants);

            MPI_Type_create_struct(GlobalSetupItemsNumber, blocklengths, displacements, types, &GlobalSetupMPI);
            MPI_Type_commit(&GlobalSetupMPI);
        }
        /* done */

        /* Create MPI structure for State */
        {
            int blocklengths[STATE_ARRAY_SIZE];
            std::fill_n(blocklengths, STATE_ARRAY_SIZE, 1);
            MPI_Datatype types[STATE_ARRAY_SIZE];
            std::fill_n(types, STATE_ARRAY_SIZE, MPI_DOUBLE);

            MPI_Aint displacements[STATE_ARRAY_SIZE];

            displacements[0] = offsetof(State, V);
            displacements[1] = offsetof(State, Na_c);
            displacements[2] = offsetof(State, Na_i);
            displacements[3] = offsetof(State, m);
            displacements[4] = offsetof(State, h1);
            displacements[5] = offsetof(State, h2);
            displacements[6] = offsetof(State, Ca_d);
            displacements[7] = offsetof(State, d_L);
            displacements[8] = offsetof(State, f_L1);
            displacements[9] = offsetof(State, f_L2);
            displacements[10] = offsetof(State, K_c);
            displacements[11] = offsetof(State, K_i);
            displacements[12] = offsetof(State, r);
            displacements[13] = offsetof(State, s);
            displacements[14] = offsetof(State, a_ur);
            displacements[15] = offsetof(State, i_ur);
            displacements[16] = offsetof(State, n);
            displacements[17] = offsetof(State, pa);
            displacements[18] = offsetof(State, Ca_c);
            displacements[19] = offsetof(State, Ca_i);
            displacements[20] = offsetof(State, O_C);
            displacements[21] = offsetof(State, O_TC);
            displacements[22] = offsetof(State, O_TMgC);
            displacements[23] = offsetof(State, O_TMgMg);
            displacements[24] = offsetof(State, O);
            displacements[25] = offsetof(State, Ca_rel);
            displacements[26] = offsetof(State, Ca_up);
            displacements[27] = offsetof(State, O_Calse);
            displacements[28] = offsetof(State, F1);
            displacements[29] = offsetof(State, F2);
            displacements[30] = offsetof(State, d_ord);
            displacements[31] = offsetof(State, ff);
            displacements[32] = offsetof(State, fs);
            displacements[33] = offsetof(State, fcaf);
            displacements[34] = offsetof(State, fcas);
            displacements[35] = offsetof(State, jca);
            displacements[36] = offsetof(State, ffp);
            displacements[37] = offsetof(State, fcafp);
            displacements[38] = offsetof(State, nca);

            MPI_Type_create_struct(STATE_ARRAY_SIZE, blocklengths, displacements, types, &StateVectorMPI);
            MPI_Type_commit(&StateVectorMPI);
        }
        //MPI_Type_contiguous(STATE_ARRAY_SIZE, MPI_DOUBLE, &StateVectorMPI);
        //MPI_Type_commit(&StateVectorMPI);
        /* done */
            
    }
    
    ~Population()
    {
        MPI_Type_free(&GlobalSetupMPI);
        MPI_Type_free(&StateVectorMPI);
    }
    
    
    void scatter()
    {
        MPI_Scatter(next_generation, gs.number_organisms * gs.number_genes / size,
                    MPI_DOUBLE, (rank == 0) ? MPI_IN_PLACE : next_generation,
                    gs.number_organisms * gs.number_genes / size,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Scatter(state_struct, gs.number_baselines * gs.number_organisms / size,
                    StateVectorMPI, (rank == 0) ? MPI_IN_PLACE : state_struct,
                    gs.number_baselines * gs.number_organisms / size,
                    StateVectorMPI, 0, MPI_COMM_WORLD);   
    }
    void gather()
    {
        MPI_Gather((rank == 0) ? MPI_IN_PLACE : AP_current, (time_sum) * gs.number_organisms / size, MPI_DOUBLE,
                   AP_current, (time_sum) * gs.number_organisms / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gather((rank == 0) ? MPI_IN_PLACE : state_struct, gs.number_baselines * gs.number_organisms / size,
                   StateVectorMPI,
                   state_struct, gs.number_baselines * gs.number_organisms / size, StateVectorMPI, 0, MPI_COMM_WORLD);

        MPI_Gather((rank == 0) ? MPI_IN_PLACE : next_generation, gs.number_genes * gs.number_organisms / size,
                   MPI_DOUBLE, next_generation,
                   gs.number_genes * gs.number_organisms / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    void run_generation()
    {
    
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
    
    }


    void fitness_function(std::vector<std::pair<double, int>> & sd_n_index)
    {
        fitness_function(AP_control, AP_current, best_scaling_factor, best_scaling_shift, TIME, sd_n_index,
                             gs.number_organisms, gs.number_baselines, time_sum);
    }
    
    
    void save_elite_to_elite_buffer(int from, int to)
    {
        for (int j = 0; j < gs.number_baselines; j++)
            buf_elite_state[to * gs.number_baselines + j] =
                    state_struct[from * gs.number_baselines + j];

        for (int j = 0; j < gs.number_genes; j++)
            buf_elite_genes[to * gs.number_genes + j] =
                next_generation[from * gs.number_genes + j];
    }
    
    void save_mutant_to_mutant_buffer(int from, int to)
    {
        for (int j = 0; j < gs.number_baselines; j++)
            buf_mutant_state[to * gs.number_baselines + j]
                = state_struct[from * gs.number_baselines + j];

        for (int j = 0; j < gs.number_genes; j++)
            buf_mutant_genes[to * gs.number_genes + j] =
                next_generation[from * gs.number_genes + j];
    }
    
    void restore_elites_to_main_array()
    {
            for (int i = 0; i < gs.number_elites; i++) {
                for (int j = 0; j < gs.number_baselines; j++)
                    state_struct[i * gs.number_baselines + j] =
                            buf_elite_state[i * gs.number_baselines + j];

                for (int j = 0; j < gs.number_genes; j++)
                    next_generation[i * gs.number_genes + j] =
                            buf_elite_genes[i * gs.number_genes + j];
            }
    }
    void restore_mutants_to_main_array()
    {
            for (int i = 0; i < gs.number_mutants; i++) {
                const int index = i + gs.number_elites;
                for (int j = 0; j < gs.number_baselines; j++)
                    state_struct[index * gs.number_baselines + j] =
                            buf_mutant_state[i * gs.number_baselines + j];

                for (int j = 0; j < gs.number_genes; j++)
                    next_generation[index * gs.number_genes + j] =
                            buf_mutant_genes[i * gs.number_genes + j];
            }
    }


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



static char *read_line(char *pcBuf, int iMaxSize, FILE *fStream)
{
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
            std::cerr << "Error! Unable to read line!" << std::endl;
            exit(-1);
        }

        if (feof(fStream)) break;
        p = strtok(pcBuf, "\n\r");
    } while ((p == NULL) || (*p == '#'));

    return p;
}

struct GlobalSetup
{
    int number_organisms;
    int number_genes;
    int number_generations;
    int number_baselines;
    int number_elites;
    int INIT_FROM_BACKUP_FILE;
    int period_backup;
    int number_mutants;
};

const int GlobalSetupItemsNumber = 8;

void print_log_stdout(struct GlobalSetup gs, const double *genes,
                      const std::vector<std::pair<double, int>> &sd_n_index)
{

    double average_error = 0;
    for (int i = 0; i < gs.number_organisms; i++)
        average_error += sd_n_index[i].first;
    average_error /= gs.number_organisms;

    std::cout << "Best SD: " << sd_n_index[0].first << " mV\n";
    std::cout << "Worst SD: " << sd_n_index.back().first << " mV\n";
    std::cout << "Average SD: " << average_error << " mV\n";
    std::cout << "The fittest organism:\n";

    const int index_best = sd_n_index[0].second;

    const int number_conc_types = 3;
    std::string conc[] = {"Na_i", "Ca_rel", "K_i"};
    const int number_cond = gs.number_genes - number_conc_types * gs.number_baselines;
    std::string mult[] = {"INa", "ICaL", "Ito", "IKur", "IK1",
                          "IKs", "IKr", "IBNa", "IBCa", "INaK",
                          "ICaP", "INaCa", "Iup", "Irel"};
    const int stride = 4;

    for (int i = 0; i < number_cond; ++i) {
        std::cout << mult[i] << " " << genes[index_best * gs.number_genes + i] << "\t";
        if (i % stride == stride - 1) std::cout << std::endl;
    }
    std::cout << std::endl;

    for (int i = 0; i < number_conc_types; ++i) {
        std::cout << conc[i] << ": ";
        for (int j = 0; j < gs.number_baselines; ++j) {
            std::cout << genes[index_best * gs.number_genes + number_cond + i * gs.number_baselines + j] << " ";
        }
        std::cout << std::endl;
    }
}


void old_code(int argc, char *argv[])
{
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);



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

int main(int argc, char *argv[])
{
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    MPI_Init(&argc, &argv);

    pcg_extras::seed_seq_from<std::random_device> seed_source;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    RosenbrockFunction func(15);

    BasicPopulation pop(func, ScalarFunctionMinFitnessFunctor(), 100, 10000);

    pop.init(pcg64(seed_source));

    genetic_algorithm(pop,
                    TournamentSelectionFast(pcg64(seed_source)),
                    SBXcrossover(pcg64(seed_source)),
                    PolynomialMutation<pcg64, pcg_extras::seed_seq_from<std::random_device>>(seed_source),
                    100);

    if (rank == 0) {
        std::cout << "Error: " << func.solution_error_l2(pop.best()) << std::endl;
        
        
        //try malecar on one node
        MaleckarModel model;
        std::vector<double> state(model.state_size());
        double * constants = new double [model.constants_size()];
        model.initConsts(constants);
        model.set_constants(constants);
        model.initState(state.data());
        std::vector<double> orig_state = state;
        int is_correct;
        double t0 = 0, start_record = 999, tout = start_record + 1;
        std::vector<double> ap(1000);
        
        double st = MPI_Wtime();
        solve(model, state, is_correct, t0, start_record, tout, ap);
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
    }
    MPI_Finalize();
    return 0;
}
