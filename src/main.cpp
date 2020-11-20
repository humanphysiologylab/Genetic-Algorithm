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


template <typename SeedSource, typename Problem>
void gen_algo_call(SeedSource & seed_source, Problem & problem, json & config, std::vector<std::pair<int, double>> & error_per_gen)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    BasicPopulation popMal(problem,
                    config["n_elites"].get<unsigned>(),
                    config["n_organisms"].get<unsigned>());

    double time_population_init = MPI_Wtime();
    popMal.init(pcg64(seed_source));
    time_population_init = MPI_Wtime() - time_population_init;
    
    if (mpi_rank == 0)
        std::cout << "time_population_init, s: " << time_population_init << std::endl;
    
    if (config["mutation_type"].get<std::string>() == "Cauchy") {
        genetic_algorithm(popMal,
            TournamentSelectionFast(pcg64(seed_source)),
            SBXcrossover(pcg64(seed_source), config["crossrate"].get<double>(), config["eta_crossover"].get<int>()),
            CauchyMutation<pcg64, pcg_extras::seed_seq_from<std::random_device>>
                                            (seed_source,
                                            config["mutrate"].get<double>(),
                                            config["gamma"].get<double>(),
                                            problem.get_gamma_vector()),
            config["n_generations"].get<unsigned>());
            
    } else if (config["mutation_type"].get<std::string>() == "Poly") {
        genetic_algorithm(popMal,
            TournamentSelectionFast(pcg64(seed_source)),
            SBXcrossover(pcg64(seed_source), config["crossrate"].get<double>(), config["eta_crossover"].get<int>()),
            PolynomialMutation<pcg64, pcg_extras::seed_seq_from<std::random_device>>
                                            (seed_source,
                                            config["mutrate"].get<double>(),
                                            config["eta_mutation"].get<int>()),
            config["n_generations"].get<unsigned>());
    } else if (config["mutation_type"].get<std::string>() == "None") {
        genetic_algorithm(popMal,
            TournamentSelectionFast(pcg64(seed_source)),
            SBXcrossover(pcg64(seed_source), config["crossrate"].get<double>(), config["eta_crossover"].get<int>()),
            NoMutation(),
            config["n_generations"].get<unsigned>());
    } else {
        throw("Unknown mutation in config");
    }

    error_per_gen = popMal.get_error_per_gen();
}

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

    if (mpi_rank == 0) {
        std::cout << "time_read_config, s: " << time_read_config << std::endl;
    }
/* test Polynomial Mutation
    const double crossrates[] = {0.1, 0.5, 0.9};
    const int eta_crossovers[] = {10, 20, 200};
    const double mutrates[] = {0.1, 0.5, 0.9};
    const int eta_mutations[] = {10, 20, 200};
    for (double crossrate: crossrates) {
        for (int eta_crossover: eta_crossovers) {
            for (double mutrate: mutrates) {
                for (int eta_mutation: eta_mutations) {
    config["crossrate"] = crossrate;
    config["eta_crossover"] = eta_crossover;
    config["mutrate"] = mutrate;
    config["eta_mutation"] = eta_mutation;
     
    const int num_tries = 5;
    for (int i = 0; i < num_tries; i++) {
        std::vector<std::pair<int, double>> error_per_gen;
        gen_algo_call(seed_source, problem, config, error_per_gen);
        //write error_per_gen to file
        if (mpi_rank == 0) {
            std::ostringstream crossrate_string, mutrate_string;
            crossrate_string << std::noshowpoint << crossrate;
            mutrate_string << std::noshowpoint << mutrate;
            
            std::string filename = crossrate_string.str() + 
                                   "_" + std::to_string(eta_crossover) + "_" +
                                   mutrate_string.str() + "_" + 
                                   std::to_string(eta_mutation) + "_" + std::to_string(i);
            
            std::ofstream file(filename);
            for (const auto & p: error_per_gen)
                file << p.first << " " << p.second << std::endl;
        }
    }
                }
            }
        }
    }
*/

        std::vector<std::pair<int, double>> error_per_gen;
        gen_algo_call(seed_source, problem, config, error_per_gen);
        //write error_per_gen to file
        if (mpi_rank == 0) {
            std::string filename = "convergence_log.txt";
            std::ofstream file(filename);
            for (const auto & p: error_per_gen)
                file << p.first << " " << p.second << std::endl;
        }


    
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
        std::cout << "relative state" << std::endl;
        for (const auto & sit: res.statesResult)
            std::cout << sit.first << " " << sit.second << std::endl;
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
        double t0 = 0, start_record = 999, tout = start_record + 1;
        std::vector<double> ap(1000);
        
        double st = MPI_Wtime();
        solver.solve(model, state, is_correct, t0, start_record, tout, ap);
        st = MPI_Wtime() - st;
        std::cout << "TIME: " << st << std::endl;
        st = MPI_Wtime();
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
        
    }
    */
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
