// Main GA module
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
#include <omp.h>
#include "pcg_random.hpp"

#include "gradient_descent.h"
#include "genetic_algorithm.h"
#include "test_functions.h"
#include "polynomial_mutation.h"
#include "tournament_selection.h"
#include "sbx_crossover.h"
#include "cauchy_mutation.h"

#include "maleckar_model.h"
#include "kernik_clancy_model.h"

#include "cellml_ode_solver.h"
#include "optimization_problem.h"
#include "nelder_mead.h"
#include "cell_chain.h"

#include <json.hpp>
#include "mcmc.h"
#include "pso.h"

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
    popMal.init_selective(pcg64(seed_source), config["initial_selective_multiplier"].get<int>());
    time_population_init = MPI_Wtime() - time_population_init;
    
    if (mpi_rank == 0)
        std::cout << "time_population_init, s: " << time_population_init << std::endl;

    if (config["mutation_type"].get<std::string>() == "Cauchy") {
        genetic_algorithm(popMal,
            TournamentSelectionFast(pcg64(seed_source)),
            SBXcrossover(pcg64(seed_source), config["crossrate"].get<double>(), config["eta_crossover"].get<int>()),
            CauchyMutation<pcg64, SeedSource>
                                            (seed_source,
                                            config["mutrate"].get<double>(),
                                            config["gamma"].get<double>(),
                                            problem.get_gamma_vector()),
            config["n_generations"].get<unsigned>());
            
    } else if (config["mutation_type"].get<std::string>() == "Poly") {
        genetic_algorithm(popMal,
            TournamentSelectionFast(pcg64(seed_source)),
            SBXcrossover(pcg64(seed_source), config["crossrate"].get<double>(), config["eta_crossover"].get<int>()),
            PolynomialMutation<pcg64, SeedSource>
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
//const int seed_source = 42;
    std::ifstream configFile;
    configFile.open(configFilename);
    if (!configFile.is_open() && mpi_rank == 0) {
        std::cerr << "Cannot open main config file" << std::endl;
        throw 1;
    }
    json config;
    configFile >> config;
    configFile.close();


    //MaleckarModel model;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    KernikClancyModel model;

    ODESolver solver;
    //MOCKODESolver solver;
    MinimizeAPbaselines obj;
    ODEoptimization problem(model, solver, obj);

    double time_read_config = MPI_Wtime();
    try {
        problem.read_config(configFilename);
        problem.read_baselines(config);
    } catch(const char * err) {
        std::cout << "catch in main:" << std::endl;
        std::cout << err << std::endl;
        throw;
    }
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
    return;
        std::vector<std::pair<int, double>> error_per_gen;
        gen_algo_call(seed_source, problem, config, error_per_gen);
        if (mpi_rank == 0) {
            std::string filename = "convergence_GA.txt";
            std::ofstream file(filename);
            for (const auto & p: error_per_gen)
                file << p.first << " " << p.second << std::endl;
        }
        error_per_gen.clear();
        if (mpi_rank == 0) {
            //addition GD
            error_per_gen = weirdSteepestGradientDescent(pcg64(seed_source), problem, config["n_generations_GD"].get<int>(), 1e-1, problem.get_results_optimizer_format());
        }
        
     // error_per_gen = simpleGradientDescent(pcg64(seed_source), problem, config["n_generations"].get<int>());
     //  error_per_gen = weirdSteepestGradientDescent(pcg64(seed_source), problem, config["n_generations"].get<int>(), 1e-3);

        //write error_per_gen to file
        if (mpi_rank == 0) {
            std::string filename = "convergence_GD.txt";
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
        std::cout << std::endl << "Relative to CL = 1000 state" << std::endl;
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



void script_genetic_algorithm(json & config)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    //const int seed_source = 42;
    
    //MaleckarModel model;
    KernikClancyModel model;
    ODESolver solver;

    //MinimizeAPbaselines obj;
    ScaleMinimizeAPbaselines obj;
    
    //MaximizeAPinnerProduct obj;
    //LeastSquaresMinimizeAPbaselines obj;
    ODEoptimization problem(model, solver, obj);
    
    try {
        problem.read_config(config);
        problem.read_baselines(config);
    } catch (const char * err) {
        std::cout << "catch in main:" << std::endl;
        std::cout << err << std::endl;
        throw;
    }
    
    std::vector<std::pair<int, double>> error_per_gen;
    gen_algo_call(seed_source, problem, config, error_per_gen);
    if (mpi_rank == 0) {
        std::string filename = "convergence_GA.txt";
        std::ofstream file(filename);
        for (const auto & p: error_per_gen)
            file << p.first << " " << p.second << std::endl;
        problem.export_gen_algo_tables();
    }
    if (mpi_rank == 0) {
        problem.dump_ap(problem.get_results_optimizer_format(), 0);
        using Results = decltype(problem)::Results;
        using BaselineResult = decltype(problem)::BaselineResult;
        Results results = problem.get_relative_results();
        BaselineResult res = results[0];
        std::cout << "Printing relative to default results" << std::endl;
        for (const auto & cit: res.constantsResult)
            std::cout << cit.first << " " << cit.second << std::endl;
        std::cout << std::endl << "Relative to CL = 1000 state" << std::endl;
        for (const auto & sit: res.statesResult)
            std::cout << sit.first << " " << sit.second << std::endl;
    }
}






void script_PSO(json & config)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    pcg_extras::seed_seq_from<std::random_device> seed_source;
    //const int seed_source = 42;

    //MaleckarModel model;
    KernikClancyModel model;
    ODESolver solver;

    MinimizeAPbaselines obj;
    //ScaleMinimizeAPbaselines obj;

    //MaximizeAPinnerProduct obj;
    //LeastSquaresMinimizeAPbaselines obj;
    ODEoptimization problem(model, solver, obj);

    try {
        problem.read_config(config);
        problem.read_baselines(config);
    } catch (const char * err) {
        std::cout << "catch in main:" << std::endl;
        std::cout << err << std::endl;
        throw;
    }

    std::vector<std::pair<int, double>> error_per_gen;


    PSO_population<decltype(problem), pcg64, decltype(seed_source)> pop(problem,
                    config["n_organisms"].get<unsigned>(), seed_source);

    double time_population_init = MPI_Wtime();
    pop.init(pcg64(seed_source));
    time_population_init = MPI_Wtime() - time_population_init;

    if (mpi_rank == 0)
        std::cout << "time_population_init, s: " << time_population_init << std::endl;

    particle_swarm_optimization(pop, config["n_generations"].get<unsigned>());

    error_per_gen = pop.get_error_per_gen();

    if (mpi_rank == 0) {
        std::string filename = "convergence_GA.txt";
        std::ofstream file(filename);
        for (const auto & p: error_per_gen)
            file << p.first << " " << p.second << std::endl;
        //TODO
        //problem.export_gen_algo_tables();
    }
    if (mpi_rank == 0) {
        problem.dump_ap(problem.get_results_optimizer_format(), 0);
        using Results = decltype(problem)::Results;
        using BaselineResult = decltype(problem)::BaselineResult;
        Results results = problem.get_relative_results();
        BaselineResult res = results[0];
        std::cout << "Printing relative to default results" << std::endl;
        for (const auto & cit: res.constantsResult)
            std::cout << cit.first << " " << cit.second << std::endl;
        std::cout << std::endl << "Relative to CL = 1000 state" << std::endl;
        for (const auto & sit: res.statesResult)
            std::cout << sit.first << " " << sit.second << std::endl;
    }
}








void script_nelder_mead(json & config)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    pcg_extras::seed_seq_from<std::random_device> seed_source;
//const int seed_source = 42;


   // MaleckarModel model;
    KernikClancyModel model;

    ODESolver solver;
//MaximizeAPinnerProduct obj;
    MinimizeAPbaselines obj;
  // LeastSquaresMinimizeAPbaselines obj;
  //  ODEoptimizationTrackVersion problem(model, solver, obj);///////////////////////////
    ODEoptimization problem(model, solver, obj);

    double time_read_config = MPI_Wtime();
    try {
        problem.read_config(config);
        problem.read_baselines(config);
    } catch(const char * err) {
        std::cout << "catch in main:" << std::endl;
        std::cout << err << std::endl;
        throw;
    }
    time_read_config = MPI_Wtime() - time_read_config;

    if (mpi_rank == 0) {
        std::cout << "time_read_config, s: " << time_read_config << std::endl;

        std::vector<std::pair<int, double>> error_per_gen = nelder_mead(problem, config["NM_limit_calls"].get<int>(), 1e-14, 1, config["NM_simplex_step"].get<double>());
        auto res1 = problem.get_results_optimizer_format();

        problem.dump_ap(res1.begin(), 5, 1000);


       // problem.unfreeze_global_variable("i_stim_Amplitude", 5, 100, res1);
       // problem.beats = 100;
       // std::cout << "stage 2" << std::endl;
       // std::vector<std::pair<int, double>> error_per_gen2 = nelder_mead(problem, config["NM_limit_calls"].get<int>(), 1e-14, 1, config["NM_simplex_step"].get<double>()/10, res1);
       // error_per_gen.insert(error_per_gen.end(), error_per_gen2.begin(), error_per_gen2.end());
       // auto res2 = problem.get_results_optimizer_format();
       // problem.dump_ap(res2.begin(), 10);
        
        std::string filename = "convergence_NM.txt";
        std::ofstream file(filename);
        for (const auto & p: error_per_gen)
            file << p.first << " " << p.second << std::endl;

        using Results = decltype(problem)::Results;
        using BaselineResult = decltype(problem)::BaselineResult;
        Results results = problem.get_relative_results();
        BaselineResult res = results[0];
        std::cout << "Printing relative to default results" << std::endl;
        for (const auto & cit: res.constantsResult)
            std::cout << cit.first << " " << cit.second << std::endl;
        std::cout << std::endl << "Relative to CL = 1000 state" << std::endl;
        for (const auto & sit: res.statesResult)
            std::cout << sit.first << " " << sit.second << std::endl;
    }
}

void script_direct_problem(json & config)
{
  //  MaleckarModel model;
    KernikClancyModel model;

    ODESolver solver;
    LeastSquaresMinimizeAPbaselines obj;
    ODEoptimization problem(model, solver, obj);

    try {
        problem.read_config(config);
    } catch(const std::string & err) {
        std::cout << "catch in main:" << std::endl;
        std::cout << err << std::endl;
        throw;
    }

    problem.run_direct_and_dump(config["start_record_time"].get<double>(),
                                config["max_time"].get<double>(),
                                config["dump_period"].get<double>(),
                                config["dump_filename"].get<std::string>(),
                                config["dump_vars"].get<std::vector<std::string>>());

}

void script_direct_problem_chain(json & config)
{
	CellChainModel<KernikClancyModel> model(config["chain_len"].get<int>(), config["main_cell_index"].get<int>());

	ODESolver solver;
    LeastSquaresMinimizeAPbaselines obj;
    ODEoptimization problem(model, solver, obj);

    try {
        problem.read_config(config);
    } catch(const std::string & err) {
        std::cout << "catch in main:" << std::endl;
        std::cout << err << std::endl;
        throw;
    }

    problem.run_direct_and_dump(config["start_record_time"].get<double>(),
                                config["max_time"].get<double>(),
                                config["dump_period"].get<double>(),
                                config["dump_filename"].get<std::string>(),
                                config["dump_vars"].get<std::vector<std::string>>());
}


void script_gradient_descent(json & config)
{
    
}

void script_test_function(json & config)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
        std::cout << "Problem with a test function will be solved now" << std::endl;

    pcg_extras::seed_seq_from<std::random_device> seed_source;
   // int seed_source = 42;
    using FuncToOptimize = RosenbrockFunction<8>;
    FuncToOptimize func;
    FuncOptimization<FuncToOptimize, MinimizeFunc> optim(func);
/*
    {
        PSO_population<decltype(optim), pcg64, decltype(seed_source)> pop(optim,
                    config["n_organisms"].get<unsigned>(), seed_source);

        double time_population_init = MPI_Wtime();
        pop.init(pcg64(seed_source));
        time_population_init = MPI_Wtime() - time_population_init;

        if (rank == 0)
            std::cout << "time_population_init, s: " << time_population_init << std::endl;

        particle_swarm_optimization(pop, config["n_generations"].get<unsigned>());
    }
*/

  //  auto res = nelder_mead(pcg64(seed_source), optim, 10000000, 1e-14, 100, 1e-1);
  //  for (auto & p: res) {
    //    std::cout << p.first << " " << p.second << std::endl;
   // }
//    BasicPopulation pop(optim, 10, 100);

  //  pop.init(pcg64(seed_source));

    //genetic_algorithm(pop,
      //              TournamentSelectionFast(pcg64(seed_source)),
        //            SBXcrossover(pcg64(seed_source), 0.1, 10),
                 /*
                    CauchyMutation<pcg64, decltype(seed_source)>
                                            (seed_source,
                                            0.1,
                                            0.01,
                                            optim.get_gamma_vector()),
                   */
                    //NoMutation(),
          //         PolynomialMutation<pcg64, decltype(seed_source)>(seed_source, 0.1, 200),
            //        100);


 //   simpleGradientDescent(pcg64(seed_source), optim, 1000, 1e-6);
//weirdSteepestGradientDescent(pcg64(seed_source), optim, 3000, 1e-3);
    if (rank == 0) {
        auto res = optim.get_result();
        std::cout << "Parameter error: " << func.solution_error(res) << std::endl;
    }
}

void script_track_minimum(json & config)
{
    
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

//    pcg_extras::seed_seq_from<std::random_device> seed_source;
const int seed_source = 42;


   // MaleckarModel model;
    KernikClancyModel model;

    ODESolver solver;
//MaximizeAPinnerProduct obj;
    MinimizeAPbaselines obj;
  // LeastSquaresMinimizeAPbaselines obj;
  //  ODEoptimizationTrackVersion problem(model, solver, obj);///////////////////////////
    ODEoptimization problem(model, solver, obj);

    double time_read_config = MPI_Wtime();
    try {
        problem.read_config(config);
        problem.read_baselines(config);
    } catch(const char * err) {
        std::cout << "catch in main:" << std::endl;
        std::cout << err << std::endl;
        throw;
    }
    time_read_config = MPI_Wtime() - time_read_config;
return;
    if (mpi_rank == 0) {
        std::cout << "time_read_config, s: " << time_read_config << std::endl;

        const int global_steps = 10;
        const int local_steps = 20;
        const double r_eps = 1e-1;
        //std::vector<std::pair<int, double>> error_per_gen = weirdSteepestGradientDescentTrack(pcg64(seed_source), problem, r_eps, global_steps, local_steps);////////
     //   std::vector<std::pair<int, double>> error_per_gen = weirdSteepestGradientDescent(pcg64(seed_source), problem, 15, r_eps);

     //   error_per_gen = weirdSteepestGradientDescent(pcg64(seed_source), problem, 10, 1e-7, problem.get_results_optimizer_format());
       // error_per_gen = weirdSteepestGradientDescent(pcg64(seed_source), problem, 10, 1e-9, problem.get_results_optimizer_format());
        //error_per_gen = weirdSteepestGradientDescent(pcg64(seed_source), problem, 10, 1e-3, problem.get_results_optimizer_format());
       // error_per_gen = weirdSteepestGradientDescent(pcg64(seed_source), problem, 10, 1e-4, problem.get_results_optimizer_format());

        std::vector<std::pair<int, double>> error_per_gen = nelder_mead(problem, 500, 1e-14, 1, 1e-1);
        auto res1 = problem.get_results_optimizer_format();
        
        
        problem.dump_ap(res1.begin(), 5);
        /*
        problem.unfreeze_global_variable("i_stim_Amplitude", 5, 100, res1);
        std::vector<std::pair<int, double>> error_per_gen2 = nelder_mead(problem, 500, 1e-14, 1, 1e-1, res1);
        error_per_gen.insert(error_per_gen.end(), error_per_gen2.begin(), error_per_gen2.end());
        auto res2 = problem.get_results_optimizer_format();
        problem.dump_ap(res2.begin(), 10);
        */
        //std::string filename = "convergence_TRACK.txt";
        std::string filename = "convergence_CG.txt";
        std::ofstream file(filename);
        for (const auto & p: error_per_gen)
            file << p.first << " " << p.second << std::endl;

        using Results = decltype(problem)::Results;
        using BaselineResult = decltype(problem)::BaselineResult;
        Results results = problem.get_relative_results();
        BaselineResult res = results[0];
        std::cout << "Printing relative to default results" << std::endl;
        for (const auto & cit: res.constantsResult)
            std::cout << cit.first << " " << cit.second << std::endl;
        std::cout << std::endl << "Relative to CL = 1000 state" << std::endl;
        for (const auto & sit: res.statesResult)
            std::cout << sit.first << " " << sit.second << std::endl;
    }
}


void script_mcmc(json & config)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


#ifdef MCMC_ENABLED
    pcg_extras::seed_seq_from<std::random_device> seed_source;
//const int seed_source = 42;


   // MaleckarModel model;
    KernikClancyModel model;

    ODESolver solver;
//MaximizeAPinnerProduct obj;
    MinimizeAPbaselines obj;
  // LeastSquaresMinimizeAPbaselines obj;
  //  ODEoptimizationTrackVersion problem(model, solver, obj);///////////////////////////
    ODEoptimization problem(model, solver, obj);

    double time_read_config = MPI_Wtime();
    try {
        problem.read_config(config);
        problem.read_baselines(config);
    } catch(const char * err) {
        std::cout << "catch in main:" << std::endl;
        std::cout << err << std::endl;
        throw;
    }
    time_read_config = MPI_Wtime() - time_read_config;


    std::vector<double> res(problem.get_number_parameters());
    problem.initial_guess_for_optimizer(res.begin());
    mcmc(problem, config["output_name"].get<std::string>(), res);

return;
    if (mpi_rank == 0) {
        std::cout << "time_read_config, s: " << time_read_config << std::endl;

        std::vector<std::pair<int, double>> error_per_gen = nelder_mead(problem, config["NM_limit_calls"].get<int>(), 1e-14, 1, config["NM_simplex_step"].get<double>());
        auto res1 = problem.get_results_optimizer_format();

        problem.dump_ap(res1.begin(), 5);
        /*
        problem.unfreeze_global_variable("i_stim_Amplitude", 5, 100, res1);
        std::vector<std::pair<int, double>> error_per_gen2 = nelder_mead(problem, 500, 1e-14, 1, 1e-1, res1);
        error_per_gen.insert(error_per_gen.end(), error_per_gen2.begin(), error_per_gen2.end());
        auto res2 = problem.get_results_optimizer_format();
        problem.dump_ap(res2.begin(), 10);
        */
        std::string filename = "convergence_NM.txt";
        std::ofstream file(filename);
        for (const auto & p: error_per_gen)
            file << p.first << " " << p.second << std::endl;

        using Results = decltype(problem)::Results;
        using BaselineResult = decltype(problem)::BaselineResult;
        Results results = problem.get_relative_results();
        BaselineResult res = results[0];
        std::cout << "Printing relative to default results" << std::endl;
        for (const auto & cit: res.constantsResult)
            std::cout << cit.first << " " << cit.second << std::endl;
        std::cout << std::endl << "Relative to CL = 1000 state" << std::endl;
        for (const auto & sit: res.statesResult)
            std::cout << sit.first << " " << sit.second << std::endl;



        mcmc(problem, config["output_name"].get<std::string>(), res1);

    }
#else
    if (mpi_rank == 0)
        std::cerr << "Build with no MCMC support" << std::endl;
#endif

}


/**
 * @brief Entry point
 *
 * Execution of the program starts here.
 * The program expects configuration file to run.
 * Please refer to example config files.
 *
 * @param argv[1] Path to config file
 * 
 */
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (argc != 2) {
        if (mpi_rank == 0)
            std::cout << "Usage: ./ga config.json" << std::endl;
        MPI_Finalize();
        return 0;
    }

    std::ifstream configFile(argv[1]);
    if (!configFile.is_open()) {
        std::cerr << "Node " << mpi_rank << ": Cannot open config file " << argv[1] << std::endl;
        MPI_Finalize();
        return 0;
    }
    json config;
    configFile >> config;

    if (mpi_rank == 0) {
        printf("Number of MPI nodes: %d\n", mpi_size);
        printf("Number of OpenMP threads at root: %d\n", omp_get_max_threads());
    }

    const std::string sname = config["script"].get<std::string>();
    try {
        if (sname == "Genetic Algorithm") {
            script_genetic_algorithm(config);
        } else if (sname == "PSO") {
            script_PSO(config);
        }else if (sname == "Nelder Mead") {
            script_nelder_mead(config);
        } else if (sname == "Direct Problem") {
            script_direct_problem(config);
        } else if (sname == "Direct Problem Chain") {
			script_direct_problem_chain(config);
        } else if (sname == "Gradient Descent") {
            script_gradient_descent(config);
        } else if (sname == "Test Function") {
            script_test_function(config);
        } else if (sname == "Track Minimum") {
            script_track_minimum(config);
        } else if (sname == "MCMC") {
            script_mcmc(config);
        } else {
            std::cout << "Unknown script name: " << sname << std::endl;
        }
    } catch(const std::string & e) {
        std::cerr << "Node " << mpi_rank << ": Exception string caught in main:\n" << e << std::endl;
    }
    MPI_Finalize();
    return 0;
}
