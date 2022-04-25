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
#include <iomanip>
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
#include "optimization_problem_rest.h"
#include "objective.h"
#include "nelder_mead.h"
#include "cell_chain.h"

#include <json.hpp>
#include "mcmc.h"
#include "pso.h"

/**
 * @brief
 *
 * The function generates unique_ptr to mutation object for corresponding @p name
 *
 * @param[in] seed_source
 * @param[in] name possible names: Poly, Cauchy, None
 * @param[in] config
 * @param[in] problem
 *
 * @return unique_ptr to mutation object
 *
 * @throw std::logic_error If mutation type is unknown
 */
template <typename SeedSource, typename Problem>
std::unique_ptr<BaseMutation> new_mutation(SeedSource & seed_source, const std::string & name, json & config, const Problem & problem)
{
    if (name == "Poly")
        return std::make_unique<PolynomialMutation<pcg64, SeedSource>>
                                                (seed_source,
                                                config["mutrate"].get<double>(),
                                                config["eta_mutation"].get<int>());
    if (name == "Cauchy")
       return std::make_unique<CauchyMutation<pcg64, SeedSource>>
                                                (seed_source,
                                                config["mutrate"].get<double>(),
                                                config["gamma"].get<double>(),
                                                problem.get_gamma_vector());
    if (name == "None")
        return std::make_unique<NoMutation>();
    throw std::logic_error("Unknown mutation type");
}

/**
 * @brief Genetic algorithm for abstract optimization problem
 *
 * @param[in] seed_source Seed source for pcg64 pseudorandom number generator
 * @param[in,out] problem Initialized problem to optimize with genetic algorithm
 * @param[in] config Json config with genetic algorithm parameters
 * @param[out] error_per_gen Vector storing error per generation
 *
 */
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

    auto mutation = new_mutation(seed_source, config["mutation_type"].get<std::string>(), config, problem);
    auto selection = TournamentSelectionFast(pcg64(seed_source));
    auto crossover = SBXcrossover(pcg64(seed_source), config["crossrate"].get<double>(), config["eta_crossover"].get<int>());
    genetic_algorithm(popMal,
            selection,
            crossover,
            *mutation,
            config["n_generations"].get<unsigned>());
    error_per_gen = popMal.get_error_per_gen();
}

template <typename SeedSource, typename Problem>
void PSO_call(SeedSource & seed_source, Problem & problem, json & config, std::vector<std::pair<int, double>> & error_per_gen)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    PSO_population<Problem, pcg64, SeedSource> pop(problem,
                    config["n_organisms"].get<unsigned>(), seed_source);

    double time_population_init = MPI_Wtime();
    pop.init(pcg64(seed_source));
    time_population_init = MPI_Wtime() - time_population_init;

    if (mpi_rank == 0)
        std::cout << "time_population_init, s: " << time_population_init << std::endl;

    particle_swarm_optimization(pop, config["n_generations"].get<unsigned>());

    error_per_gen = pop.get_error_per_gen();
}


#if 0
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
#endif


template<typename Problem>
void nelder_mead_call(Problem & problem, json & config, std::vector<std::pair<int, double>> & error_per_gen,
        std::vector<double> init_vector = std::vector<double>())
{
    if (init_vector.size() > 0)
        error_per_gen = nelder_mead(problem, config["NM_limit_calls"].get<int>(), 1e-14, 1, config["NM_simplex_step"].get<double>(), init_vector);
    else
        error_per_gen = nelder_mead(problem, config["NM_limit_calls"].get<int>(), 1e-14, 1, config["NM_simplex_step"].get<double>());

    // problem.unfreeze_global_variable("i_stim_Amplitude", 5, 100, res1);
   // problem.beats = 100;
   // std::cout << "stage 2" << std::endl;
   // std::vector<std::pair<int, double>> error_per_gen2 = nelder_mead(problem, config["NM_limit_calls"].get<int>(), 1e-14, 1, config["NM_simplex_step"].get<double>()/10, res1);
   // error_per_gen.insert(error_per_gen.end(), error_per_gen2.begin(), error_per_gen2.end());
   // auto res2 = problem.get_results_optimizer_format();
   // problem.dump_ap(res2.begin(), 10);
}

template<typename Problem>
void dump_table_ode_problem(Problem & problem, const std::vector<double> & res, const std::string & filename)
{
    int param_num = problem.get_number_parameters();
    if (res.size() % (param_num + 1) != 0) {
        throw("dump_table: incorrect res vector size");
    }
    std::ofstream file(filename);
    if (!file.is_open())
        throw("dump_table: cannot create file");
    //create header
    std::vector<std::string> header_columns = problem.get_param_names();
    for (int j = 0; j < param_num; j++)
        file << header_columns[j] << " ";
    file << "loss" << std::endl;
    file << std::scientific << std::setprecision(12);

    std::vector<double> model_params(param_num);
    for (int i = 0; i < res.size() / (param_num + 1); i++) {
        problem.optimizer_model_scale(res.begin() + i * (param_num + 1), model_params.begin());
        for (int j = 0; j < param_num; j++)
            file << model_params[j] << " ";
        file << res[i * (param_num + 1) + param_num] << std::endl;
    }
}

template<typename Problem>
void dump_table(Problem & problem, const std::vector<double> & res, const std::string & filename)
{
    int param_num = problem.get_number_parameters();
    if (res.size() % (param_num + 1) != 0) {
        throw("dump_table: incorrect res vector size");
    }
    std::ofstream file(filename);
    if (!file.is_open())
        throw("dump_table: cannot create file");
    //create header
    for (int j = 0; j < param_num; j++)
        file << "p" << j << " ";
    file << "loss" << std::endl;
    file << std::scientific << std::setprecision(12);

    for (size_t i = 0; i < res.size() / (param_num + 1); i++) {
        for (int j = 0; j < param_num; j++)
            file << res[i * (param_num + 1) + j] << " ";
        file << res[i * (param_num + 1) + param_num] << std::endl;
    }
}

template<typename Problem>
void gd_call(Problem & problem, json & config, std::vector<std::pair<int, double>> & error_per_gen,
        std::vector<double> init_vector = std::vector<double>())
{
    if (config["MultipleStarts"].get<bool>()) {
        std::vector<double> res = multipleStartsScript(config, problem);
        dump_table(problem, res, config["MS_output_filename"].get<std::string>());
    } else {
        const int gd_max_step = config["GD_max_step"].get<int>();
        if (config["GD_type"].get<std::string>() == "Simple") {
            error_per_gen = simpleGradientDescent(problem, gd_max_step,
                config["GD_r_eps"].get<double>(), config["GD_learning_rate"].get<double>(),
                init_vector);
        } else if (config["GD_type"].get<std::string>() == "Momentum") {
            error_per_gen = MomentumGradientDescent(problem, gd_max_step,
                config["GD_r_eps"].get<double>(), config["GD_alpha"].get<double>(),
                config["GD_beta"].get<double>(),
                init_vector);
        } else if (config["GD_type"].get<std::string>() == "RMSprop") {
            error_per_gen = RMSprop(problem, gd_max_step,
                    config["GD_r_eps"].get<double>(),
                    config["GD_learning_rate"].get<double>(),
                    config["GD_ema_coef"].get<double>(),
                    init_vector);
        } else if (config["GD_type"].get<std::string>() == "Adam") {
            error_per_gen = Adam(problem, gd_max_step,
                    config["GD_r_eps"].get<double>(),
                    config["GD_learning_rate"].get<double>(),
                    config["GD_beta1"].get<double>(),
                    config["GD_beta2"].get<double>(),
                    init_vector);
        } else {
            throw("gd_call: unknown gradient descent type");
        }
    }
}

/**
 * @brief MultipleStarts script
 *
 *
 */
template <typename Problem>
std::vector<double> multipleStartsScript(json & config, Problem & problem)
{
    std::unique_ptr<BaseCoreGD> gd;
    auto string_gd_type = config["GD_type"].get<std::string>();
    int max_steps = config["GD_max_step"].get<int>();
    double r_eps = config["GD_r_eps"].get<double>();
    double learning_rate = config["GD_learning_rate"].get<double>();
    if (string_gd_type == "SimpleGD")
        gd = std::make_unique<CoreSimpleGradientDescent>(max_steps, r_eps, learning_rate);
    else if (string_gd_type == "Adam")
        gd = std::make_unique<CoreAdam>(max_steps, r_eps, learning_rate,
            config["GD_beta1"].get<double>(),
            config["GD_beta2"].get<double>());
    else
        throw(std::logic_error("Unknown gradient descent type for multiple starts run"));
    auto res = MultipleStartsGD(problem, *gd,
            config["MultipleStarts_number"].get<int>(),
            config["SobolStartIndex"].get<int>());
    return res;
}

/**
 * @brief Script to find model parameters
 *
 *
 * Optimizer is set in @p config
 *
 * @param[in] config Json config
 */
void script_general_optimizer(json & config)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    pcg_extras::seed_seq_from<std::random_device> seed_source;
    //const int seed_source = 42; /// @todo Move it to config file

    //MaleckarModel model; /// @todo Move it to config file
    KernikClancyModel model;

    ODESolver solver;

    ODEoptimization problem(model, solver);

    try {
        problem.read_config(config);
        problem.read_baselines(config);
    } catch (const std::string & err) {
        std::cerr << "Cannot read config for the problem: " << std::endl;
        std::cerr << err << std::endl;
        std::cerr << "Optimization will not be run" << std::endl;
        return;
    }

    std::vector<std::pair<int, double>> error_per_gen;

    const std::string sname = config["script"].get<std::string>();

    if (sname == "Genetic Algorithm") {
        if (mpi_rank == 0)
            std::cout << "Genetic algorithm starting now" << std::endl;

        gen_algo_call(seed_source, problem, config, error_per_gen);

        if (mpi_rank == 0)
            std::cout << "Genetic algorithm is complete" << std::endl;
    } else if (sname == "PSO") {
        if (mpi_rank == 0)
            std::cout << "PSO starting now" << std::endl;

        PSO_call(seed_source, problem, config, error_per_gen);

        if (mpi_rank == 0)
            std::cout << "PSO is complete" << std::endl;
    } else if (sname == "Nelder Mead") {
        if (mpi_rank == 0)
            std::cout << "Nelder-Mead starting now" << std::endl;

        nelder_mead_call(problem, config, error_per_gen);

        if (mpi_rank == 0)
            std::cout << "Nelder-Mead is complete" << std::endl;
    } else if (sname == "Gradient Descent") {
         if (mpi_rank == 0)
            std::cout << "Gradient Descent starting now" << std::endl;

        gd_call(problem, config, error_per_gen);

        if (mpi_rank == 0)
            std::cout << "Gradient Descent is complete" << std::endl;
    } else if (sname == "NM-GD") {
        if (mpi_rank == 0)
            std::cout << "Nelder-Mead and Gradient Descent starting now" << std::endl;

        nelder_mead_call(problem, config, error_per_gen);
        gd_call(problem, config, error_per_gen, problem.get_results_optimizer_format());

        if (mpi_rank == 0)
            std::cout << "Nelder-Mead and Gradient Descent is complete" << std::endl;

    } else if (sname == "MultipleStarts") {
        std::vector<double> res = multipleStartsScript(config, problem);
        if (mpi_rank == 0)
            dump_table_ode_problem(problem, res, config["MS_output_filename"].get<std::string>());
        return;
    } else {
        throw(std::logic_error("Unknown optimizer type"));
    }
    if (mpi_rank == 0) {
        std::cout << "Saving output files to disk..." << std::endl;
        std::string convergence_filename = "convergence_GA.txt";
        std::ofstream file(convergence_filename);
        for (const auto & p: error_per_gen)
            file << p.first << " " << p.second << std::endl;

        if (sname == "Genetic Algorithm")
            problem.export_gen_algo_tables();

        problem.dump_ap(problem.get_results_optimizer_format(), 0);

        std::cout << "USE THE FOLLOWING OUTPUT ONLY AS SOME HINT OF THE FINAL RESULT" << std::endl;
        using Results = decltype(problem)::Results;
        using BaselineResult = decltype(problem)::BaselineResult;
        Results results = problem.get_relative_results();
        BaselineResult res = results[0];
        std::cout << "Printing relative to default parameters" << std::endl;
        for (const auto & cit: res.constantsResult)
            std::cout << cit.first << " " << cit.second << std::endl;
        std::cout << std::endl << "Printing states relative to some baseline" << std::endl;
        for (const auto & sit: res.statesResult)
            std::cout << sit.first << " " << sit.second << std::endl;
    }
}



void script_direct_problem(json & config)
{
  //  MaleckarModel model;
    KernikClancyModel model;

    ODESolver solver;
    ODEoptimization problem(model, solver);

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
    CellChainModel<KernikClancyModel> model(config["chain_len"].get<int>(), config["main_cell_index"].get<int>(), config["gj_conduct"].get<double>());
    ODESolver solver;
    ODEoptimization problem(model, solver);

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
///@todo
}

void script_test_function(json & config)
{
    ///@todo
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
        std::cout << "Problem with a test function will be solved now" << std::endl;

    pcg_extras::seed_seq_from<std::random_device> seed_source;
   // int seed_source = 42;
    //using FuncToOptimize = SphereFunction<2>;
    //std::vector<double> init_guess(10, 1);
    //using FuncToOptimize = RosenbrockFunction<10>;
    //std::vector<double> init_guess(10, 3);


    //This guy is really the worst
    using FuncToOptimize = RastriginFunction<2>;
    std::vector<double> init_guess(10, 3);

    //using FuncToOptimize = StyblinskiTangFunction<10>;
    //std::vector<double> init_guess(10, 0);

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

    std::vector<std::pair<int, double>> error_per_gen;
    //nelder_mead_call(optim, config, error_per_gen, init_guess);
    gd_call(optim, config, error_per_gen, init_guess);
    if (rank == 0) {
        auto res = optim.get_result();
        std::cout << "Parameter error: " << func.solution_error(res) << std::endl;
    }
}



void script_track_minimum(json & config)
{
    ///@todo
#if 0
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

#endif
}

void script_mcmc(json & config)
{
    ///@todo
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
        std::cerr << "This is a build with no MCMC support" << std::endl;
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
        if (sname == "Direct Problem") {
            script_direct_problem(config);
        } else if (sname == "Direct Problem Chain") {
			script_direct_problem_chain(config);
        } else if (sname == "Test Function") {
            script_test_function(config);
        } else if (sname == "Track Minimum") {
            script_track_minimum(config);
        } else if (sname == "MCMC") {
            script_mcmc(config);
        } else {
            script_general_optimizer(config);
        }
    } catch(const std::string & e) {
        std::cerr << "Node " << mpi_rank << ": Exception string caught in main:\n" << e << std::endl;
    } catch(const char * s) {
        std::cerr << s << std::endl;
    }
    MPI_Finalize();
    return 0;
}
