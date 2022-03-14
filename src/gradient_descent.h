#ifndef GRADIENT_DESCENT
#define GRADIENT_DESCENT

#include <mpi.h>
#include <omp.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <cassert>
#include <algorithm>
#include "penalty.h"
#include "sobol_sequence.h"
#include "utils.h"

template <typename OptimizationProblem>
std::vector<std::pair<int, double>> simpleGradientDescent(OptimizationProblem & problem, int max_steps, double r_eps, double learning_rate, std::vector<double> init_vector = std::vector<double>())
{
    std::vector<std::pair<int, double>> error_per_gen;
    int param_num = problem.get_number_parameters();
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    std::vector<double> sol(param_num);
    if (init_vector.size() == 0) {
        init_vector.resize(param_num);
        int init_status = problem.initial_guess_for_optimizer(init_vector.begin());
    }
    if (init_vector.size() != param_num)
        throw("MomentumGradientDescent: init_vector is incorrect");



/* this is not really what we need because gradient descent should start exactly
 * from the point provided as an initial guess by the problem
 * It is not optimization method's business where is this guess from
    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j]) {
            sol[j] = std::uniform_real_distribution<double>(min_v[j], max_v[j])(rg);
        } else {
            assert(init_status != -1);
            sol[j] = init_vector[j];
        }
    }
*/

    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j]) {
            if (not (init_vector[j] >= min_v[j] && init_vector[j] <= max_v[j]))
                throw("simpleGradientDescent: initial guess is out of boundaries");
        } else {
            assert(init_status != -1);
        }
        sol[j] = init_vector[j];
    }

    std::vector<int> mut_pos;
    for (int i = 0 ; i < param_num; i++) {
        if (is_mutation_applicable[i])
            mut_pos.push_back(i);
    }

    std::vector<double> df(mut_pos.size()), dfdx(mut_pos.size());
    for (int step = 0; step < max_steps; step++) {
        //const double f = fitn(problem, sol, is_mutation_applicable, min_v, max_v);
        const double f = problem.get_objective_value(sol.begin());
        std::cout << step << " f: " << f << std::endl;
        error_per_gen.push_back({step, f});
        #pragma omp parallel for
        for (size_t i = 0; i < mut_pos.size(); i++) {
            /*
            std::vector<double> params = sol;
            const double eps = r_eps * std::abs(params[mut_pos[i]]);
            params[mut_pos[i]] += eps ;
            df[i] = fitn(problem, params, is_mutation_applicable, min_v, max_v) - f;
            dfdx[i] = df[i] / eps;
            */
            std::vector<double> params_fwd = sol, params_back = sol;
            double eps = r_eps * std::max(std::abs(sol[mut_pos[i]]), 1e-12); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = problem.get_objective_value(params_fwd.begin()) -
                    problem.get_objective_value(params_back.begin());
            //df[i] = fitn(problem, params_fwd, is_mutation_applicable, min_v, max_v)
            //       -
            //        fitn(problem, params_back, is_mutation_applicable, min_v, max_v);
            dfdx[i] = df[i] / (2*eps);
        }

        for (size_t i = 0; i < mut_pos.size(); i++) {
            const int pos = mut_pos[i];
            double ss = sol[pos] - learning_rate * dfdx[i];
            if (min_v[pos] > ss || max_v[pos] < ss) {
                std::cout << "out of boundary" << std::endl;
                std::cout << "sol[" << pos << "]: " << sol[pos] << std::endl;
                std::cout << "dfdx[" << i << "]: " << dfdx[i] << std::endl;
                std::cout << "ss[" << pos << "]: " << ss << std::endl;
                std::cout << "min max: " << min_v[pos] << " " << max_v[pos] <<  std::endl;
            }
            if (!std::isnan(ss))
                sol[pos] = ss;
        }
    }

    problem.submit_result(sol);
    return error_per_gen;
}


template <typename OptimizationProblem>
std::vector<std::pair<int, double>> MomentumGradientDescent(OptimizationProblem & problem, int max_steps, double r_eps, double alpha, double beta, std::vector<double> init_vector = std::vector<double>())
{
    std::vector<std::pair<int, double>> error_per_gen;
    int param_num = problem.get_number_parameters();
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    std::vector<double> sol(param_num);
    if (init_vector.size() == 0) {
        init_vector.resize(param_num);
        int init_status = problem.initial_guess_for_optimizer(init_vector.begin());
    }
    if (init_vector.size() != param_num)
        throw("MomentumGradientDescent: init_vector is incorrect");


/* this is not really what we need because gradient descent should start exactly
 * from the point provided as an initial guess by the problem
 * It is not optimization method's business where is this guess from
    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j]) {
            sol[j] = std::uniform_real_distribution<double>(min_v[j], max_v[j])(rg);
        } else {
            assert(init_status != -1);
            sol[j] = init_vector[j];
        }
    }
*/

    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j]) {
            if (not (init_vector[j] >= min_v[j] && init_vector[j] <= max_v[j]))
                throw("simpleGradientDescent: initial guess is out of boundaries");
        } else {
            assert(init_status != -1);
        }
        sol[j] = init_vector[j];
    }

    std::vector<int> mut_pos;
    for (int i = 0 ; i < param_num; i++) {
        if (is_mutation_applicable[i])
            mut_pos.push_back(i);
    }

    std::vector<double> df(mut_pos.size()), dfdx(mut_pos.size()), v_prev(mut_pos.size());
    for (int step = 0; step < max_steps; step++) {
        //const double f = fitn(problem, sol, is_mutation_applicable, min_v, max_v);
        const double f = problem.get_objective_value(sol.begin());
        std::cout << step << " f: " << f << std::endl;
        error_per_gen.push_back({step, f});
        #pragma omp parallel for
        for (size_t i = 0; i < mut_pos.size(); i++) {
            /*
            std::vector<double> params = sol;
            const double eps = r_eps * std::abs(params[mut_pos[i]]);
            params[mut_pos[i]] += eps ;
            df[i] = fitn(problem, params, is_mutation_applicable, min_v, max_v) - f;
            dfdx[i] = df[i] / eps;
            */
            std::vector<double> params_fwd = sol, params_back = sol;
            double eps = r_eps * std::max(std::abs(sol[mut_pos[i]]), 1e-12); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = problem.get_objective_value(params_fwd.begin()) -
                    problem.get_objective_value(params_back.begin());
            dfdx[i] = df[i] / (2*eps);
        }

        for (size_t i = 0; i < mut_pos.size(); i++) {
            const int pos = mut_pos[i];
            double v = beta * v_prev[i] - dfdx[i];
            double x_new = sol[pos] + alpha * v;
            if (min_v[pos] > x_new || max_v[pos] < x_new) {
                std::cout << "out of boundary" << std::endl;
                std::cout << "sol[" << pos << "]: " << sol[pos] << std::endl;
                std::cout << "dfdx[" << i << "]: " << dfdx[i] << std::endl;
                std::cout << "ss[" << pos << "]: " << x_new << std::endl;
                std::cout << "min max: " << min_v[pos] << " " << max_v[pos] <<  std::endl;
            }
            if (!std::isnan(x_new)) {
                sol[pos] = x_new;
                v_prev[i] = v;
            }
        }
    }

    problem.submit_result(sol);
    return error_per_gen;
}


template <typename OptimizationProblem>
std::vector<std::pair<int, double>> RMSprop(OptimizationProblem & problem, int max_steps, double r_eps, double learning_rate, double ema_coef, std::vector<double> init_vector = std::vector<double>())
{
    std::vector<std::pair<int, double>> error_per_gen;
    int param_num = problem.get_number_parameters();
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    std::vector<double> sol(param_num);
    if (init_vector.size() == 0) {
        init_vector.resize(param_num);
        int init_status = problem.initial_guess_for_optimizer(init_vector.begin());
    }
    if (init_vector.size() != param_num)
        throw("MomentumGradientDescent: init_vector is incorrect");

    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j]) {
            if (not (init_vector[j] >= min_v[j] && init_vector[j] <= max_v[j]))
                throw("simpleGradientDescent: initial guess is out of boundaries");
        } else {
            assert(init_status != -1);
        }
        sol[j] = init_vector[j];
    }

    std::vector<int> mut_pos;
    for (int i = 0 ; i < param_num; i++) {
        if (is_mutation_applicable[i])
            mut_pos.push_back(i);
    }

    std::vector<double> df(mut_pos.size()), dfdx(mut_pos.size());
    std::vector<double> ema_sq_grad(mut_pos.size()); // squared gradient exponential moving average
    for (int step = 0; step < max_steps; step++) {
        //const double f = fitn(problem, sol, is_mutation_applicable, min_v, max_v);
        const double f = problem.get_objective_value(sol.begin());
        std::cout << step << " f: " << f << std::endl;
        error_per_gen.push_back({step, f});
        #pragma omp parallel for
        for (size_t i = 0; i < mut_pos.size(); i++) {
            /*
            std::vector<double> params = sol;
            const double eps = r_eps * std::abs(params[mut_pos[i]]);
            params[mut_pos[i]] += eps ;
            df[i] = fitn(problem, params, is_mutation_applicable, min_v, max_v) - f;
            dfdx[i] = df[i] / eps;
            */
            std::vector<double> params_fwd = sol, params_back = sol;
            double eps = r_eps * std::max(std::abs(sol[mut_pos[i]]), 1e-12); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = problem.get_objective_value(params_fwd.begin()) -
                    problem.get_objective_value(params_back.begin());
            dfdx[i] = df[i] / (2*eps);
        }

        for (size_t i = 0; i < mut_pos.size(); i++) {
            const int pos = mut_pos[i];
            ema_sq_grad[i] = ema(std::pow(dfdx[i], 2), ema_sq_grad[i], ema_coef);
            double x_new = sol[pos] - learning_rate * dfdx[i] / (std::sqrt(ema_sq_grad[i]) + 1e-8);
            if (min_v[pos] > x_new || max_v[pos] < x_new) {
                std::cout << "out of boundary" << std::endl;
                std::cout << "sol[" << pos << "]: " << sol[pos] << std::endl;
                std::cout << "dfdx[" << i << "]: " << dfdx[i] << std::endl;
                std::cout << "ss[" << pos << "]: " << x_new << std::endl;
                std::cout << "min max: " << min_v[pos] << " " << max_v[pos] <<  std::endl;
            }
            if (!std::isnan(x_new)) {
                sol[pos] = x_new;
            }
        }
    }

    problem.submit_result(sol);
    return error_per_gen;
}



template <typename OptimizationProblem>
std::vector<std::pair<int, double>> Adam(OptimizationProblem & problem, int max_steps, double r_eps, double learning_rate, double beta1, double beta2, std::vector<double> init_vector = std::vector<double>())
{
    std::vector<std::pair<int, double>> error_per_gen;
    int param_num = problem.get_number_parameters();
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    std::vector<double> sol(param_num);
    if (init_vector.size() == 0) {
        init_vector.resize(param_num);
        int init_status = problem.initial_guess_for_optimizer(init_vector.begin());
    }
    if (init_vector.size() != param_num)
        throw("MomentumGradientDescent: init_vector is incorrect");

    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j]) {
            if (not (init_vector[j] >= min_v[j] && init_vector[j] <= max_v[j]))
                throw("simpleGradientDescent: initial guess is out of boundaries");
        } else {
            assert(init_status != -1);
        }
        sol[j] = init_vector[j];
    }

    std::vector<int> mut_pos;
    for (int i = 0 ; i < param_num; i++) {
        if (is_mutation_applicable[i])
            mut_pos.push_back(i);
    }

    std::vector<double> df(mut_pos.size()), dfdx(mut_pos.size());
    std::vector<double> ema_sq_grad(mut_pos.size()); // squared gradient exponential moving average
    std::vector<double> ema_grad(mut_pos.size()); // ema for gradient
    for (int step = 0; step < max_steps; step++) {
        //const double f = fitn(problem, sol, is_mutation_applicable, min_v, max_v);
        const double f = problem.get_objective_value(sol.begin());
        std::cout << step << " f: " << f << std::endl;
        error_per_gen.push_back({step, f});
        #pragma omp parallel for
        for (size_t i = 0; i < mut_pos.size(); i++) {

            std::vector<double> params = sol;
            const double eps = r_eps * std::max(std::abs(params[mut_pos[i]]), 1e-12);
            params[mut_pos[i]] += eps;
            df[i] = problem.get_objective_value(params.begin()) - f;
            dfdx[i] = df[i] / eps;

            /*
            std::vector<double> params_fwd = sol, params_back = sol;
            double eps = r_eps * std::max(std::abs(sol[mut_pos[i]]), 1e-12); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = problem.get_objective_value(params_fwd.begin()) -
                    problem.get_objective_value(params_back.begin());
            dfdx[i] = df[i] / (2*eps);
            */
        }

        for (size_t i = 0; i < mut_pos.size(); i++) {
            const int pos = mut_pos[i];
            ema_grad[i] = ema(dfdx[i], ema_grad[i], beta1);
            ema_sq_grad[i] = ema(std::pow(dfdx[i], 2), ema_sq_grad[i], beta2);

            double corrected_ema_grad = ema_grad[i] / (1 - std::pow(beta1, step + 1));
            double corrected_ema_sq_grad = ema_sq_grad[i] / (1 - std::pow(beta2, step + 1));

            double x_new = sol[pos] - learning_rate * corrected_ema_grad / (std::sqrt(corrected_ema_sq_grad) + 1e-8);
            if (min_v[pos] > x_new || max_v[pos] < x_new) {
                std::cout << "out of boundary" << std::endl;
                std::cout << "sol[" << pos << "]: " << sol[pos] << std::endl;
                std::cout << "dfdx[" << i << "]: " << dfdx[i] << std::endl;
                std::cout << "ss[" << pos << "]: " << x_new << std::endl;
                std::cout << "min max: " << min_v[pos] << " " << max_v[pos] <<  std::endl;
            }
            if (!std::isnan(x_new)) {
                sol[pos] = x_new;
            }
        }
    }

    problem.submit_result(sol);
    return error_per_gen;
}


template <typename OptimizationProblem>
std::vector<double> CoreAdam(OptimizationProblem & problem, int max_steps, double r_eps, double learning_rate, double beta1, double beta2, const std::vector<double> & initial_guess, const std::vector<int> & is_mutation_applicable)
{
    //store vector and loss value on each step
    std::vector<double> res;
    int param_num = initial_guess.size();
    res.reserve((param_num + 1) * max_steps);

    std::vector<int> mut_pos;
    for (int i = 0 ; i < param_num; i++) {
        if (is_mutation_applicable[i])
            mut_pos.push_back(i);
    }
    std::vector<double> sol = initial_guess;

    std::vector<double> df(mut_pos.size()), dfdx(mut_pos.size());
    std::vector<double> ema_sq_grad(mut_pos.size()); // squared gradient exponential moving average
    std::vector<double> ema_grad(mut_pos.size()); // ema for gradient
    for (int step = -1; step < max_steps; step++) {
        //const double f = fitn(problem, sol, is_mutation_applicable, min_v, max_v);
        const double f = problem.get_objective_value(sol.begin());
        if (step == -1) //warmup
            continue;
        //std::cout << step << " f: " << f << std::endl;
        res.insert(res.end(), sol.begin(), sol.end());
        res.insert(res.end(), f);
        #pragma omp parallel for
        for (size_t i = 0; i < mut_pos.size(); i++) {

            std::vector<double> params = sol;
            const double eps = r_eps * std::max(std::abs(params[mut_pos[i]]), 1e-12);
            params[mut_pos[i]] += eps;
            df[i] = problem.get_objective_value(params.begin()) - f;
            dfdx[i] = df[i] / eps;

            /*
            std::vector<double> params_fwd = sol, params_back = sol;
            double eps = r_eps * std::max(std::abs(sol[mut_pos[i]]), 1e-12); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = problem.get_objective_value(params_fwd.begin()) -
                    problem.get_objective_value(params_back.begin());
            dfdx[i] = df[i] / (2*eps);
            */
        }

        for (size_t i = 0; i < mut_pos.size(); i++) {
            const int pos = mut_pos[i];
            ema_grad[i] = ema(dfdx[i], ema_grad[i], beta1);
            ema_sq_grad[i] = ema(std::pow(dfdx[i], 2), ema_sq_grad[i], beta2);

            double corrected_ema_grad = ema_grad[i] / (1 - std::pow(beta1, step + 1));
            double corrected_ema_sq_grad = ema_sq_grad[i] / (1 - std::pow(beta2, step + 1));

            double x_new = sol[pos] - learning_rate * corrected_ema_grad / (std::sqrt(corrected_ema_sq_grad) + 1e-12);
            /*
            if (min_v[pos] > x_new || max_v[pos] < x_new) {
                std::cout << "out of boundary" << std::endl;
                std::cout << "sol[" << pos << "]: " << sol[pos] << std::endl;
                std::cout << "dfdx[" << i << "]: " << dfdx[i] << std::endl;
                std::cout << "ss[" << pos << "]: " << x_new << std::endl;
                std::cout << "min max: " << min_v[pos] << " " << max_v[pos] <<  std::endl;
            }
            */
            if (!std::isnan(x_new)) {
                sol[pos] = x_new;
            }
        }
    }

    return res;
}

template <typename OptimizationProblem>
std::vector<double> MultipleStartsAdam(OptimizationProblem & problem, int max_steps, double r_eps, double learning_rate, double beta1, double beta2, int starts_number, int sobol_start_index = 0)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    int param_num = problem.get_number_parameters();
    std::vector<double> min_v(param_num), max_v(param_num);

    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    std::vector<int> mutable_params_indices;
    for (int i = 0; i < is_mutation_applicable.size(); i++) {
        if (is_mutation_applicable[i])
            mutable_params_indices.push_back(i);
    }

    std::vector<double> main_tight_min(param_num),
                        main_tight_max(param_num);
    problem.get_tight_boundaries(main_tight_min, main_tight_max);

    std::vector<double> tight_min(mutable_params_indices.size()),
                        tight_max(mutable_params_indices.size());
    for (int i = 0; i < mutable_params_indices.size(); i++) {
        tight_min[i] = main_tight_min[mutable_params_indices[i]];
        tight_max[i] = main_tight_max[mutable_params_indices[i]];
    }

    Sobol sobol01(mutable_params_indices.size());

    std::vector<double> init_vector(param_num);
    int init_status = problem.initial_guess_for_optimizer(init_vector.begin());

    std::vector<double> total_res;
    total_res.reserve(starts_number * (param_num + 1) * max_steps);

    //each mpi process has its own sobol subsequence
    //skip first sobol_start_index elements since they could have been checked
    //in previous runs of the program
    for (int i = 0; i < sobol_start_index + starts_number * mpi_rank; i++)
        sobol01();

    for (int i = 0; i < starts_number; i++) {
        if (mpi_rank == 0)
            std::cout << i + 1 << "/" << starts_number << std::endl;
        auto mutable_params = sobol01();
        scale_vector(mutable_params, tight_min, tight_max);
        // now fill params_complete with params
        auto params_complete = init_vector;
        for (int i = 0; i < mutable_params_indices.size(); i++)
            params_complete[mutable_params_indices[i]] = mutable_params[i];
        auto res = CoreAdam(problem, max_steps, r_eps, learning_rate, beta1, beta2, params_complete, is_mutation_applicable);
        total_res.insert(total_res.end(), res.begin(), res.end());
    }
    if (mpi_rank == 0)
        total_res.resize(mpi_size * total_res.size());
    MPI_Gather(
    (mpi_rank == 0) ? MPI_IN_PLACE : total_res.data(),
    starts_number * (param_num + 1) * max_steps,
    MPI_DOUBLE,
    total_res.data(),
    starts_number * (param_num + 1) * max_steps,
    MPI_DOUBLE,
    0,
    MPI_COMM_WORLD);

    return total_res;
}

/*

template <typename InitializedRandomGenerator, typename OptimizationProblem>
std::vector<std::pair<int, double>> weirdSteepestGradientDescent(InitializedRandomGenerator rg, OptimizationProblem & problem, int max_steps, double r_eps)
{
    int param_num = problem.get_number_parameters();
    std::vector<double> init_vector(param_num);
    int init_status = problem.initial_guess(init_vector.begin());
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    if (init_status == -1) {
        throw ("initial guess is required for gradient descent");
    }
    return weirdSteepestGradientDescent(rg, problem, max_steps, r_eps, init_vector);
}


template <typename InitializedRandomGenerator, typename OptimizationProblem>
std::vector<std::pair<int, double>> weirdSteepestGradientDescentTrack(InitializedRandomGenerator rg, OptimizationProblem & problem, double r_eps, int global_steps, int local_steps)
{
    int param_num = problem.get_number_parameters();
    std::vector<double> init_vector(param_num);
    int init_status = problem.initial_guess(init_vector.begin());
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    if (init_status == -1) {
        std::cout << "NOOOOOOOOOOOOOOOO" << std::endl;
        for (int j = 0; j < param_num; j++) {
            if (is_mutation_applicable[j])
                init_vector[j] = std::uniform_real_distribution<double>(min_v[j], max_v[j])(rg);
        }
    }
    return weirdSteepestGradientDescentTrack(rg, problem, r_eps, global_steps, local_steps, init_vector);
}

template <typename InitializedRandomGenerator, typename OptimizationProblem>
std::vector<std::pair<int, double>> weirdSteepestGradientDescentTrack(InitializedRandomGenerator rg,
                            OptimizationProblem & problem, double r_eps, int global_steps,
                            int local_steps, const std::vector<double> & init_vector)
{
    std::vector<std::pair<int, double>> error_per_gen;
    int param_num = problem.get_number_parameters();
    std::cout << param_num << std::endl;
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    if (init_vector.size() != param_num)
        throw 1;
    std::vector<double> sol = init_vector;

    std::vector<int> mut_pos;
    for (int i = 0 ; i < param_num; i++) {
        if (is_mutation_applicable[i])
            mut_pos.push_back(i);
    }
    std::cout << mut_pos.size() << std::endl;

    problem.start_track(sol.begin());

    std::vector<double> df(mut_pos.size()), dfdx(mut_pos.size());
    std::cout << "Initial fitness: " << fitn(problem, sol, is_mutation_applicable, min_v, max_v) << std::endl;
    //global_steps = -1;////////////////////////////////////////////////////////////////////////////////////////////////

    //global cycle
    for (int global_step = 0; global_step <= global_steps; global_step++) {
        problem.set_alpha((double)global_step / global_steps);
        std::cout << global_step << "/" << global_steps << std::endl;
        //local cycle
        for (int step = 0; step < local_steps; step++) {
            const double f = fitn(problem, sol, is_mutation_applicable, min_v, max_v);
            error_per_gen.push_back({step, f});
           // std::cout << step << " f: " << f << std::endl;
            double grad_time = MPI_Wtime();
            #pragma omp parallel for
            for (int i = 0; i < mut_pos.size(); i++) {
                std::vector<double> params_fwd = sol, params_back = sol;
                double eps = r_eps * std::abs(sol[mut_pos[i]]); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
                params_fwd[mut_pos[i]] += eps;
                params_back[mut_pos[i]] -= eps;
                df[i] = fitn(problem, params_fwd, is_mutation_applicable, min_v, max_v)
                        -
                        fitn(problem, params_back, is_mutation_applicable, min_v, max_v);
                dfdx[i] = df[i] / (2*eps);
            }
            std::cout << "grad_time: " << MPI_Wtime() - grad_time << std::endl;
            std::cout << "fitness: " << f << std::endl;
            double best_step_time = MPI_Wtime();
            const int min_pow = 0, max_pow = 8;
            std::vector<double> best_power_f(max_pow - min_pow);
            #pragma omp parallel for
            for (int power = 0; power < best_power_f.size(); power++) {
                const double gd_step = std::pow(10, -(power + min_pow));
                std::vector<double> params = sol;
                for (int i = 0; i < mut_pos.size(); i++) {
                    const int pos = mut_pos[i];
                    params[pos] -= std::abs(params[pos]) * gd_step * dfdx[i];
                }
                best_power_f[power] = fitn(problem, params, is_mutation_applicable, min_v, max_v);
            }
            std::cout << "best_step_time: " << MPI_Wtime() - best_step_time << std::endl;

            int best_power = std::min_element(best_power_f.begin(), best_power_f.end()) - best_power_f.begin();
           // std::cout << "best power " << best_power + min_pow << "f : " << f << " new_f: " << best_power_f[best_power] << std::endl;
            if (f < best_power_f[best_power] + 1e-2 * std::abs(f)) break;

            const double gd_step = std::pow(10, -(best_power + min_pow));

            for (int i = 0; i < mut_pos.size(); i++) {
                const int pos = mut_pos[i];
                double ss = sol[pos] -  std::abs(sol[pos]) * gd_step * dfdx[i];
                if (min_v[pos] > ss || max_v[pos] < ss) {
                    std::cout << "df: " << df[i] << std::endl;
                    std::cout << sol[pos] << " " << dfdx[i] << " " << ss << " " << min_v[pos] << " " << max_v[pos] << std::endl;
                }
                sol[pos] = ss;
            }
        }
        problem.dump_ap(sol.begin(), global_step);
    }
    problem.is_AP_normalized = 0;
    problem.dump_ap(sol.begin(), global_steps);
    problem.submit_result(sol);
    return error_per_gen;
}



template <typename InitializedRandomGenerator, typename OptimizationProblem>
std::vector<std::pair<int, double>> weirdSteepestGradientDescent(InitializedRandomGenerator rg, OptimizationProblem & problem, int max_steps, double r_eps, const std::vector<double> & init_vector)
{
    std::vector<std::pair<int, double>> error_per_gen;
    int param_num = problem.get_number_parameters();
    std::cout << param_num << std::endl;
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    if (init_vector.size() != param_num)
        throw 1;
    std::vector<double> sol = init_vector;

    std::vector<int> mut_pos;
    for (int i = 0 ; i < param_num; i++) {
        if (is_mutation_applicable[i])
            mut_pos.push_back(i);
    }
    std::cout << mut_pos.size() << std::endl;

    std::vector<double> df(mut_pos.size()), dfdx(mut_pos.size());
    for (int step = 0; step < max_steps; step++) {
        const double f = fitn(problem, sol, is_mutation_applicable, min_v, max_v);
        std::cout << step << " f: " << f << std::endl;
        error_per_gen.push_back({step, f});
        #pragma omp parallel for
        for (int i = 0; i < mut_pos.size(); i++) {
            std::vector<double> params_fwd = sol, params_back = sol;
            //TODO
            const double eps = r_eps * std::abs(sol[mut_pos[i]]); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            //const double eps = r_eps * (max_v[mut_pos[i]] - min_v[mut_pos[i]]);
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = fitn(problem, params_fwd, is_mutation_applicable, min_v, max_v)
                    -
                    fitn(problem, params_back, is_mutation_applicable, min_v, max_v);
            dfdx[i] = df[i] / (2*eps);
        }

        const int min_pow = -2, max_pow = 5;
        std::vector<double> best_power_f(max_pow - min_pow);
        #pragma omp parallel for
        for (int power = 0; power < best_power_f.size(); power++) {
            const double gd_step = std::pow(10, -(power + min_pow));
            std::vector<double> params = sol;
            for (int i = 0; i < mut_pos.size(); i++) {
                const int pos = mut_pos[i];
                //TODO
                params[pos] -= std::abs(params[pos]) * gd_step * dfdx[i];
                //params[pos] -= (max_v[pos] - min_v[pos]) * gd_step * dfdx[i];
            }
            best_power_f[power] = fitn(problem, params, is_mutation_applicable, min_v, max_v);
        }

        int best_power = std::min_element(best_power_f.begin(), best_power_f.end()) - best_power_f.begin();
       // std::cout << "best power " << best_power + min_pow << "f : " << f << " new_f: " << best_power_f[best_power] << std::endl;
        if (f < best_power_f[best_power]) continue;/////////////////////////////////////////////////////////////////////////////////////////////////////////

        const double gd_step = std::pow(10, -(best_power + min_pow));

        for (int i = 0; i < mut_pos.size(); i++) {
            const int pos = mut_pos[i];
            //TODO
            double ss = sol[pos] -  std::abs(sol[pos]) * gd_step * dfdx[i];

           // double ss = sol[pos] -  (max_v[pos] - min_v[pos]) * gd_step * dfdx[i];
            if (min_v[pos] > ss || max_v[pos] < ss) {
                std::cout << "df: " << df[i] << std::endl;
                std::cout << sol[pos] << " " << dfdx[i] << " " << ss << " " << min_v[pos] << " " << max_v[pos] << std::endl;
            }
            sol[pos] = ss;
        }
    }
    problem.dump_ap(sol.begin(), 10);

    problem.genetic_algorithm_result(sol);
    return error_per_gen;
}
*/
#endif
