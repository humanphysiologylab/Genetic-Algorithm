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

template <typename InitializedRandomGenerator, typename OptimizationProblem>
std::vector<std::pair<int, double>> simpleGradientDescent(InitializedRandomGenerator rg, OptimizationProblem & problem, int max_steps, double r_eps)
{
    std::vector<std::pair<int, double>> error_per_gen;
    int param_num = problem.get_number_parameters();
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);

    std::vector<double> init_vector(param_num), sol(param_num);
    int init_status = problem.initial_guess(init_vector.begin());


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
            assert(init_vector[j] >= min_v[j] && init_vector[j] <= max_v[j]);
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
        const double f = fitn(problem, sol, is_mutation_applicable, min_v, max_v);
        std::cout << step << " f: " << f << std::endl;
        error_per_gen.push_back({step, f});
        #pragma omp parallel for
        for (int i = 0; i < mut_pos.size(); i++) {
            /*
            std::vector<double> params = sol;
            const double eps = r_eps * std::abs(params[mut_pos[i]]);
            params[mut_pos[i]] += eps ;
            df[i] = fitn(problem, params, is_mutation_applicable, min_v, max_v) - f;
            dfdx[i] = df[i] / eps;
            */
            std::vector<double> params_fwd = sol, params_back = sol;
            //TODO what if sol[mut_pos[i]] is 0 or close to 0?
            double eps = r_eps * std::abs(sol[mut_pos[i]]); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = fitn(problem, params_fwd, is_mutation_applicable, min_v, max_v)
                    -
                    fitn(problem, params_back, is_mutation_applicable, min_v, max_v);
            dfdx[i] = df[i] / (2*eps);
        }

        const double gd_step = 0.005;//TODO maybe put it in config

        for (int i = 0; i < mut_pos.size(); i++) {
            const int pos = mut_pos[i];
            double ss = sol[pos] - gd_step * std::abs(sol[pos]) * dfdx[i];
            if (min_v[pos] > ss || max_v[pos] < ss) {
                std::cout << "df: " << df[i] << std::endl;
                std::cout << sol[pos] << " " << dfdx[i] << " " << ss << " " << min_v[pos] << " " << max_v[pos] << std::endl;
            }
            sol[pos] = ss;
        }
    }

    problem.genetic_algorithm_result(sol);
    return error_per_gen;
}




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
    problem.genetic_algorithm_result(sol);
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
#endif
