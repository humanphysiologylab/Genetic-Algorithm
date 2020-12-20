#ifndef BASIC_POPULATION
#define BASIC_POPULATION

#include <mpi.h>
#include <omp.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <cassert>
#include <algorithm>

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

    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j]) {
            sol[j] = std::uniform_real_distribution<double>(min_v[j], max_v[j])(rg);
        } else {
            assert(init_status != -1);
            sol[j] = init_vector[j];
        }
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
            double eps = r_eps * std::abs(sol[mut_pos[i]]); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = fitn(problem, params_fwd, is_mutation_applicable, min_v, max_v)
                    - 
                    fitn(problem, params_back, is_mutation_applicable, min_v, max_v);
            dfdx[i] = df[i] / (2*eps);
        }

        const double gd_step = 0.005;

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




template <typename OptimizationProblem>
double fitn(OptimizationProblem & problem, std::vector<double> & params, const std::vector<int> & is_mutation_applicable,
                const std::vector<double> & min_v, const std::vector<double> & max_v)
{
    double penalty = 0;
    for (int i = 0; i < is_mutation_applicable.size(); i++) {
        if (!is_mutation_applicable[i]) continue;
 
        if (min_v[i] > params[i])
            penalty += std::pow(20 * (params[i] - min_v[i]) / min_v[i], 2);
        if (max_v[i] < params[i])
            penalty += std::pow(20 * (params[i] - max_v[i]) / max_v[i], 2);
    }
    //penalty = 0;////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return penalty + problem.genetic_algorithm_calls(params.begin());
}



template <typename InitializedRandomGenerator, typename OptimizationProblem>
std::vector<std::pair<int, double>> weirdSteepestGradientDescent(InitializedRandomGenerator rg, OptimizationProblem & problem, int max_steps, double r_eps)
{
    int param_num = problem.get_number_parameters();
    std::vector<double> init_vector(param_num);
    problem.initial_guess(init_vector.begin());
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);
    
    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j])
            init_vector[j] = std::uniform_real_distribution<double>(min_v[j], max_v[j])(rg);
    }
    return weirdSteepestGradientDescent(rg, problem, max_steps, r_eps, init_vector);
}


template <typename InitializedRandomGenerator, typename OptimizationProblem>
std::vector<std::pair<int, double>> weirdSteepestGradientDescentTrack(InitializedRandomGenerator rg, OptimizationProblem & problem, double r_eps, int global_steps, int local_steps)
{
    int param_num = problem.get_number_parameters();
    std::vector<double> init_vector(param_num);
    problem.initial_guess(init_vector.begin());
    std::vector<double> min_v(param_num), max_v(param_num);
    std::vector<int> is_mutation_applicable(param_num);
    int boundaries_status = problem.get_boundaries(min_v, max_v, is_mutation_applicable);
    
    for (int j = 0; j < param_num; j++) {
        if (is_mutation_applicable[j])
            init_vector[j] = std::uniform_real_distribution<double>(min_v[j], max_v[j])(rg);
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
        throw;
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
            if (f < best_power_f[best_power] + 1e-2 * f) break;

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
        throw;
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
       // std::cout << step << " f: " << f << std::endl;
        error_per_gen.push_back({step, f});
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

        const int min_pow = -2, max_pow = 8;
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

        int best_power = std::min_element(best_power_f.begin(), best_power_f.end()) - best_power_f.begin();
       // std::cout << "best power " << best_power + min_pow << "f : " << f << " new_f: " << best_power_f[best_power] << std::endl;
        if (f < best_power_f[best_power]) continue;

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

    problem.genetic_algorithm_result(sol);
    return error_per_gen;
}


template <typename OptimizationProblem>
class BasicPopulation
{

public:
    OptimizationProblem & problem;

    int number_organisms;
    int number_elites;
    int number_mutants;

    int mpi_rank;
    int mpi_size;

    std::vector<double> min_gene, max_gene;
    std::vector<int> is_mutation_applicable;

    const int * get_is_mutation_applicable() const
    {
        return is_mutation_applicable.data();
    }
    
    std::vector<double> all_genes;
    int genes_per_organism;
    int get_genes_per_organism() const
    {
        return genes_per_organism;
    }

    double * elite_genes_buffer;
    int elite_genes_buffer_size;

    double * mutant_genes_buffer;
    int mutant_genes_buffer_size;

    std::vector<double> fitness_values;

    BasicPopulation(OptimizationProblem & problem, int number_elites_p, int number_mutants_p)
    : problem(problem),
      number_organisms(number_elites_p + number_mutants_p),
      number_elites(number_elites_p),
      number_mutants(number_mutants_p),
      genes_per_organism(problem.get_number_parameters())
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        if (number_mutants % 2 != 0) {
            //printf("Even number of mutants required!\nWe will try to add one more\n");
            number_mutants++;
            number_organisms++;
        }
        if (number_organisms % mpi_size != 0) {
            //printf("Number of nodes should divide number of organisms completely\n");
            //printf("We will add more mutants and elites for that\n");
            int additional_organisms = mpi_size - number_organisms % mpi_size;
            number_organisms += additional_organisms;
            number_mutants += additional_organisms;
            if (number_mutants % 2 != 0) {
                number_mutants--;
                number_elites++;
            }
        }


        min_gene.resize(genes_per_organism, nan(""));
        max_gene.resize(genes_per_organism, nan(""));
        is_mutation_applicable.resize(genes_per_organism, -1);

        int boundaries_status = problem.get_boundaries(min_gene, max_gene, is_mutation_applicable);

        if (boundaries_status == -1) {
            if (mpi_rank == 0)
                std::cerr << "non-constrained optimization problems are not supported"
                    << std::endl;
            throw;
        }

        fitness_values.resize(number_organisms,  nan(""));
        int all_genes_size = number_organisms * genes_per_organism;
        all_genes.resize(all_genes_size,  nan(""));

        elite_genes_buffer_size = number_elites * genes_per_organism;
        elite_genes_buffer = new double [elite_genes_buffer_size];

        mutant_genes_buffer_size = number_mutants * genes_per_organism;
        mutant_genes_buffer = new double [mutant_genes_buffer_size];

        for (int i = 0; i < elite_genes_buffer_size; i++)
            elite_genes_buffer[i] = nan("");
        for (int i = 0; i < mutant_genes_buffer_size; i++)
            mutant_genes_buffer[i] = nan("");
    }

    ~BasicPopulation()
    {
        delete [] elite_genes_buffer;
        delete [] mutant_genes_buffer;
    }

    template<typename InitializedRandomGenerator>
    void init(InitializedRandomGenerator rg)
    {
        std::vector<double> init_vector(genes_per_organism, nan(""));
        int init_status = problem.initial_guess(init_vector.begin());
        
        for (int i = 0; i < number_organisms; i++) {
            for (int j = 0; j < genes_per_organism; j++) {
                if (is_mutation_applicable[j]) {
                    all_genes[j + i * genes_per_organism] =
                        std::uniform_real_distribution<double>(min_gene[j], max_gene[j])(rg);
                } else {
                    assert(init_status != -1);
                    //or maybe it is fine to put zero if init_status == -1
                    all_genes[j + i * genes_per_organism] = init_vector[j];
                }
            }
        }
    }

    void scatter()
    {
        MPI_Scatter(all_genes.data(), all_genes.size() / mpi_size,
                    MPI_DOUBLE, (mpi_rank == 0) ? MPI_IN_PLACE : all_genes.data(),
                    all_genes.size() / mpi_size,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void gather()
    {
        MPI_Gather((mpi_rank == 0) ? MPI_IN_PLACE : fitness_values.data(), number_organisms / mpi_size, MPI_DOUBLE,
                   fitness_values.data(), number_organisms / mpi_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gather((mpi_rank == 0) ? MPI_IN_PLACE : all_genes.data(), all_genes.size() / mpi_size,
                   MPI_DOUBLE, all_genes.data(),
                   all_genes.size() / mpi_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void run_generation()
    {
        #pragma omp parallel for
        for (int i = 0; i < number_organisms / mpi_size; i++) {
            fitness_values[i] = problem.genetic_algorithm_calls(all_genes.begin() + i * genes_per_organism);
        }
    }

    void fitness_function(std::vector<std::pair<double, int>> & sd_n_index)
    {
        for (int i = 0; i < number_organisms; i++) {
            sd_n_index[i].first = fitness_values[i];
            sd_n_index[i].second = i;
        }
    }

protected:

    void copy_genes(double * from, double * to)
    {
        for (int i = 0; i < genes_per_organism; i++)
            to[i] = from[i];
    }

public:
    void save_elite_to_elite_buffer(int from, int to)
    {
        copy_genes(all_genes.data() + from * genes_per_organism, elite_genes_buffer + to * genes_per_organism);
    }
    void save_mutant_to_mutant_buffer(int from, int to)
    {
        copy_genes(all_genes.data() + from * genes_per_organism, mutant_genes_buffer + to * genes_per_organism);
    }
    void restore_elites_to_main_array()
    {
        for (int i = 0; i < number_elites; i++)
            copy_genes(elite_genes_buffer + i * genes_per_organism, all_genes.data() + i * genes_per_organism);
    }
    void restore_mutants_to_main_array()
    {
        for (int i = 0; i < number_mutants; i++)
            copy_genes(mutant_genes_buffer + i * genes_per_organism, all_genes.data() + (number_elites + i) * genes_per_organism);
    }
    std::vector<double> best() const
    {
        return std::vector<double>(all_genes.begin(), all_genes.begin() + genes_per_organism);
    }
    
    std::vector<std::pair<int, double>> error_per_gen;
    std::vector<std::pair<int, double>> get_error_per_gen() const
    {
        return error_per_gen;
    }
    void log(const std::vector<std::pair<double, int>> & sd_n_index, int gen)
    {
        error_per_gen.push_back({gen, sd_n_index[0].first});
        if (mpi_rank != 0) return;
        if (gen % 1 != 0) return ;
        std::cout << "Generation: " << gen << std::endl
                  << "Best: " << sd_n_index[0].first << std::endl
                  << "Worst: " << sd_n_index.back().first << std::endl;
        
       // auto best_genes = best();
      //  std::cout << "Genes:";
       // for (auto &g: best_genes)
       //     std::cout << " " << g;
       // std::cout << std::endl << std::endl;
    }
    double * get_mutant_buffer_genes()
    {
        return mutant_genes_buffer;
    }
    double * get_min_gene_value()
    {
        return min_gene.data();
    }
    double * get_max_gene_value()
    {
        return max_gene.data();
    }
    int get_number_genes()
    {
        return genes_per_organism;
    }
    void done()
    {
        if (mpi_rank != 0) return;
        problem.genetic_algorithm_result(all_genes);
    }
};

#endif
