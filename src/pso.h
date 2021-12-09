#ifndef PSO
#define PSO

#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <omp.h>
#include <vector>
#include <cmath>


template <typename OptimizationProblem, typename RandomGenerator, typename Seed>
class PSO_population
{

public:
    std::vector<RandomGenerator> random_generators;

    OptimizationProblem & problem;

    int number_organisms;

    int mpi_rank;
    int mpi_size;

    std::vector<double> min_gene, max_gene;
    std::vector<int> is_mutation_applicable;

    const int * get_is_mutation_applicable() const
    {
        return is_mutation_applicable.data();
    }

    std::vector<double> all_genes, all_best_genes, velocities;
    std::vector<double> best_global_genes;
    int genes_per_organism;
    int get_genes_per_organism() const
    {
        return genes_per_organism;
    }

    std::vector<double> all_best_genes_fitness_values;
    double best_global_genes_fitness_value;
    std::vector<double> init_vector;
    int init_status;
    PSO_population(OptimizationProblem & problem, int number_organisms, Seed & seed)
    : problem(problem),
      number_organisms(number_organisms),
      genes_per_organism(problem.get_number_parameters())
    {
        const int openmp_threads = omp_get_max_threads();
        for (int i = 0; i < openmp_threads; i++)
            random_generators.push_back(RandomGenerator(seed));

        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        if (number_organisms % mpi_size != 0) {
            //printf("Number of nodes should divide number of organisms completely\n");
            //printf("We will add more mutants and elites for that\n");
            int additional_organisms = mpi_size - number_organisms % mpi_size;
            number_organisms += additional_organisms;
        }

        min_gene.resize(genes_per_organism, nan(""));
        max_gene.resize(genes_per_organism, nan(""));
        is_mutation_applicable.resize(genes_per_organism, -1);

        int boundaries_status = problem.get_boundaries(min_gene, max_gene, is_mutation_applicable);

        if (boundaries_status == -1) {
            if (mpi_rank == 0)
                std::cerr << "non-constrained optimization problems are not supported"
                    << std::endl;
            throw 1;
        }

        best_global_genes.resize(genes_per_organism,  nan(""));
        all_best_genes_fitness_values.resize(number_organisms,  nan(""));

        int all_genes_size = number_organisms * genes_per_organism;
        all_genes.resize(all_genes_size,  nan(""));
        all_best_genes.resize(all_genes_size,  nan(""));

        //initial zero velocities
        velocities.resize(all_genes_size,  0);
        init_vector.resize(genes_per_organism);
    }

    template<typename InitializedRandomGenerator>
    void reset_organism(InitializedRandomGenerator & rg, std::vector<double>::iterator genes,
        std::vector<double>::iterator vel)
    {
        if (1) {
            for (int j = 0; j < genes_per_organism; j++) {
                if (is_mutation_applicable[j]) {
                    vel[j] = std::uniform_real_distribution<double>(- (max_gene[j] - min_gene[j])/2, (max_gene[j] - min_gene[j])/2)(rg);
                } else {
                    vel[j] = 0;
                }
                genes[j] = init_vector[j];
            }
        } else {
            for (int j = 0; j < genes_per_organism; j++) {
                vel[j] = 0;
                if (is_mutation_applicable[j]) {
                    genes[j] = std::uniform_real_distribution<double>(min_gene[j], max_gene[j])(rg);
                } else {
                    assert(init_status != -1);
                    //or maybe it is fine to put zero if init_status == -1
                    genes[j] = init_vector[j];
                }
            }
        }
    }

    template<typename InitializedRandomGenerator>
    void init(InitializedRandomGenerator rg)
    {
        init_status = problem.initial_guess_for_optimizer(init_vector.begin());

        //lets have at least one guy from initial guess, maybe it is not that bad
        //but it will be only here and not when some organism has to be reset
        for (int j = 0; j < genes_per_organism; j++) {
            all_best_genes[j] = all_genes[j] = init_vector[j];
        }
        //the rest starting from 1
        for (int i = 1; i < number_organisms; i++) {
            reset_organism(rg, all_genes.begin() +  i * genes_per_organism,
                velocities.begin() + i * genes_per_organism);

            for (int j = 0; j < genes_per_organism; j++)
                all_best_genes[j + i * genes_per_organism] = all_genes[j + i * genes_per_organism];
        }
        run_func(true);
        global_sync(true);
    }

    void global_sync(bool first_call = false)
    {
        //we only need to find the best among the whole swarm
        //first, find local smallest value
        auto smallest_local_it = std::min_element(all_best_genes_fitness_values.begin(),
                                    all_best_genes_fitness_values.begin() + number_organisms/mpi_size);

        struct {
            double value;
            int index;
        } in, out;
        in.value = *smallest_local_it;
        //std::cout << "V" <<  in.value << std::endl;
        in.index = mpi_rank;
        //now find process number and value
        MPI_Allreduce(
            &in,
            &out,
            1,
            MPI_DOUBLE_INT,
            MPI_MINLOC,
            MPI_COMM_WORLD);

        if (0) {
            //lets update best_genes each time to get steady state
            if (!first_call && best_global_genes_fitness_value < out.value) {
                // best_global_genes_fitness_value stays the same
                return;
            }
        }
        best_global_genes_fitness_value = out.value;
        const int broadcaster = out.index;
        if (mpi_rank == broadcaster) {
            const int pos = smallest_local_it - all_best_genes_fitness_values.begin();
            for (int j = 0; j < genes_per_organism; j++) {
                best_global_genes[j] = all_best_genes[pos * genes_per_organism + j];
                //std::cout << best_global_genes[j] << std::endl;
            }
        }
        //broadcaster broadcasts new best_global_genes
        MPI_Bcast(
            best_global_genes.data(),
            genes_per_organism,
            MPI_DOUBLE,
            broadcaster,
            MPI_COMM_WORLD);
    }

    void run_generation()
    {
        // Loewe et al, 2016
        const double phi1 = 2.05;
        const double phi2 = 2.05;
        const double phi = phi1 + phi2;
        const double chi = 2.0 / (phi - 2 + std::sqrt(phi * phi  - 4 * phi));
        std::uniform_real_distribution<double> Uphi1d(0, phi1), Uphi2d(0, phi2);

        #pragma omp parallel for
        for (int i = 0; i < number_organisms / mpi_size; i++) {
            RandomGenerator & rg = random_generators[omp_get_thread_num()];
            auto genes = all_genes.begin() + i * genes_per_organism;
            auto best_genes = all_best_genes.begin() + i * genes_per_organism;
            auto vel = velocities.begin() + i * genes_per_organism;
            //first, update velocities
            for (int j = 0; j < genes_per_organism; j++) {
                if (!is_mutation_applicable[j])
                    continue;
                const double Uphi1 = Uphi1d(rg);
                const double Uphi2 = Uphi2d(rg);
                vel[j] = chi * (vel[j] + Uphi1 * (best_genes[j] - genes[j])
                        + Uphi2 * (best_global_genes[j] - genes[j]));
            }
            //update genes
            for (int j = 0; j < genes_per_organism; j++) {
                if (!is_mutation_applicable[j])
                    continue;
                genes[j] += vel[j];
                //std::cout << genes[j] << std::endl;
            }
            //std::cout << std::endl;
        }
        run_func();
    }
    void run_func(bool first_call = false)
    {
        #pragma omp parallel for
        for (int i = 0; i < number_organisms / mpi_size; i++) {
            auto genes = all_genes.begin() + i * genes_per_organism;
            auto best_genes = all_best_genes.begin() + i * genes_per_organism;
            auto vel = velocities.begin() + i * genes_per_organism;
            //call optimized function
            const double new_val = problem.get_objective_value(genes);
            if (1) {
                //lets update best_genes each time to get steady state
                all_best_genes_fitness_values[i] = problem.get_objective_value(best_genes);
                if (all_best_genes_fitness_values[i] > 100) {
                    reset_organism(random_generators[omp_get_thread_num()], genes, vel);
                    for (int j = 0; j < genes_per_organism; j++)
                        best_genes[j] = genes[j];
                }
            }
            if (new_val > 100) { //hardcoded here, TODO
                //reset only genes and vel
                reset_organism(random_generators[omp_get_thread_num()], genes, vel);
            }
            //std::cout << "***FIT: " << new_val << std::endl;
            //update best parameter set for this organism
            if (first_call || new_val < all_best_genes_fitness_values[i]) {
                all_best_genes_fitness_values[i] = new_val;
                for (int j = 0; j < genes_per_organism; j++)
                    best_genes[j] = genes[j];
            }
        }
    }

    void fitness_function_values(std::vector<std::pair<double, int>> & sd_n_index)
    {
        return;
        //TODO
        for (int i = 0; i < number_organisms; i++) {
            sd_n_index[i].first = all_best_genes_fitness_values[i];
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

    std::vector<double> best() const
    {
        return best_global_genes;
    }

    std::vector<std::pair<int, double>> error_per_gen;
    std::vector<std::pair<int, double>> get_error_per_gen() const
    {
        return error_per_gen;
    }
    void collect_statistics()
    {
        //TODO
    }
    void log(const std::vector<std::pair<double, int>> & sd_n_index, int gen, int total_gen)
    {
        //TODO
        error_per_gen.push_back({gen, best_global_genes_fitness_value});
        if (mpi_rank != 0) return;
        if (gen % 1 != 0) return ;
        std::cout << "Generation: " << gen << std::endl
                  << "Best: " << best_global_genes_fitness_value << std::endl;
        /*
        error_per_gen.push_back({gen, sd_n_index[0].first});
        if (mpi_rank != 0) return;
        if (gen % 1 != 0) return ;
        std::cout << "Generation: " << gen << std::endl
                  << "Best: " << sd_n_index[0].first << std::endl
                  << "Worst: " << sd_n_index.back().first << std::endl;
        //uncomment me when collect_statistics done
        //problem.gen_algo_stats(sd_n_index, all_genes, gen, total_gen);
        */

       // auto best_genes = best();
      //  std::cout << "Genes:";
       // for (auto &g: best_genes)
       //     std::cout << " " << g;
       // std::cout << std::endl << std::endl;
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
        problem.submit_result(best_global_genes);
    }
};


template <typename Pop>
void particle_swarm_optimization(Pop & pop, const int generations)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Number of organisms: %d\n", pop.number_organisms);
        printf("Number of optimized parameters: %d\n", pop.get_genes_per_organism());
        printf("Number of generations: %d\n", generations);
        printf("Number of MPI nodes: %d\n", size);
        printf("Number of cores at root: %d\n", omp_get_max_threads());
    }
    //timer
    const double start_time = MPI_Wtime();

    //main PSO cycle
    for (int index_generation = 0; index_generation < generations; index_generation++) {
        double total_time = MPI_Wtime();

        double run_generation_time = MPI_Wtime();
        pop.run_generation();
        run_generation_time = MPI_Wtime() - run_generation_time;

        double global_sync_time = MPI_Wtime();
        pop.global_sync();
        global_sync_time = MPI_Wtime() - global_sync_time;

        // to dump log, collect all data at main node
        // the algorithm itself does not really need it
        pop.collect_statistics();

        if (rank == 0) {
            //store pairs of (error, index in state_struct) sorted by error in increasing order
            //thus, first elements are for elite organisms
            //TODO

            std::vector<std::pair<double, int>> sd_n_index(pop.number_organisms);
            /*
            pop.fitness_function_values(sd_n_index);

            double sort_time = MPI_Wtime();
            //sort by error increasing
            std::sort(sd_n_index.begin(), sd_n_index.end(),
                      [](const std::pair<double, int> &left_element, const std::pair<double, int> &right_element) {
                          return left_element.first < right_element.first;
                      });
            sort_time = MPI_Wtime() - sort_time;

            if (!(sd_n_index[0].first <= sd_n_index.back().first)) {
                std::cout << sd_n_index[0].first << " " << sd_n_index.back().first << std::endl;
            }
            assert(sd_n_index[0].first <= sd_n_index.back().first);
            */
            double log_time = MPI_Wtime();
            pop.log(sd_n_index, index_generation, generations);
            log_time = MPI_Wtime() - log_time;

            total_time = MPI_Wtime() - total_time;

            printf("\nGeneration %d\n", index_generation);
            printf("total_time         %9.3f %3d%%\n", total_time, 100);
            printf("run_gen_time       %9.3f %3d%%\n", run_generation_time, (int) (run_generation_time / total_time * 100));
            printf("global_sync_time   %9.3f %3d%%\n", global_sync_time, (int) (global_sync_time / total_time * 100));
            printf("log_time           %9.3f %3d%%\n", log_time, (int) (log_time / total_time * 100));
            //printf("sort_time          %9.3f %3d%%\n", sort_time, (int) (sort_time / total_time * 100));
        }
    }
    const double end_time = MPI_Wtime();
    if (rank == 0)
        printf("PSO time: %f sec\n", end_time - start_time);

    pop.done();
}

#endif
