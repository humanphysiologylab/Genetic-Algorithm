#ifndef GENETIC_ALGORITHM
#define GENETIC_ALGORITHM

#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <omp.h>
#include <vector>
#include <cmath>


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
            throw 1;
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
        int init_status = problem.initial_guess_for_optimizer(init_vector.begin());
        
        //lets have at least one guy from initial guess, maybe it is not that bad
        for (int j = 0; j < genes_per_organism; j++) {
            all_genes[j] = init_vector[j];
        }
        //the rest starting from 1
        for (int i = 1; i < number_organisms; i++) {
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
    void log(const std::vector<std::pair<double, int>> & sd_n_index, int gen, int total_gen)
    {
        error_per_gen.push_back({gen, sd_n_index[0].first});
        if (mpi_rank != 0) return;
        if (gen % 1 != 0) return ;
        std::cout << "Generation: " << gen << std::endl
                  << "Best: " << sd_n_index[0].first << std::endl
                  << "Worst: " << sd_n_index.back().first << std::endl;
        problem.gen_algo_stats(sd_n_index, all_genes, gen, total_gen);
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



template <typename Pop, typename Selection, typename Crossover, typename Mutation>
void genetic_algorithm(Pop & pop, Selection  selection, Crossover  crossover, Mutation  mutation, const int generations)
{
    /* Initialize population before calling genetic algorithm
     * 
     * 
     * 
     */

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Number of organisms: %d\n", pop.number_organisms);
        printf("Number of elite organisms: %d\n", pop.number_elites);
        printf("Number of optimized parameters: %d\n", pop.get_genes_per_organism());
        printf("Number of generations: %d\n", generations);
        printf("Number of MPI nodes: %d\n", size);
        printf("Number of cores at root: %d\n", omp_get_max_threads());
    }

    //timer
    const double start_time = MPI_Wtime();


    //main GA cycle
    for (int index_generation = 0; index_generation < generations; index_generation++) {
        double total_time = MPI_Wtime();
        double scatter_time = total_time;
        pop.scatter();
        scatter_time = MPI_Wtime() - scatter_time;

        double run_generation_time = MPI_Wtime();
        pop.run_generation();
        run_generation_time = MPI_Wtime() - run_generation_time;

        double gather_time = MPI_Wtime();
        pop.gather();
        gather_time = MPI_Wtime() - gather_time;


        //store pairs of (error, index in state_struct) sorted by error in increasing order
        //thus, first elements are for elite organisms
        std::vector<std::pair<double, int>> sd_n_index(pop.number_organisms);
         
        double fitness_time = MPI_Wtime();
        pop.fitness_function(sd_n_index);
        fitness_time = MPI_Wtime() - fitness_time;


        if (rank == 0) {

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

            double save_elite_time = MPI_Wtime();
            //save elites to elite buffer and later copy it to main array
            for (int i = 0; i < pop.number_elites; i++) {
                const int elite_index = sd_n_index[i].second;
                pop.save_elite_to_elite_buffer(elite_index, i); //from -> to
            }
            save_elite_time = MPI_Wtime() - save_elite_time;
            
            double log_time = MPI_Wtime();
            pop.log(sd_n_index, index_generation, generations);
            log_time = MPI_Wtime() - log_time;
            
            
            /*Genetic Operators for mutants*/

            int mpool[pop.number_mutants]; //mpool: mutant_index -> next_generation_index
            //so let's find it
            double selection_time = MPI_Wtime();
            selection(mpool, pop.number_mutants, sd_n_index, pop.number_elites, (index_generation < 100 ? 8: 2 ) * sd_n_index[0].first);
            selection_time = MPI_Wtime() - selection_time;
            //Do not use sd_n_index anymore!

            double save_mutant_time = MPI_Wtime();
            //only now copy states from next_generation to buf_mutant_state according to the shuffled mpool!
            for (int i = 0; i < pop.number_mutants; i++) {
                const int mutant_index = mpool[i];
                pop.save_mutant_to_mutant_buffer(mutant_index, i); //from -> to
            }
            save_mutant_time = MPI_Wtime() - save_mutant_time;
            //no need of mpool anymore


            double crossover_time = MPI_Wtime();
            crossover(pop.get_mutant_buffer_genes(), pop.get_min_gene_value(), pop.get_max_gene_value(), pop.number_mutants,
                          pop.get_number_genes(), pop.get_is_mutation_applicable());
            crossover_time = MPI_Wtime() - crossover_time;


            double mutation_time = MPI_Wtime();
            mutation(pop.get_mutant_buffer_genes(), pop.get_min_gene_value(), pop.get_max_gene_value(), pop.number_mutants, pop.get_number_genes(), pop.get_is_mutation_applicable());
            mutation_time = MPI_Wtime() - mutation_time;

            /*Genetic Operators for mutants DONE*/
            double restore_organisms_time = MPI_Wtime();
            //now copy elite to main arrays
            pop.restore_elites_to_main_array();

            //and copy mutants to main arrays
            pop.restore_mutants_to_main_array();
            restore_organisms_time = MPI_Wtime() - restore_organisms_time;

            total_time = MPI_Wtime() - total_time;

            printf("\nGeneration %d\n", index_generation);
            printf("total_time         %9.3f %3d%%\n", total_time, 100);
            printf("scatter_time       %9.3f %3d%%\n", scatter_time, (int) (scatter_time / total_time * 100));
            printf("run_gen_time       %9.3f %3d%%\n", run_generation_time, (int) (run_generation_time / total_time * 100));
            printf("gather_time        %9.3f %3d%%\n", gather_time, (int) (gather_time / total_time * 100));
            printf("fitness_time       %9.3f %3d%%\n", fitness_time, (int) (fitness_time / total_time * 100));
            printf("selection_time     %9.3f %3d%%\n", selection_time, (int) (selection_time / total_time * 100));
            printf("crossover_time     %9.3f %3d%%\n", crossover_time, (int) (crossover_time / total_time * 100));
            printf("mutation_time      %9.3f %3d%%\n", mutation_time, (int) (mutation_time / total_time * 100));
            printf("log_time           %9.3f %3d%%\n", log_time, (int) (log_time / total_time * 100));
            printf("restore_org_time   %9.3f %3d%%\n", restore_organisms_time, (int) (restore_organisms_time / total_time * 100));
            printf("save_mutant_time   %9.3f %3d%%\n", save_mutant_time, (int) (save_mutant_time / total_time * 100));
            printf("save_elite_time    %9.3f %3d%%\n", save_elite_time, (int) (save_elite_time / total_time * 100));
            printf("sort_time          %9.3f %3d%%\n", sort_time, (int) (sort_time / total_time * 100));

        }
    }
    const double end_time = MPI_Wtime();
    if (rank == 0)
        printf("GA time: %f sec\n", end_time - start_time);

    pop.done();
}

#endif
