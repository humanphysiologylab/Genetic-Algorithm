#ifndef BASIC_POPULATION
#define BASIC_POPULATION

#include <mpi.h>
#include <omp.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <cassert>

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

    double *min_gene;
    double *max_gene;
    int    *is_mutation_applicable;

    double *all_genes;
    int genes_per_organism;
    int all_genes_size;

    double *elite_genes_buffer;
    int elite_genes_buffer_size;

    double *mutant_genes_buffer;
    int mutant_genes_buffer_size;

    double *fitness_values;

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


        min_gene = new double [genes_per_organism];
        max_gene = new double [genes_per_organism];
        is_mutation_applicable = new int [genes_per_organism];

        int boundaries_status = problem.get_boundaries(min_gene, max_gene, is_mutation_applicable);

        if (boundaries_status == -1) {
            std::cerr << "non-constrained optimization problems are not supported"
                << std::endl;
            throw;
        }

        all_genes_size = number_organisms * genes_per_organism;
        all_genes = new double [all_genes_size];

        elite_genes_buffer_size = number_elites * genes_per_organism;
        elite_genes_buffer = new double [elite_genes_buffer_size];

        mutant_genes_buffer_size = number_mutants * genes_per_organism;
        mutant_genes_buffer = new double [mutant_genes_buffer_size];

        for (int i = 0; i < elite_genes_buffer_size; i++)
            elite_genes_buffer[i] = nan("");
        for (int i = 0; i < mutant_genes_buffer_size; i++)
            mutant_genes_buffer[i] = nan("");

        fitness_values = new double [number_organisms];
    }

    ~BasicPopulation()
    {
        delete [] min_gene;
        delete [] max_gene;
        delete [] is_mutation_applicable;
        delete [] all_genes;
        delete [] elite_genes_buffer;
        delete [] mutant_genes_buffer;
        delete [] fitness_values;
    }

    template<typename InitializedRandomGenerator>
    void init(InitializedRandomGenerator rg)
    {
        std::vector<double> init_vector(genes_per_organism);
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
        MPI_Scatter(all_genes, all_genes_size / mpi_size,
                    MPI_DOUBLE, (mpi_rank == 0) ? MPI_IN_PLACE : all_genes,
                    all_genes_size / mpi_size,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void gather()
    {
        MPI_Gather((mpi_rank == 0) ? MPI_IN_PLACE : fitness_values, number_organisms / mpi_size, MPI_DOUBLE,
                   fitness_values, number_organisms / mpi_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gather((mpi_rank == 0) ? MPI_IN_PLACE : all_genes, all_genes_size / mpi_size,
                   MPI_DOUBLE, all_genes,
                   all_genes_size / mpi_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void run_generation()
    {
        #pragma omp parallel for
        for (int i = 0; i < number_organisms / mpi_size; i++) {
            fitness_values[i] = problem.genetic_algorithm_calls(all_genes + i * genes_per_organism);
        }
    }

    void fitness_function(std::vector<std::pair<double, int>> & sd_n_index)
    {
        if (mpi_rank != 0)
            return;
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
        copy_genes(all_genes + from * genes_per_organism, elite_genes_buffer + to * genes_per_organism);
    }
    void save_mutant_to_mutant_buffer(int from, int to)
    {
        copy_genes(all_genes + from * genes_per_organism, mutant_genes_buffer + to * genes_per_organism);
    }
    void restore_elites_to_main_array()
    {
        for (int i = 0; i < number_elites; i++)
            copy_genes(elite_genes_buffer + i * genes_per_organism, all_genes + i * genes_per_organism);
    }
    void restore_mutants_to_main_array()
    {
        for (int i = 0; i < number_mutants; i++)
            copy_genes(mutant_genes_buffer + i * genes_per_organism, all_genes + (number_elites + i) * genes_per_organism);
    }
    std::vector<double> best() const
    {
        return std::vector<double>(all_genes, all_genes + genes_per_organism);
    }
    void log(const std::vector<std::pair<double, int>> & sd_n_index, int gen)
    {
        if (gen % 1 != 0) return ;
        std::cout << "Generation: " << gen << std::endl
                  << "Value: " << sd_n_index[0].first << std::endl;

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
        return min_gene;
    }
    double * get_max_gene_value()
    {
        return max_gene;
    }
    int get_number_genes()
    {
        return genes_per_organism;
    }
    void done()
    {
        problem.genetic_algorithm_result(all_genes);
    }
};

#endif
