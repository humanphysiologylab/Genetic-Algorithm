#ifndef BASIC_POPULATION
#define BASIC_POPULATION

#include <mpi.h>
#include <omp.h>
#include <vector>
#include <cmath>
#include <iostream>
#include "random_number_generator.h"





template <typename FunctionFunctor, typename FitnessFunctor>
class BasicPopulation
{

public:
    FunctionFunctor function;
    FitnessFunctor fitness;

    int number_organisms;
    int number_elites;
    int number_mutants;

    int mpi_rank;
    int mpi_size;

    double *max_gene;
    double *min_gene;

    double *all_genes;
    int genes_per_organism;
    int all_genes_size;

    double *elite_genes_buffer;
    int elite_genes_buffer_size;

    double *mutant_genes_buffer;
    int mutant_genes_buffer_size;

    double *all_y;
    int y_per_organism;
    int all_y_size;

    BasicPopulation(FunctionFunctor function, FitnessFunctor fitness, int number_elites_p, int number_mutants_p)
    : function(function),
      fitness(fitness),
      number_organisms(number_elites_p + number_mutants_p),
      number_elites(number_elites_p),
      number_mutants(number_mutants_p),
      genes_per_organism(function.xdim),
      y_per_organism(function.ydim)
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
        
        max_gene = new double [genes_per_organism];
        min_gene = new double [genes_per_organism];

        for (int i = 0; i < genes_per_organism; i++) {
            max_gene[i] = function.max_x[i];
            min_gene[i] = function.min_x[i];
        }
        
        all_genes_size = number_organisms * genes_per_organism;
        all_genes = new double [all_genes_size];
        
        elite_genes_buffer_size = number_elites * genes_per_organism;
        elite_genes_buffer = new double [elite_genes_buffer_size];
        
        mutant_genes_buffer_size = number_mutants * genes_per_organism;
        mutant_genes_buffer = new double [mutant_genes_buffer_size];
        
        all_y_size = number_organisms * y_per_organism;
        all_y = new double [all_y_size];

        for (int i = 0; i < elite_genes_buffer_size; i++)
            elite_genes_buffer[i] = nan("");
        for (int i = 0; i < mutant_genes_buffer_size; i++)
            mutant_genes_buffer[i] = nan("");
        for (int i = 0; i < all_y_size; i++)
            all_y[i] = nan("");
    }

    ~BasicPopulation()
    {
        delete [] max_gene;
        delete [] min_gene;
        delete [] all_genes;
        delete [] elite_genes_buffer;
        delete [] mutant_genes_buffer;
        delete [] all_y;
    }

    void init()
    {
        /* Basic initialization:
         * all genes are random
         * 
         */
        long seed = (long) time(NULL);
        long seed_negative = -seed;
        ran2(&seed_negative);

        for (int i = 0; i < number_organisms; i++)
            for (int j = 0; j < genes_per_organism; j++)
            all_genes[j + i * genes_per_organism] = min_gene[j] + (max_gene[j] - min_gene[j]) * ran2(&seed);

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
        MPI_Gather((mpi_rank == 0) ? MPI_IN_PLACE : all_y, all_y_size / mpi_size, MPI_DOUBLE,
                   all_y, all_y_size / mpi_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gather((mpi_rank == 0) ? MPI_IN_PLACE : all_genes, all_genes_size / mpi_size,
                   MPI_DOUBLE, all_genes,
                   all_genes_size / mpi_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void run_generation()
    {
        #pragma omp parallel for
        for (int i = 0; i < number_organisms / mpi_size; i++) {
            function(all_genes + i * genes_per_organism, all_y + i * y_per_organism);
        }
    }

    void fitness_function(std::vector<std::pair<double, int>> & sd_n_index)
    {
        for (int i = 0; i < number_organisms; i++) {
            sd_n_index[i].first = fitness(all_genes + i * genes_per_organism, all_y + i * y_per_organism);
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
    
    void log(std::vector<std::pair<double, int>> & sd_n_index, int gen)
    {
        if (gen % 10 != 0) return ;
        std::cout << "Generation: " << gen << std::endl
                  << "Value: " << sd_n_index[0].first << std::endl;
        auto best_genes = best();
        std::cout << "Genes:";
        for (auto &g: best_genes)
            std::cout << " " << g;
        std::cout << std::endl << std::endl;
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
    
};

#endif
