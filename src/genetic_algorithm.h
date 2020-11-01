#ifndef GENETIC_ALGORITHM
#define GENETIC_ALGORITHM

#include <mpi.h>
#include <iostream>
#include <algorithm>

#include "polynomial_mutation.h"
#include "sbx_crossover.h"
#include "tournament_selection.h"


template <typename Pop>
void genetic_algorithm(Pop & pop, const int generations)
{
    /* Initialize population before calling genetic algorithm
     * 
     * 
     * 
     */

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

        if (rank == 0) {
            //store pairs of (error, index in state_struct) sorted by error in increasing order
            //thus, first elements are for elite organisms
            std::vector<std::pair<double, int>> sd_n_index(pop.number_organisms);

            //TODO fitness_function should be evaluated at each node
            //And there is not need to send AP to the root

            double fitness_time = MPI_Wtime();

            pop.fitness_function(sd_n_index);

            //now sort by error increasing
            std::sort(sd_n_index.begin(), sd_n_index.end(),
                      [](const std::pair<double, int> &left_element, const std::pair<double, int> &right_element) {
                          return left_element.first < right_element.first;
                      });

            fitness_time = MPI_Wtime() - fitness_time;

            if (!(sd_n_index[0].first <= sd_n_index.back().first)) {
                std::cout << sd_n_index[0].first << " " << sd_n_index.back().first << std::endl;
            }
            assert(sd_n_index[0].first <= sd_n_index.back().first);

            /*Genetic Operators for mutants*/

            int mpool[pop.number_mutants]; //mpool: mutant_index -> next_generation_index
            //so let's find it
            double tournament_time = MPI_Wtime();
            tournament_selection(mpool, pop.number_mutants, sd_n_index, pop.number_elites);

            //we also need to shuffle mpool because it is sorted!
            for (int i = 0; i < pop.number_mutants; i++) {
                int j = i + rand() % (pop.number_mutants - i);
                std::swap(mpool[i], mpool[j]);
            }
            tournament_time = MPI_Wtime() - tournament_time;


            //save elites to elite buffer and later copy it to main array
            for (int i = 0; i < pop.number_elites; i++) {
                const int elite_index = sd_n_index[i].second;
                pop.save_elite_to_elite_buffer(elite_index, i); //from -> to
            }
            //only now copy states from next_generation to buf_mutant_state according to the shuffled mpool!
            for (int i = 0; i < pop.number_mutants; i++) {
                const int mutant_index = mpool[i];
                pop.save_mutant_to_mutant_buffer(mutant_index, i); //from -> to
            }
            //no need of mpool anymore

            double log_time = MPI_Wtime();
            pop.log(sd_n_index, index_generation);
            log_time = MPI_Wtime() - log_time;


            double crossover_time = MPI_Wtime();
            sbx_crossover(pop.get_mutant_buffer_genes(), pop.get_min_gene_value(), pop.get_max_gene_value(), pop.number_mutants,
                          pop.get_number_genes());
            crossover_time = MPI_Wtime() - crossover_time;


            double mutation_time = MPI_Wtime();
            polynomial_mutation(pop.get_mutant_buffer_genes(), pop.get_min_gene_value(), pop.get_max_gene_value(), pop.number_mutants, pop.get_number_genes());
            mutation_time = MPI_Wtime() - mutation_time;

            /*Genetic Operators for mutants DONE*/

            //now copy elite to main arrays
            pop.restore_elites_to_main_array();

            //and copy mutants to main arrays
            pop.restore_mutants_to_main_array();


            total_time = MPI_Wtime() - total_time;
/*
            printf("\nGeneration %d\n", index_generation);
            printf("total_time         %9.3f %3d%%\n", total_time, 100);
            printf("scatter_time       %9.3f %3d%%\n", scatter_time, (int) (scatter_time / total_time * 100));
            printf("run_gen_time       %9.3f %3d%%\n", run_generation_time, (int) (run_generation_time / total_time * 100));
            printf("gather_time        %9.3f %3d%%\n", gather_time, (int) (gather_time / total_time * 100));
            printf("fitness_time       %9.3f %3d%%\n", fitness_time, (int) (fitness_time / total_time * 100));
            printf("tournament_time    %9.3f %3d%%\n", tournament_time, (int) (tournament_time / total_time * 100));
            printf("crossover_time     %9.3f %3d%%\n", crossover_time, (int) (crossover_time / total_time * 100));
            printf("mutation_time      %9.3f %3d%%\n", mutation_time, (int) (mutation_time / total_time * 100));
            printf("log_time           %9.3f %3d%%\n", log_time, (int) (log_time / total_time * 100));
*/
        }
    }

    const double end_time = MPI_Wtime();
    if (rank == 0)
        printf("GA time: %f sec\n", end_time - start_time);
}

#endif
