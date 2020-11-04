#ifndef GENETIC_ALGORITHM
#define GENETIC_ALGORITHM

#include <mpi.h>
#include <iostream>
#include <algorithm>

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
            fitness_time = MPI_Wtime() - fitness_time;

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
            pop.log(sd_n_index, index_generation);
            log_time = MPI_Wtime() - log_time;
            
            
            /*Genetic Operators for mutants*/

            int mpool[pop.number_mutants]; //mpool: mutant_index -> next_generation_index
            //so let's find it
            double selection_time = MPI_Wtime();
            selection(mpool, pop.number_mutants, sd_n_index, pop.number_elites);
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
                          pop.get_number_genes());
            crossover_time = MPI_Wtime() - crossover_time;


            double mutation_time = MPI_Wtime();
            mutation(pop.get_mutant_buffer_genes(), pop.get_min_gene_value(), pop.get_max_gene_value(), pop.number_mutants, pop.get_number_genes());
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
}

#endif
