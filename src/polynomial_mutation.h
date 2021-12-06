#ifndef POLYNOMIAL_MUTATION
#define POLYNOMIAL_MUTATION

#include <cmath>
#include <random>
#include <omp.h>
#include <vector>
#include "base_mutation.h"

template <typename RandomGenerator, typename Seed>
class PolynomialMutation:
    public BaseMutation
{
    std::vector<RandomGenerator> random_generators;
    const int eta_m;
    const double pmut_real;

    void real_mutate_ind(double * genes, const double * min_realvar, const double * max_realvar, int genes_number, const int *is_mutation_applicable)
    {
        /* Routine for real polynomial mutation of an individual
         *
         * adopted real_mutate_ind from NSGA-II: Non-dominated Sorting Genetic Algorithm - II
         *
         * Authors: Dr. Kalyanmoy Deb, Sameer Agrawal, Amrit Pratap, T Meyarivan
         * Paper Title: A Fast and Elitist multi-objective Genetic Algorithm: NSGA-II
         * Journal: IEEE Transactions on Evolutionary Computation (IEEE-TEC)
         * Year: 2002
         * Volume: 6
         * Number: 2
         * Pages: 182-197
         */

        std::uniform_real_distribution<double> ran(0, 1);

        RandomGenerator & rg = random_generators[omp_get_thread_num()];

        for (int j = 0; j < genes_number; j++) {
            if (!is_mutation_applicable[j]) continue;

            if (ran(rg) <= pmut_real) {
                double y = genes[j];
                double deltaq;
                const double yl = min_realvar[j];
                const double yu = max_realvar[j];
                const double rnd = ran(rg);
                const double mut_pow = 1.0 / (eta_m + 1.0);
                if (rnd <= 0.5) {
                    const double delta1 = (y - yl) / (yu - yl);
                    const double xy = 1.0 - delta1;
                    const double val = 2.0 * rnd + (1.0 - 2.0 * rnd) * std::pow(xy, (eta_m + 1.0));
                    deltaq = std::pow(val, mut_pow) - 1.0;
                } else {
                    const double delta2 = (yu - y) / (yu - yl);
                    const double xy = 1.0 - delta2;
                    const double val = 2.0 * (1.0 - rnd) + (2.0 * rnd - 1.0) * std::pow(xy, (eta_m + 1.0));
                    deltaq = 1.0 - std::pow(val, mut_pow);
                }
                y = y + deltaq * (yu - yl);
                if (y < yl)
                    y = yl;
                else if (y > yu)
                    y = yu;
                genes[j] = y;
            }
        }
    }
public:
    PolynomialMutation(Seed & seed, double mutrate = 0.1, int eta_m = 20)
    : eta_m(eta_m), pmut_real(mutrate)
    {
        const int openmp_threads = omp_get_max_threads();
        for (int i = 0; i < openmp_threads; i++)
            random_generators.push_back(RandomGenerator(seed));

    }

    void operator()(double *population_genes, const double * min_value, const double * max_value, int population_size, int genes_number, const int * is_mutation_applicable)
    {
        #pragma omp parallel for
        for (int i = 0; i < population_size; i++) {
            real_mutate_ind(population_genes + i * genes_number, min_value, max_value, genes_number, is_mutation_applicable);
        }
    }

};


class NoMutation:
    public BaseMutation
{
public:
    void operator()(double *population_genes, const double * min_value, const double * max_value, int population_size, int genes_number, const int * is_mutation_applicable);
};
#endif
