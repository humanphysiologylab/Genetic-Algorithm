#ifndef GA_CAUCHY_MUTATION_H
#define GA_CAUCHY_MUTATION_H

#include <cmath>
#include <random>
#include <omp.h>
#include <vector>
#include <cassert>
#include <iostream>
#include "base_mutation.h"

template <typename RandomGenerator, typename Seed>
class CauchyMutation:
    public BaseMutation
{
    std::vector<RandomGenerator> random_generators;
    double mutrate; //probability of mutation
    double gamma;
    std::vector<double> v_gamma;
    int num_mutation_applicable_genes;
    
    void generate_uniform_vector_on_sphere(double * uniform_vector, int size)
    {
        RandomGenerator & rg = random_generators[omp_get_thread_num()];
           
        std::normal_distribution<double> nd(0, 1);
        for (int i = 0; i < size; i++)
            uniform_vector[i] = nd(rg);
        
        double len = 0;
        for (int i = 0; i < size; i++)
            len += std::pow(uniform_vector[i], 2);
        
        len = std::sqrt(len);
        for (int i = 0; i < size; i++)
            uniform_vector[i] /= len;
    }

    void transform_gene_forward(double & min, double & max, double & gene, double v_gamma)
    {
        /*
         * input: default min, max, gene, gamma
         * output: overwritten min, max, gene
         */
        max = (max - min) / v_gamma;
        gene = (gene - min) / v_gamma;
        min = 0;
    }
    void transform_gene_forward_log(double & min, double & max, double & gene, double v_gamma)
    {
        /*
         * input: default min, max, gene, gamma
         * output: overwritten min, max, gene
         */
        assert(min > 0);
        max = (std::log10(max) - std::log10(min)) / v_gamma;
        gene = (std::log10(gene) - std::log10(min)) / v_gamma;
        min = 0;
    }
    void transform_gene_backward(double min, double max, double & gene, double v_gamma)
    {
        /*
         * input: default min, max, gene to transform back, v_gamma
         * output: transformed back gene
         */
        gene = gene * v_gamma + min;
    }
    void transform_gene_backward_log(double min, double max, double & gene, double v_gamma)
    {
        /*
         * input: default min, max, gene to transform back, v_gamma
         * output: transformed back gene
         */
        gene = std::pow(10, gene * v_gamma + std::log10(min));
    }
    void cauchy_mutation(double * genes_input, double * genes_output,
                        const double * min_value,
                        const double * max_value,
                        int number_genes, const int * is_mutation_applicable)
    {
        /* Cauchy mutation
         * borders are mirrors, bounce-bounce
         * genes_input can be equal to genes_output
         */
        std::uniform_real_distribution<double> ran(0, 1);
        RandomGenerator & rg = random_generators[omp_get_thread_num()];

        const double mut_prob = ran(rg);

        if (mut_prob <= 1 - mutrate) {
            for (int i = 0; i < number_genes; i++)
                genes_output[i] = genes_input[i];
        } else {
            double uniform_vector[num_mutation_applicable_genes];
            generate_uniform_vector_on_sphere(uniform_vector, num_mutation_applicable_genes);

            int mut_applicable_iter = 0;

            const double shift = gamma * std::tan((M_PI / 2 * ran(rg)));
            for (int j = 0; j < number_genes; j++) {
                if (!is_mutation_applicable[j]) {
                     genes_output[j] = genes_input[j];
                     continue;
                }
                double min = min_value[j], max = max_value[j];
                double gene = std::max(min, std::min(max, genes_input[j]));
                if (is_mutation_applicable[j] == 1)
                    transform_gene_forward(min, max, gene, v_gamma[j]);
                else if (is_mutation_applicable[j] == 2)
                    transform_gene_forward_log(min, max, gene, v_gamma[j]);
                else
                    throw("Wrong transform type");
                const double L = max - min;
                const double a = gene - min;
                const double b = max - gene;

                double x = shift * uniform_vector[mut_applicable_iter];

                assert(L > 0);
                if (x >= 0)
                    x = std::fmod(x, 2 * L);
                else
                    x = 2 * L + std::fmod(x, 2 * L);

                const double y = std::abs(std::abs(x - b) - L) - a;
                gene += y;
                if (is_mutation_applicable[j] == 1)
                    transform_gene_backward(min_value[j], max_value[j], gene, v_gamma[j]);
                else if (is_mutation_applicable[j] == 2)
                    transform_gene_backward_log(min_value[j], max_value[j], gene, v_gamma[j]);

                genes_output[j] = gene;
                if (std::isnan(gene)) {
                    std::cout << "Gene " << j << "is nan after mutation" << std::endl;
                    throw("Bad mutation");
                }
                mut_applicable_iter++;
            }
        }
    }
public:
    CauchyMutation(Seed & seed, double mutrate = 0.1, double gamma = 1, std::vector<double> v_gamma = std::vector<double>())
    : mutrate(mutrate), gamma(gamma), v_gamma(v_gamma)
    {
        assert(v_gamma.size() != 0);
        const int openmp_threads = omp_get_max_threads();
        for (int i = 0; i < openmp_threads; i++)
            random_generators.push_back(RandomGenerator(seed));
        num_mutation_applicable_genes = static_cast<int>(v_gamma.size());
        for (double g: v_gamma) {
            if (g == 0)
                num_mutation_applicable_genes--;
        }
    }
    
    void operator()(double * population_genes, const double * min_value, const double * max_value, int population_size, int genes_number, const int * is_mutation_applicable)
    {
        #pragma omp parallel for
        for (int i = 0; i < population_size; i++) {
            cauchy_mutation(population_genes + i * genes_number,
                            population_genes + i * genes_number,
                            min_value,
                            max_value,
                            genes_number,
                            is_mutation_applicable);
        }
    }
};

#endif //GA_CAUCHY_MUTATION_H
