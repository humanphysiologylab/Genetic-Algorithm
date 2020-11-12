#ifndef GA_CAUCHY_MUTATION_H
#define GA_CAUCHY_MUTATION_H

#include <cmath>
#include <random>
#include <omp.h>
#include <vector>
#include <cassert>

template <typename RandomGenerator, typename Seed>
class CauchyMutation
{
    std::vector<RandomGenerator> random_generators;
    double mutrate; //probability of mutation
    double gamma;
    
    void box_muller_transform(double * uniform_vector, int number_genes)
    {
        std::uniform_real_distribution<double> ran(0, 1);
        RandomGenerator & rg = random_generators[omp_get_thread_num()];
       
        double vector_len = 0;
        for (int i = 0; i < number_genes; i++) {
            double s = 0;
             //double mu = 0, epsilon = 1.;
            double u = 0;
            do {
                u = 2 * ran(rg) - 1;
                const double v = ran(rg);
                s = u * u + v * v;
            } while ((s == 0) || (s >= 1));

            const double sq = std::sqrt( - 2 * std::log(s) / s);
            const double z1 = u * sq;

            uniform_vector[i] = z1; // mu + z1 * epsilon;
            vector_len += uniform_vector[i] * uniform_vector[i];
        }
        vector_len = std::sqrt(vector_len);
        for (int i = 0; i < number_genes; i++)
            uniform_vector[i] /= vector_len;
    }
    
    void cauchy_mutation(double * genes_input, double * genes_output,
                        const double * min_value,
                        const double * max_value,
                        int number_genes)
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
            double uniform_vector[number_genes];
            box_muller_transform(uniform_vector, number_genes);
            
            const double shift = gamma * tan((M_PI / 2 * ran(rg)));
            for (int j = 0; j < number_genes; j++) {
                const double L = max_value[j] - min_value[j];
                const double a = genes_input[j] - min_value[j];
                const double b = max_value[j] - genes_input[j];


                double x = shift * uniform_vector[j];

                assert(L > 0);
                if (x >= 0)
                    x = std::fmod(x, 2 * L);
                else
                    x = 2 * L + std::fmod(x, 2 * L);

                const double y = std::abs(std::abs(x - b) - L) - a;
                genes_output[j] = genes_input[j] + y;
            }
        }
    }
public:
    CauchyMutation(Seed & seed, double mutrate = 0.1, double gamma = 1)
    : mutrate(mutrate), gamma(gamma)
    {
        const int openmp_threads = omp_get_max_threads();
        for (int i = 0; i < openmp_threads; i++)
            random_generators.push_back(RandomGenerator(seed));
    }
    
    void operator()(double *population_genes, const double * min_value, const double * max_value, int population_size, int genes_number)
    {
        #pragma omp parallel for
        for (int i = 0; i < population_size; i++) {
            cauchy_mutation(population_genes + i * genes_number,
                            population_genes + i * genes_number,
                            min_value,
                            max_value,
                            genes_number);
        }
    }
};



void normalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                               int number_organisms, int number_genes, int number_logscale_params,
        /*out*/   double *genes_output, double *left_border_normalized, double *right_border_normalized,
        /*other*/ double left_border_value = 0, double right_border_value = 1);

void denormalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                                 int number_organisms, int number_genes, int number_logscale_params,
        /*out*/   double *genes_output,
        /*other*/ double left_border_value = 0, double right_border_value = 1);

void transform_genes(/*in*/    double *genes_input, double *left_border_input, double *right_border_input,
                               int number_organisms, int number_genes, int number_log10scale_params,
                               const double *left_border_output, const double *right_border_output,
        /*out*/   double *genes_output);

void transform_genes_back(/*in*/    double *genes_input, double *left_border_input, double *right_border_input,
                                    int number_organisms, int number_genes, int number_log10scale_params,
                                    const double *left_border_output, const double *right_border_output,
        /*out*/   double *genes_output);

#endif //GA_CAUCHY_MUTATION_H
