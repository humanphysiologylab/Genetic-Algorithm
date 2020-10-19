#include <sys/time.h>
#include <cstdlib>
#include <math.h>
#include "cauchy_mutation.h"
#include "random_number_generator.h"

const double pi = 3.141592653589793238462643;
const double MUTRATE = 0.9; //probability of mutation


void normalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                               int number_organisms, int number_genes, int number_conductancies,
                     /*out*/   double *genes_output, double *left_border_normalized, double *right_border_normalized,
                     /*other*/ double left_border_value, double right_border_value) {


    for (int i = 0; i < number_organisms; ++i) {
        for (int j = 0; j < number_genes; ++j) {

            if (j < number_conductancies) {
                genes_output[i * number_genes + j] = log10(genes_input[i * number_genes + j]);
                left_border_normalized[j] = log10(left_border[j]);
                right_border_normalized[j] = log10(right_border[j]);
            } else {
                genes_output[i * number_genes + j] = genes_input[i * number_genes + j];
                left_border_normalized[j] = left_border[j]; right_border_normalized[j] = right_border[j];
            }

            genes_output[i * number_genes + j] = left_border_value +
                                                 (genes_output[i * number_genes + j] - left_border_normalized[j]) /
                                                 (right_border_normalized[j] - left_border_normalized[j]) *
                                                 (right_border_value - left_border_value);
            left_border_normalized[j] = left_border_value; right_border_normalized[j] = right_border_value;
            //std::cout << genes_input[i * number_genes + j] << " -- > " << genes_output[i * number_genes + j] << "\n";
        }
    }
}

void denormalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                                 int number_organisms, int number_genes, int number_conductancies,
                       /*out*/   double *genes_output,
                       /*other*/ double left_border_value, double right_border_value) {


    for (int i = 0; i < number_organisms; ++i) {
        for (int j = 0; j < number_conductancies; ++j) {

            genes_output[i * number_genes + j] = log10(left_border[j]) +
                                                 (genes_input[i * number_genes + j] - left_border_value) /
                                                 (right_border_value - left_border_value) *
                                                 (log10(right_border[j]) - log10(left_border[j]));

            genes_output[i * number_genes + j] = pow(10, genes_output[i * number_genes + j]);
            //std::cout << genes_input[i * number_genes + j] << " -- > " << genes_output[i * number_genes + j] << "\n";

        }
        for (int j = number_conductancies; j < number_genes; ++j) {

            genes_output[i * number_genes + j] = left_border[j] +
                                                 (genes_input[i * number_genes + j] - left_border_value) /
                                                 (right_border_value - left_border_value) *
                                                 (right_border[j] - left_border[j]);

            //std::cout << genes_input[i * number_genes + j] << " -- > " << genes_output[i * number_genes + j] << "\n";
        }
    }
}


void box_muller_transform(double *uniform_vector, int NUMBER_GENES, long seed) {
    double u, v, s, z1, sq, rand_ud;
    double mu, epsilon;
    double vector_len = 0.;
    int item, zi;

    mu = 0.;
    epsilon = 1.;

    for (zi = 0; zi < NUMBER_GENES; zi++) {
        do {
            rand_ud = ran2(&seed);
            if (rand_ud > 0.5) {
                u = ran2(&seed);
            } else {
                u = -ran2(&seed);
            }

            rand_ud = ran2(&seed);
            if (rand_ud > 0.5) {
                v = ran2(&seed);
            } else {
                v = -ran2(&seed);
            }

            s = sqrt(u * u + v * v);
        } while ((s == 0) || (s >= 1));

        sq = sqrt(-2 * log(s) / s);

        z1 = u * sq;

        uniform_vector[zi] = mu + z1 * epsilon;

    }

    for (item = 0; item < NUMBER_GENES; item++) {
        vector_len += uniform_vector[item] * uniform_vector[item];
    }
    vector_len = sqrt(vector_len);

    for (item = 0; item < NUMBER_GENES; item++) {
        uniform_vector[item] = uniform_vector[item] / vector_len;
    }
}

void
cauchy_mutation(double *input, double *output, double *left_border, double *right_border, int NUMBER_ORGANISMS,
                int NUMBER_GENES) {

    double gamma = 0.05;
    int mu = 0;

    struct timeval tv;
    gettimeofday(&tv, NULL);
    long seed = long(tv.tv_usec);
    long seed_negative = -seed;
    ran2(&seed_negative);

    for (int i = 0; i < NUMBER_ORGANISMS; i++) {
        double mut_prob = ran2(&seed);

        if (mut_prob <= 1 - MUTRATE) {
            for (int j = 0; j < NUMBER_GENES; j++) {
                output[i * NUMBER_GENES + j] = input[i * NUMBER_GENES + j];
            }

        } else {

            double distance_min1 = 1000000.;
            double distance_min2 = 1000000.;
            double distance = 0;

            double uniform_vector[NUMBER_GENES];
            box_muller_transform(uniform_vector, NUMBER_GENES, seed);

            /* right border */
            for (int j = 0; j < NUMBER_GENES; j++) {
                if (uniform_vector[j] > 0) {
                    distance = (right_border[j] - input[i * NUMBER_GENES + j]) / uniform_vector[j];
                } else {
                    distance = (left_border[j] - input[i * NUMBER_GENES + j]) / uniform_vector[j];
                }
                if (distance < distance_min1) {
                    distance_min1 = distance;
                }
            }

            /* left border*/
            for (int j = 0; j < NUMBER_GENES; j++) {
                if (uniform_vector[j] > 0) {
                    distance = (input[i * NUMBER_GENES + j] - left_border[j]) / uniform_vector[j];
                } else {
                    distance = (input[i * NUMBER_GENES + j] - right_border[j]) / uniform_vector[j];
                }
                if (distance < distance_min2) {
                    distance_min2 = distance;
                }
            }


            double r_min = (atan(-distance_min2 / gamma)) / pi + 0.5;
            double r_max = (atan(distance_min1 / gamma)) / pi + 0.5;

            double rnd = r_min + (r_max - r_min) * ran2(&seed);

            //std::cout << "r_min = " << r_min << "; r_max = " << r_max << "; rnd = " << rnd << ";\n";
            double shift = mu + gamma * tan(pi * (rnd - 1. / 2));

            for (int j = 0; j < NUMBER_GENES; j++) {
                output[i * NUMBER_GENES + j] =
                        input[i * NUMBER_GENES + j] + shift * uniform_vector[j];

            }
        }

    }
}



