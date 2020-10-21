#include "cauchy_mutation.h"


const double pi = 3.141592653589793238462643;
const double MUTRATE = 0.9; //probability of mutation


void normalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                               int number_organisms, int number_genes, int number_logscale_params,
        /*out*/   double *genes_output, double *left_border_normalized, double *right_border_normalized,
        /*other*/ double left_border_value, double right_border_value) {


    for (int i = 0; i < number_organisms; ++i) {
        for (int j = 0; j < number_genes; ++j) {

            if (j < number_logscale_params) {
                genes_output[i * number_genes + j] = log10(genes_input[i * number_genes + j]);
                left_border_normalized[j] = log10(left_border[j]);
                right_border_normalized[j] = log10(right_border[j]);
            } else {
                genes_output[i * number_genes + j] = genes_input[i * number_genes + j];
                left_border_normalized[j] = left_border[j];
                right_border_normalized[j] = right_border[j];
            }

            genes_output[i * number_genes + j] = left_border_value +
                                                 (genes_output[i * number_genes + j] - left_border_normalized[j]) /
                                                 (right_border_normalized[j] - left_border_normalized[j]) *
                                                 (right_border_value - left_border_value);
            left_border_normalized[j] = left_border_value;
            right_border_normalized[j] = right_border_value;
            //std::cout << genes_input[i * number_genes + j] << " -- > " << genes_output[i * number_genes + j] << "\n";
        }
    }
}


void denormalize_genes(/*in*/    double *genes_input, double *left_border, double *right_border,
                                 int number_organisms, int number_genes, int number_logscale_params,
        /*out*/   double *genes_output,
        /*other*/ double left_border_value, double right_border_value) {


    for (int i = 0; i < number_organisms; ++i) {
        for (int j = 0; j < number_logscale_params; ++j) {

            genes_output[i * number_genes + j] = log10(left_border[j]) +
                                                 (genes_input[i * number_genes + j] - left_border_value) /
                                                 (right_border_value - left_border_value) *
                                                 (log10(right_border[j]) - log10(left_border[j]));

            genes_output[i * number_genes + j] = pow(10, genes_output[i * number_genes + j]);
            //std::cout << genes_input[i * number_genes + j] << " -- > " << genes_output[i * number_genes + j] << "\n";

        }
        for (int j = number_logscale_params; j < number_genes; ++j) {

            genes_output[i * number_genes + j] = left_border[j] +
                                                 (genes_input[i * number_genes + j] - left_border_value) /
                                                 (right_border_value - left_border_value) *
                                                 (right_border[j] - left_border[j]);

            //std::cout << genes_input[i * number_genes + j] << " -- > " << genes_output[i * number_genes + j] << "\n";
        }
    }
}


void transform_genes(/*in*/    double *genes_input, double *left_border_input, double *right_border_input,
                               int number_organisms, int number_genes, int number_log10scale_params,
                               const double *left_border_output, const double *right_border_output,
        /*out*/   double *genes_output) {

    double l = 0, r = 0;

    for (int i = 0; i < number_organisms; ++i) {
        for (int j = 0; j < number_genes; ++j) {

            if (j < number_log10scale_params) {
                genes_output[i * number_genes + j] = log10(genes_input[i * number_genes + j]);
                l = log10(left_border_input[j]);
                r = log10(right_border_input[j]);
            } else {
                genes_output[i * number_genes + j] = genes_input[i * number_genes + j];
                l = left_border_input[j];
                r = right_border_input[j];
            }

            genes_output[i * number_genes + j] = left_border_output[j] +
                                                 (genes_output[i * number_genes + j] - l) / (r - l) *
                                                 (right_border_output[j] - left_border_output[j]);

            /*
            std::cout << "[" << left_border_input[j] << "; " << genes_input[i * number_genes + j] << "; "
                      << right_border_input[j] << "] -- > [" << left_border_output[j] << "; "
                      << genes_output[i * number_genes + j] << "; " << right_border_output[j] << "]\n";
                      */

        }
    }
}


void transform_genes_back(/*in*/    double *genes_input, double *left_border_input, double *right_border_input,
                                    int number_organisms, int number_genes, int number_log10scale_params,
                                    const double *left_border_output, const double *right_border_output,
        /*out*/   double *genes_output) {

    for (int i = 0; i < number_organisms; ++i) {
        for (int j = 0; j < number_log10scale_params; ++j) {

            genes_output[i * number_genes + j] = log10(left_border_output[j]) +
                                                 (genes_input[i * number_genes + j] - left_border_input[j]) /
                                                 (right_border_input[j] - left_border_input[j]) *
                                                 (log10(right_border_output[j]) - log10(left_border_output[j]));

            genes_output[i * number_genes + j] = pow(10, genes_output[i * number_genes + j]);

            /*
            std::cout << "[" << left_border_input[j] << "; " << genes_input[i * number_genes + j] << "; "
                      << right_border_input[j] << "] -- > [" << left_border_output[j] << "; "
                      << genes_output[i * number_genes + j] << "; " << right_border_output[j] << "]\n";
                      */

        }
        for (int j = number_log10scale_params; j < number_genes; ++j) {

            genes_output[i * number_genes + j] = left_border_output[j] +
                                                 (genes_input[i * number_genes + j] - left_border_input[j]) /
                                                 (right_border_input[j] - left_border_input[j]) *
                                                 (right_border_output[j] - left_border_output[j]);

            /*
            std::cout << "[" << left_border_input[j] << "; " << genes_input[i * number_genes + j] << "; "
                      << right_border_input[j] << "] -- > [" << left_border_output[j] << "; "
                      << genes_output[i * number_genes + j] << "; " << right_border_output[j] << "]\n";
                      */

        }
    }
}


void box_muller_transform(double *uniform_vector, int NUMBER_GENES, long seed) {

    double u = 0, v = 0, s = 0;
    //double mu = 0, epsilon = 1.;
    double vector_len = 0.;

    for (int i = 0; i < NUMBER_GENES; i++) {

        do {
            u = ran2(&seed) * 2 - 1;
            v = ran2(&seed) * 2 - 1;
            s = u * u + v * v;
        } while ((s == 0) || (s >= 1));

        double sq = sqrt(-2 * log(s) / s);
        double z1 = u * sq;

        uniform_vector[i] = z1; // mu + z1 * epsilon;
        vector_len += uniform_vector[i] * uniform_vector[i];

    }

    vector_len = sqrt(vector_len);

    for (int i = 0; i < NUMBER_GENES; i++) {
        uniform_vector[i] = uniform_vector[i] / vector_len;
    }
}


void cauchy_mutation(double *genes_output, double *genes_input, double *left_border, double *right_border,
                     int NUMBER_ORGANISMS, int NUMBER_GENES) {
    // borders are mirrors, bounce-bounce

    struct timeval tv;
    gettimeofday(&tv, NULL);
    long seed = long(tv.tv_usec);
    long seed_negative = -seed;
    ran2(&seed_negative);

    double uniform_vector[NUMBER_GENES];

    double gamma = 1.; // mast be consistent with main()

    for (int i = 0; i < NUMBER_ORGANISMS; i++) {

        double mut_prob = ran2(&seed);

        if (mut_prob <= 1 - MUTRATE) {
            for (int j = 0; j < NUMBER_GENES; j++) {
                genes_output[i * NUMBER_GENES + j] = genes_input[i * NUMBER_GENES + j];
            }

        } else {

            box_muller_transform(uniform_vector, NUMBER_GENES, seed);
            double shift = gamma * tan((pi / 2 * ran2(&seed)));

            double a, b, x, y, L;
            for (int j = 0; j < NUMBER_GENES; j++) {
                L = right_border[j] - left_border[j];
                a = genes_input[i * NUMBER_GENES + j] - left_border[j];
                b = right_border[j] - genes_input[i * NUMBER_GENES + j];
                x = shift * uniform_vector[j];
                while (x < 0) x += 2 * L;
                x = fmod(x, 2 * L);
                y = fabs(fabs(x - b) - L) - a;
                genes_output[i * NUMBER_GENES + j] = genes_input[i * NUMBER_GENES + j] + y;
            }

        }
    }
}