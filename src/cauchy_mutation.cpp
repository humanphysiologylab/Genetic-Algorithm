#include <sys/time.h>
#include <cstdlib>
#include <cmath>
#include "cauchy_mutation.h"


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



