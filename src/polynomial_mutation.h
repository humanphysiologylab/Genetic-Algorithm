#ifndef POLYNOMIAL_MUTATION
#define POLYNOMIAL_MUTATION

#include <cmath>
#include <random>

template <typename InitializedRandomGenerator>
class PolynomialMutation
{
    InitializedRandomGenerator rg;
    const int eta_m = 20;//i dk is it good
    const double pmut_real = 0.1;//i dk is it good
    
    void real_mutate_ind (double *genes, double * min_realvar, double * max_realvar, int genes_number)
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

        for (int j = 0; j < genes_number; j++) {
            if (ran(rg) <= pmut_real) {
                double y = genes[j];
                double  deltaq;
                const double yl = min_realvar[j];
                const double yu = max_realvar[j];
                const double delta1 = (y-yl)/(yu-yl);
                const double delta2 = (yu-y)/(yu-yl);
                const double rnd = ran(rg);
                const double mut_pow = 1.0/(eta_m+1.0);
                if (rnd <= 0.5) {
                    const double xy = 1.0-delta1;
                    const double val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
                    deltaq =  pow(val,mut_pow) - 1.0;
                } else {
                    const double xy = 1.0-delta2;
                    const double val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
                    deltaq = 1.0 - (pow(val,mut_pow));
                }
                y = y + deltaq*(yu-yl);
                if (y<yl)
                    y = yl;
                if (y>yu)
                    y = yu;
                genes[j] = y;
            }
        }
    }
public:
    PolynomialMutation(InitializedRandomGenerator rg)
    : rg(rg)
    {}

    void operator()(double *population_genes, double * min_value, double * max_value, int population_size, int genes_number)
    {
        for (int i = 0; i < population_size; i++) {
            real_mutate_ind (population_genes, min_value, max_value, genes_number);
            population_genes += genes_number;
        }   
    }
    
};
#endif



