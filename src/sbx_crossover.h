#ifndef GA_SBX_CROSSOVER_H
#define GA_SBX_CROSSOVER_H

#include <cmath>
#include <cassert>
#include <algorithm>
#include <random>

template<typename InitializedRandomGenerator>
class SBXcrossover
{
    InitializedRandomGenerator rg;
    const int crossrate = 0.9; //probability of crossover
    const int eta_c = 10; //The order of the polynomial for the SBX crossover

public:
    SBXcrossover(InitializedRandomGenerator rg)
    : rg(rg)
    {}
    
    void operator()(double *next_generation, double *left_border, double *right_border,
                   int number_organisms, int number_genes)
    {
        /*
        adopted realcross from NSGA-II: Non-dominated Sorting Genetic Algorithm - II

        Authors: Dr. Kalyanmoy Deb, Sameer Agrawal, Amrit Pratap, T Meyarivan
        Paper Title: A Fast and Elitist multi-objective Genetic Algorithm: NSGA-II
        Journal: IEEE Transactions on Evolutionary Computation (IEEE-TEC)
        Year: 2002
        Volume: 6
        Number: 2
        Pages: 182-197
        */

        std::uniform_real_distribution<double> ran(0, 1);
        for (int i = 0; i < number_organisms / 2; i++) {
            if (ran(rg) <= crossrate) {
                for (int j = 0; j < number_genes; j++) {
                    if (ran(rg) <= 0.5) {
                        double y1 = next_generation[2 * i * number_genes + j],
                               y2 = next_generation[(2 * i + 1) * number_genes + j];

                        if (fabs(y1 - y2) > 1e-14) {

                            if (y1 > y2)
                                std::swap(y1, y2);

                            const double yl = left_border[j];
                            const double yu = right_border[j];
                            const double rand = ran(rg);
                            double beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                            double alpha = 2.0 - pow(beta,-(eta_c+1.0));
                            double betaq;
                            if (rand <= (1.0/alpha))
                            {
                                betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                            }
                            else
                            {
                                betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                            }
                            double c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                            beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                            alpha = 2.0 - pow(beta,-(eta_c+1.0));
                            if (rand <= (1.0/alpha))
                            {
                                betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                            }
                            else
                            {
                                betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                            }
                            double c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                            if (c1<yl)
                                c1=yl;
                            if (c2<yl)
                                c2=yl;
                            if (c1>yu)
                                c1=yu;
                            if (c2>yu)
                                c2=yu;

                            if (ran(rg) <= 0.5) {
                                next_generation[2 * i * number_genes + j] = c2;
                                next_generation[(2 * i + 1) * number_genes + j] = c1;
                            } else {
                                next_generation[2 * i * number_genes + j] = c1;
                                next_generation[(2 * i + 1) * number_genes + j] = c2;
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif //GA_SBX_CROSSOVER_H
