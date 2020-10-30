#include <ctime>
#include <math.h>
#include <cassert>
#include <algorithm>
#include "sbx_crossover.h"
#include "random_number_generator.h"

const double CROSSRATE = 0.9; //probability of crossover

/*
void sbx_crossover(double *next_generation, double *after_cross, int *mpool, double *left_border, double *right_border,
                   int NUMBER_ORGANISMS, int NUMBER_GENES)
{
    // mpool should be already random!

    const int etaC = 10; //The order of the polynomial for the SBX crossover
    long seed, seed_negative;

    seed = (long) time(NULL);
    seed_negative = -seed;
    ran2(&seed_negative);

    for (int i = 0; i < NUMBER_ORGANISMS / 2; i++) {
        const int num_1 = mpool[2 * i];
        const int num_2 = mpool[2 * i + 1];
        const double cr_prob = ran2(&seed);

        if (cr_prob <= CROSSRATE) {
            for (int ii = 0; ii < NUMBER_GENES; ii++) {
                const double sw_prob = ran2(&seed);


                if (sw_prob < 0.5) // "NO"
                {
                    after_cross[2 * i * NUMBER_GENES + ii] = next_generation[num_1 * NUMBER_GENES + ii];
                    after_cross[(2 * i + 1) * NUMBER_GENES + ii] = next_generation[num_2 * NUMBER_GENES + ii];
                } else //"YES"
                {
                    
                    double p1 = next_generation[num_1 * NUMBER_GENES + ii],
                           p2 = next_generation[num_2 * NUMBER_GENES + ii];
                    if (p1 > p2)
                        std::swap(p1, p2);
                    
                    assert(p1 <= p2);
                    // Calculate beta
                    double beta;
                    if (right_border[ii] - p2 < p1 - left_border[ii]) {
                        beta = (p2 - p1) / (2.0 * right_border[ii] - p2 - p1);
                    } else {
                        beta = (p2 - p1) / (p2 + p1 - 2.0 * left_border[ii]);
                    }
                  
                    const double alpha = (2.0 - pow(beta, (double) (etaC + 1))) * ran2(&seed);
                    double beta_q;
                    if (alpha <= 1.0) {
                        beta_q = pow(alpha, 1.0 / ((double) (etaC + 1)));
                    } else {
                        beta_q = pow(1.0 / (2.0 - alpha), 1.0 / ((double) (etaC + 1)));
                    }
                    const double c1 = 0.5 * ((p1 + p2) + beta_q * (p1 - p2));
                    const double c2 = 0.5 * ((p1 + p2) + beta_q * (p2 - p1));

                    after_cross[2 * i * NUMBER_GENES + ii] = c1;
                    after_cross[(2 * i + 1) * NUMBER_GENES + ii] = c2;
                }

            }
        } else {
            for (int jj = 0; jj < NUMBER_GENES; jj++) {
                after_cross[2 * i * NUMBER_GENES + jj] = next_generation[num_1 * NUMBER_GENES + jj];
                after_cross[(2 * i + 1) * NUMBER_GENES + jj] = next_generation[num_2 * NUMBER_GENES + jj];
            }
        }
    }
}

*/

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

void sbx_crossover(double *next_generation, double *left_border, double *right_border,
                   int NUMBER_ORGANISMS, int NUMBER_GENES)
{
    const int eta_c = 10; //The order of the polynomial for the SBX crossover
    long seed, seed_negative;

    seed = (long) time(NULL);
    seed_negative = -seed;
    ran2(&seed_negative);
    
    for (int i = 0; i < NUMBER_ORGANISMS / 2; i++) {
        if (ran2(&seed) <= CROSSRATE) {
            for (int j = 0; j < NUMBER_GENES; j++) {
                if (ran2(&seed) <= 0.5) {
                    double y1 = next_generation[2 * i * NUMBER_GENES + j],
                           y2 = next_generation[(2 * i + 1) * NUMBER_GENES + j];

                    if (fabs(y1 - y2) > 1e-14) {

                        if (y1 > y2)
                            std::swap(y1, y2);

                        const double yl = left_border[j];
                        const double yu = right_border[j];
                        const double rand = ran2(&seed);
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

                        if (ran2(&seed)<=0.5) {
                            next_generation[2 * i * NUMBER_GENES + j] = c2;
                            next_generation[(2 * i + 1) * NUMBER_GENES + j] = c1;
                        } else {
                            next_generation[2 * i * NUMBER_GENES + j] = c1;
                            next_generation[(2 * i + 1) * NUMBER_GENES + j] = c2;
                        }
                    }
                }
            }
        }
    }
}
