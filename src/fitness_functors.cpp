#include "fitness_functors.h"

double ScalarFunctionMinFitnessFunctor::operator()(double * x, double * y) const
{
    /* goal is to minimize scalar function
     * 
     */
    return y[0];
}
