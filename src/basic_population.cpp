#include "basic_population.h"
#include <cassert>

SphereFunction::SphereFunction(int xdim_)
{
    xdim = xdim_;
    ydim = 1;

    min_x = std::vector<double>(xdim, -5);
    max_x = std::vector<double>(xdim, 5);
}

void SphereFunction::operator()(double * x, double *y)
{
    double res = 0;
    for (int i = 0; i < xdim; i++)
        res += x[i] * x[i];
    y[0] = res - 3.14;
}

std::vector<double> SphereFunction::solution()
{
    return std::vector<double>(xdim, 0);
}

RosenbrockFunction::RosenbrockFunction(int xdim_)
{
    assert(xdim_ > 1);
    xdim = xdim_;
    ydim = 1;

    min_x = std::vector<double>(xdim, -5);
    max_x = std::vector<double>(xdim, 5);
}

void RosenbrockFunction::operator()(double * x, double *y)
{
    double res = 0;
    for (int i = 0; i < xdim - 1; i++)
        res += 100 * pow(x[i + 1] - x[i] * x[i], 2) + pow(1 - x[i], 2);
    y[0] = res;
}

std::vector<double> RosenbrockFunction::solution()
{
    return std::vector<double>(xdim, 1);
}

RastriginFunction::RastriginFunction(int xdim_)
{
    xdim = xdim_;
    ydim = 1;

    min_x = std::vector<double>(xdim, -5.12);
    max_x = std::vector<double>(xdim, 5.12);
}

void RastriginFunction::operator()(double * x, double *y)
{
    const double A = 10;
    double res = A * xdim;
    for (int i = 0; i < xdim; i++)
        res += x[i] * x[i] - A * cos(2 * M_PI * x[i]);
    y[0] = res;
}

std::vector<double> RastriginFunction::solution()
{
    return std::vector<double>(xdim, 0);
}


double ScalarFunctionMinFitnessFunctor::operator()(double * x, double *y) const
{
    /* goal is to minimize scalar function
     * 
     */
    return y[0];
}
