#include "test_functions.h"
#include <cassert>
#include <cmath>


double Basic1dFunctionFunctor::operator()(const std::vector<double> & x) const
{
    double y;
    std::vector<double> xx = x;
    operator()(xx.data(), &y);
    return y;
}

SphereFunction::SphereFunction(int xdim_)
{
    xdim = xdim_;
    ydim = 1;

    min_x = std::vector<double>(xdim, -5);
    max_x = std::vector<double>(xdim, 5);
}

void SphereFunction::operator()(double * x, double *y) const
{
    double res = 0;
    for (int i = 0; i < xdim; i++)
        res += x[i] * x[i];
    y[0] = res;
}

std::vector<double> SphereFunction::solution() const
{
    return std::vector<double>(xdim, 0);
}

double SphereFunction::ymin() const
{
    return 0;
}

RosenbrockFunction::RosenbrockFunction(int xdim_)
{
    assert(xdim_ > 1);
    xdim = xdim_;
    ydim = 1;

    min_x = std::vector<double>(xdim, -5);
    max_x = std::vector<double>(xdim, 5);
}

void RosenbrockFunction::operator()(double * x, double *y) const
{
    double res = 0;
    for (int i = 0; i < xdim - 1; i++)
        res += 100 * pow(x[i + 1] - x[i] * x[i], 2) + pow(1 - x[i], 2);
    y[0] = res;
}

std::vector<double> RosenbrockFunction::solution() const
{
    return std::vector<double>(xdim, 1);
}
double RosenbrockFunction::ymin() const
{
    return 0;
}

RastriginFunction::RastriginFunction(int xdim_)
{
    xdim = xdim_;
    ydim = 1;

    min_x = std::vector<double>(xdim, -5.12);
    max_x = std::vector<double>(xdim, 5.12);
}

void RastriginFunction::operator()(double * x, double *y) const
{
    const double A = 10;
    double res = A * xdim;
    for (int i = 0; i < xdim; i++)
        res += x[i] * x[i] - A * cos(2 * M_PI * x[i]);
    y[0] = res;
}

std::vector<double> RastriginFunction::solution() const
{
    return std::vector<double>(xdim, 0);
}

double RastriginFunction::ymin() const
{
    return 0;
}


StyblinskiTangFunction::StyblinskiTangFunction(int xdim_)
{
    xdim = xdim_;
    ydim = 1;

    min_x = std::vector<double>(xdim, -5);
    max_x = std::vector<double>(xdim, 5);
}

void StyblinskiTangFunction::operator()(double * x, double *y) const
{
    double res = 0;
    for (int i = 0; i < xdim; i++)
        res += pow(x[i], 4) - 16 * pow(x[i], 2) + 5 * x[i];
    y[0] = res / 2;
}

std::vector<double> StyblinskiTangFunction::solution() const
{
    return std::vector<double>(xdim, -2.903534);
}

double StyblinskiTangFunction::ymin() const
{
    return -39.16599 * xdim;
}
