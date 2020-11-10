#ifndef TEST_FUNCTIONS
#define TEST_FUNCTIONS

#include "basic_function_functor.h"
#include <cassert>
#include <cmath>

template <int xdim, typename numtype = double>
class SphereFunction:
    public BasicFunctionFunctor<xdim, 1, numtype>
{
public:
    using Base = BasicFunctionFunctor<xdim, 1, numtype>;
    using typename Base::argument;
    using typename Base::value;
    numtype ymin() const
    {
        return 0;
    }
    virtual argument x_min() const override
    {
        argument r;
        r.fill(-5);
        return r;
    }
    virtual argument x_max() const override
    {
        argument r;
        r.fill(5);
        return r;
    }
    virtual value operator()(const argument & x) const override
    {
        numtype res = 0;
        for (int i = 0; i < xdim; i++)
            res += x[i]*x[i];
        return {res};
    }
    virtual argument solution() const override
    {
        return argument({});
    }
};



template <int xdim, typename numtype = double>
class RosenbrockFunction:
    public BasicFunctionFunctor<xdim, 1, numtype>
{
public:
    using Base = BasicFunctionFunctor<xdim, 1, numtype>;
    using typename Base::argument;
    using typename Base::value;
    numtype ymin() const
    {
        return 0;
    }
    virtual argument x_min() const override
    {
        argument r;
        r.fill(-5);
        return r;
    }
    virtual argument x_max() const override
    {
        argument r;
        r.fill(5);
        return r;
    }
    virtual value operator()(const argument & x) const override
    {
        numtype res = 0;
        for (int i = 0; i < xdim - 1; i++)
            res += 100 * std::pow(x[i + 1] - x[i] * x[i], 2) + std::pow(1 - x[i], 2);
        return {res};
    }
    virtual argument solution() const override
    {
        argument sol;
        sol.fill(1);
        return sol;
    }
};


template <int xdim, typename numtype = double>
class RastriginFunction:
    public BasicFunctionFunctor<xdim, 1, numtype>
{
public:
    using Base = BasicFunctionFunctor<xdim, 1, numtype>;
    using typename Base::argument;
    using typename Base::value;
    numtype ymin() const
    {
        return 0;
    }
    virtual argument x_min() const override
    {
        argument r;
        r.fill(-5.12);
        return r;
    }
    virtual argument x_max() const override
    {
        argument r;
        r.fill(5.12);
        return r;
    }
    virtual value operator()(const argument & x) const override
    {
        const numtype A = 10;
        numtype res = A * xdim;
        for (int i = 0; i < xdim; i++)
            res += x[i] * x[i] - A * std::cos(2 * M_PI * x[i]);
        return {res};
    }
    virtual argument solution() const override
    {
        return argument({});
    }
};


template <int xdim, typename numtype = double>
class StyblinskiTangFunction:
    public BasicFunctionFunctor<xdim, 1, numtype>
{
public:
    using Base = BasicFunctionFunctor<xdim, 1, numtype>;
    using typename Base::argument;
    using typename Base::value;
    numtype ymin() const
    {
        return -39.16599 * xdim;
    }
    virtual argument x_min() const override
    {
        argument r;
        r.fill(-5);
        return r;
    }
    virtual argument x_max() const override
    {
        argument r;
        r.fill(5);
        return r;
    }
    virtual value operator()(const argument & x) const override
    {
        numtype res = 0;
        for (int i = 0; i < xdim; i++)
            res += std::pow(x[i], 4) - 16 * std::pow(x[i], 2) + 5 * x[i];
        res /= 2;
        return {res};
    }
    virtual argument solution() const override
    {
        argument sol;
        sol.fill(-2.903534);
        return sol;
    }
};
#endif
