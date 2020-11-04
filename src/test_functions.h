#ifndef TEST_FUNCTIONS
#define TEST_FUNCTIONS

#include "basic_function_functor.h"

class Basic1dFunctionFunctor:
    public BasicFunctionFunctor
{
public:
    virtual void operator()(double * x, double *y) const = 0;

    double operator()(const std::vector<double> & x) const;
};

class SphereFunction:
    public Basic1dFunctionFunctor
{
public:
    using Basic1dFunctionFunctor::operator ();
    SphereFunction(int xdim_);
    double ymin() const;
    virtual void operator()(double * x, double *y) const override;
    virtual std::vector<double> solution() const override;
};

class RosenbrockFunction:
    public Basic1dFunctionFunctor
{
public:
    using Basic1dFunctionFunctor::operator ();
    RosenbrockFunction(int xdim_);
    double ymin() const;
    virtual void operator()(double * x, double *y) const override;
    virtual std::vector<double> solution() const override;
};

class RastriginFunction:
    public Basic1dFunctionFunctor
{
public:
    using Basic1dFunctionFunctor::operator ();
    RastriginFunction(int xdim_);
    double ymin() const;
    virtual void operator()(double * x, double *y) const override;
    virtual std::vector<double> solution() const override;
};

class StyblinskiTangFunction:
    public Basic1dFunctionFunctor
{
public:
    using Basic1dFunctionFunctor::operator ();
    StyblinskiTangFunction(int xdim_);
    double ymin() const;
    virtual void operator()(double * x, double *y) const override;
    virtual std::vector<double> solution() const override;
};
#endif
