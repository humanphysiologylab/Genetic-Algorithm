#ifndef TEST_FUNCTIONS
#define TEST_FUNCTIONS

#include "basic_function_functor.h"

class SphereFunction:
    public BasicFunctionFunctor
{
public:
    
    SphereFunction(int xdim_);

    virtual void operator()(double * x, double *y) override;
    virtual std::vector<double> solution() override;
};

class RosenbrockFunction:
    public BasicFunctionFunctor
{
public:
    
    RosenbrockFunction(int xdim_);

    virtual void operator()(double * x, double *y) override;
    virtual std::vector<double> solution() override;
};

class RastriginFunction:
    public BasicFunctionFunctor
{
public:
    
    RastriginFunction(int xdim_);

    virtual void operator()(double * x, double *y) override;
    virtual std::vector<double> solution() override;
};

#endif
