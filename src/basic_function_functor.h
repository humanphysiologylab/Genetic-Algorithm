#ifndef BASIC_FUNCTION_FUNCTOR
#define BASIC_FUNCTION_FUNCTOR

#include <vector>

class BasicFunctionFunctor
{
public:
    int xdim, ydim;
    std::vector<double> min_x, max_x;
    virtual void operator()(double * x, double *y) = 0;
    virtual std::vector<double> solution() = 0;
};


#endif
