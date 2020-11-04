#include "basic_function_functor.h"
#include <cassert>
#include <cmath>

double BasicFunctionFunctor::solution_error_l2(const std::vector<double> & s) const
{
    double error = 0;
    auto exact_sol = solution();
    assert(exact_sol.size() == s.size());
    for (unsigned i = 0; i < s.size(); i++)
        error += pow(s[i] - exact_sol[i], 2);
    return sqrt(error / exact_sol.size());
}
