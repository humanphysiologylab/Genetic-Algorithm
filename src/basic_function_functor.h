#ifndef BASIC_FUNCTION_FUNCTOR
#define BASIC_FUNCTION_FUNCTOR

#include <array>
#include <cmath>


template <int xdim, int ydim = 1, typename numtype = double>
class BasicFunctionFunctor
{
public:
    static_assert(xdim >= 0, "function should have arguments");
    static_assert(ydim >= 1, "function should have an output");
    using argument = std::array<numtype, xdim>;
    using value = std::array<numtype, ydim>;
    using num = numtype;
    virtual argument x_min() const = 0; 
    virtual argument x_max() const = 0;
    virtual value operator()(const argument & x) const = 0;
    virtual argument solution() const = 0;
    int get_xdim() const
    {
        return xdim;
    }
    int get_ydim() const
    {
        return ydim;
    }
    template<int power = 2>
    numtype solution_error(const argument & inexact_sol) const
    {
        numtype error = 0;
        const auto exact_sol = solution();
        for (int i = 0; i < xdim; i++)
            error += pow(std::abs(inexact_sol[i] - exact_sol[i]), power);
        return pow(error / xdim, 1.0/power);
    }
};
#endif
