#include "objective.h"

std::vector<double> diff(const std::vector<double> & x)
{
    std::vector<double> res(x.size() - 1);
    for (size_t i = 0; i < res.size(); i++)
        res[i] = x[i + 1] - x[i];
    return res;
}

std::vector<double> ema(const std::vector<double> & x, double coef)
{
    std::vector<double> res(x.size());
    res[0] = x[0];
    for (size_t i = 1; i < res.size(); i++)
        res[i] = coef * x[i] + (1 - coef) * res[i - 1];
    return res;
}
