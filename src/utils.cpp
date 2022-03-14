#include "utils.h"

double ema(double f, double ema_prev, double coef)
{
    return (1 - coef) * f + coef * ema_prev;
}

std::vector<double> diff(const std::vector<double> & x)
{
    std::vector<double> res(x.size() - 1);
    for (unsigned i = 0; i < res.size(); i++)
        res[i] = x[i + 1] - x[i];
    return res;
}

std::vector<double> ema(const std::vector<double> & x, double coef)
{
    std::vector<double> res(x.size());
    res[0] = x[0];
    for (unsigned i = 1; i < res.size(); i++)
        res[i] = ema(x[i], res[i - 1], coef);
    return res;
}
