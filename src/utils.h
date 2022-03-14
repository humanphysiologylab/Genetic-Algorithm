#ifndef UTILS_H
#define UTILS_H
#include <vector>

double ema(double f, double ema_prev, double coef);

std::vector<double> diff(const std::vector<double> & x);
std::vector<double> ema(const std::vector<double> & x, double coef);

#endif
