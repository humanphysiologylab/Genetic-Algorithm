#include "gradient_descent.h"


double ema(double f, double ema_prev, double coef)
{
    return (1 - coef) * f + coef * ema_prev;
}
