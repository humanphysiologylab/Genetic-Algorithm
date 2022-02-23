#include "sobol_sequence.h"
Sobol::Sobol(unsigned int dim, double min, double max)
: sobol_generator(dim), max_value(max), min_value(min)
{
    scaler = (max_value - min_value) / (sobol_generator.max() - sobol_generator.min());
}
std::vector<double> Sobol::operator()()
{
    std::vector<double> res(sobol_generator.dimension());
    for (int i = 0; i < sobol_generator.dimension(); i++)
        res[i] = (sobol_generator() - sobol_generator.min()) * scaler + min_value;
    return res;
}


void scale_vector(std::vector<double> & to_scale, const std::vector<double> & min, const std::vector<double> & max)
{
    assert(to_scale.size() == min.size());
    assert(min.size() == max.size());
    for (int i = 0; i < to_scale.size(); i++)
        to_scale[i] = (max[i] - min[i]) * to_scale[i] + min[i];
}

