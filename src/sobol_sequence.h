#ifndef SOBOL_SEQUENCE
#define SOBOL_SEQUENCE

#include <boost/random/sobol.hpp>

void scale_vector(std::vector<double> & to_scale, const std::vector<double> & min, const std::vector<double> & max);

class Sobol
{
    boost::random::sobol sobol_generator;
    double max_value, min_value;
    double scaler;
public:
    Sobol(unsigned int dim, double min = 0, double max = 1);
    std::vector<double> operator()();
};


#endif
