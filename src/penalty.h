#ifndef PENALTY
#define PENALTY

#include <vector>
#include <cmath>

template <typename OptimizationProblem>
double fitn(OptimizationProblem & problem, std::vector<double> & params, const std::vector<int> & is_mutation_applicable,
                const std::vector<double> & min_v, const std::vector<double> & max_v)
{
    double penalty = 0;
    for (int i = 0; i < is_mutation_applicable.size(); i++) {
        if (!is_mutation_applicable[i]) continue;
 
        if (min_v[i] > params[i])
            penalty += std::pow(20 * (params[i] - min_v[i]) / min_v[i], 2);
        if (max_v[i] < params[i])
            penalty += std::pow(20 * (params[i] - max_v[i]) / max_v[i], 2);
    }
    //penalty = 0;////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return penalty + problem.genetic_algorithm_calls(params.begin());
}



template <typename OptimizationProblem>
double fitn(OptimizationProblem & problem, std::vector<double> & params, const std::vector<double> & min_v, const std::vector<double> & max_v)
{
    double penalty = 0;
    for (int i = 0; i < params.size(); i++) {
        if (min_v[i] > params[i])
            penalty += std::pow(20 * (params[i] - min_v[i]) / min_v[i], 2);
        if (max_v[i] < params[i])
            penalty += std::pow(20 * (params[i] - max_v[i]) / max_v[i], 2);
    }
    //penalty = 0;////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return penalty + problem.genetic_algorithm_calls(params.begin());
}

#endif
