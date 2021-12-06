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

template <typename OptimizationProblem>
double fitn(OptimizationProblem & problem, std::vector<double> & params, const std::vector<double> & min_v, const std::vector<double> & max_v,
            const std::vector<double> & prior, double alpha_regularization)
{
    double penalty = 0;
    double reg = 0;
    for (int i = 0; i < params.size(); i++) {
        reg += std::pow(prior[i] - params[i], 2);
        if (min_v[i] > params[i])
            penalty += std::pow(20 * (params[i] - min_v[i]) / min_v[i], 2);
        if (max_v[i] < params[i])
            penalty += std::pow(20 * (params[i] - max_v[i]) / max_v[i], 2);
    }
    reg *= alpha_regularization;
    if (std::isnan(reg)) {
        for (int i = 0; i < params.size(); i++)
            std::cout << "Prior: " << prior[i]  << " param: " << params[i] << std::endl;
        throw("Bad parameter");
    }
    //std::cout << "Reg: " << reg << std::endl;
    //penalty = 0;////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return penalty + reg + problem.genetic_algorithm_calls(params.begin());
}

#endif
