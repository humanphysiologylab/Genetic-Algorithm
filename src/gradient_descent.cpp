#include "gradient_descent.h"

BaseCoreGD::BaseCoreGD(int max_steps)
    : max_steps(max_steps)
{
}

CoreSimpleGradientDescent::CoreSimpleGradientDescent(int max_steps, double r_eps, double learning_rate)
    : BaseCoreGD(max_steps), r_eps(r_eps), learning_rate(learning_rate)
{
}
std::vector<double> CoreSimpleGradientDescent::operator()(const std::vector<double> & initial_guess, const std::vector<int> & is_mutation_applicable, std::function<double(std::vector<double> &)> func)
{
    //store vector and loss value on each step
    std::vector<double> res;
    int param_num = initial_guess.size();
    res.reserve((param_num + 1) * max_steps);

    std::vector<int> mut_pos;
    for (int i = 0 ; i < param_num; i++) {
        if (is_mutation_applicable[i])
            mut_pos.push_back(i);
    }
    std::vector<double> sol = initial_guess;

    std::vector<double> df(mut_pos.size()), dfdx(mut_pos.size());
    for (int step = -1; step < max_steps; step++) {
        const double f = func(sol);
        if (step == -1) //warmup
            continue;
        //std::cout << step << " f: " << f << std::endl;
        res.insert(res.end(), sol.begin(), sol.end());
        res.insert(res.end(), f);
        #pragma omp parallel for
        for (size_t i = 0; i < mut_pos.size(); i++) {

            std::vector<double> params = sol;
            const double eps = r_eps;///@todo * std::max(std::abs(params[mut_pos[i]]), 1e-12);
            params[mut_pos[i]] += eps;
            df[i] = func(params) - f;
            dfdx[i] = df[i] / eps;

            /*
            std::vector<double> params_fwd = sol, params_back = sol;
            double eps = r_eps * std::max(std::abs(sol[mut_pos[i]]), 1e-12); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = problem.get_objective_value(params_fwd.begin()) -
                    problem.get_objective_value(params_back.begin());
            dfdx[i] = df[i] / (2*eps);
            */
        }

        for (size_t i = 0; i < mut_pos.size(); i++) {
            const int pos = mut_pos[i];
            double x_new = sol[pos] - learning_rate * dfdx[i];
            if (!std::isnan(x_new)) {
                sol[pos] = x_new;
            }
        }
    }
    return res;
}

CoreAdam::CoreAdam(int max_steps, double r_eps, double learning_rate, double beta1, double beta2)
    : BaseCoreGD(max_steps), r_eps(r_eps), learning_rate(learning_rate), beta1(beta1), beta2(beta2)
{
}

std::vector<double> CoreAdam::operator()(const std::vector<double> & initial_guess, const std::vector<int> & is_mutation_applicable, std::function<double(std::vector<double> &)> func)
{
    //store vector and loss value on each step
    std::vector<double> res;
    int param_num = initial_guess.size();
    res.reserve((param_num + 1) * max_steps);

    std::vector<int> mut_pos;
    for (int i = 0 ; i < param_num; i++) {
        if (is_mutation_applicable[i])
            mut_pos.push_back(i);
    }
    std::vector<double> sol = initial_guess;

    std::vector<double> df(mut_pos.size()), dfdx(mut_pos.size());
    std::vector<double> ema_sq_grad(mut_pos.size()); // squared gradient exponential moving average
    std::vector<double> ema_grad(mut_pos.size()); // ema for gradient
    for (int step = -1; step < max_steps; step++) {
        //const double f = fitn(problem, sol, is_mutation_applicable, min_v, max_v);
        const double f = func(sol);
        if (step == -1) //warmup
            continue;
        //std::cout << step << " f: " << f << std::endl;
        res.insert(res.end(), sol.begin(), sol.end());
        res.insert(res.end(), f);
        #pragma omp parallel for
        for (size_t i = 0; i < mut_pos.size(); i++) {

            std::vector<double> params = sol;
            const double eps = r_eps;///@todo * std::max(std::abs(params[mut_pos[i]]), 1e-12);
            params[mut_pos[i]] += eps;
            df[i] = func(params) - f;
            dfdx[i] = df[i] / eps;

            /*
            std::vector<double> params_fwd = sol, params_back = sol;
            double eps = r_eps * std::max(std::abs(sol[mut_pos[i]]), 1e-12); //r_eps 1e-3 for rosenbrock, 1e-1 for ap models
            params_fwd[mut_pos[i]] += eps;
            params_back[mut_pos[i]] -= eps;
            df[i] = problem.get_objective_value(params_fwd.begin()) -
                    problem.get_objective_value(params_back.begin());
            dfdx[i] = df[i] / (2*eps);
            */
        }

        for (size_t i = 0; i < mut_pos.size(); i++) {
            const int pos = mut_pos[i];
            ema_grad[i] = ema(dfdx[i], ema_grad[i], beta1);
            ema_sq_grad[i] = ema(std::pow(dfdx[i], 2), ema_sq_grad[i], beta2);

            double corrected_ema_grad = ema_grad[i] / (1 - std::pow(beta1, step + 1));
            double corrected_ema_sq_grad = ema_sq_grad[i] / (1 - std::pow(beta2, step + 1));

            double x_new = sol[pos] - learning_rate * corrected_ema_grad / (std::sqrt(corrected_ema_sq_grad) + 1e-12);
            /*
            if (min_v[pos] > x_new || max_v[pos] < x_new) {
                std::cout << "out of boundary" << std::endl;
                std::cout << "sol[" << pos << "]: " << sol[pos] << std::endl;
                std::cout << "dfdx[" << i << "]: " << dfdx[i] << std::endl;
                std::cout << "ss[" << pos << "]: " << x_new << std::endl;
                std::cout << "min max: " << min_v[pos] << " " << max_v[pos] <<  std::endl;
            }
            */
            if (!std::isnan(x_new)) {
                sol[pos] = x_new;
            }
        }
    }

    return res;
}
