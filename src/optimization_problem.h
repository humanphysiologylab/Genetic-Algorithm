#ifndef OPTIMIZATION_PROBLEM
#define OPTIMIZATION_PROBLEM

#include <string>
#include <vector>
#include <cassert>


template<int power = 2>
class MinimizeAPbaselines
{
public:
    using Baseline = std::vector<double>;
    using Period = double;
    using BaselinePeriod = std::pair<Baseline, Period>;
    using BaselinePeriodVec = std::vector<BaselinePeriod>;

    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        double res = 0;
        assert(a.size() == b.size());
        for (int i = 0; i < a.size(); i++)
            res += std::pow(std::abs(a[i] - b[i]), power);
        res /= a.size();
        return res;
    }

    double dist(const BaselinePeriodVec & a, const BaselinePeriodVec & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (int i = 0; i < a.size(); i++) {
            assert(a[i].second == b[i].second);
            res += distBaselines(a[i].first, b[i].first);
        }
        return std::pow(res, 1.0/power);
    }
};

template <typename Model, typename Solver, typename Objective>
class ODEoptimization
{
    Model & model;
    Solver & solver;
    Objective & obj;
    using Baseline = typename Objective::Baseline;
    using Period = typename Objective::Period;
    using BaselinePeriod = typename Objective::BaselinePeriod;
    using BaselinePeriodVec = typename Objective::BaselinePeriodVec;
    
    BaselinePeriodVec baselines;
    
    int number_parameters = 0;
public:
    ODEoptimization(Model & model, Solver & solver, Objective & obj)
    : model(model), solver(solver), obj(obj)
    {
        //TODO
    }

    int get_number_parameters() const
    {
        return number_parameters;
    }
    void read_config()
    {
        //TODO
        number_parameters = 100;
    }
    template <typename It>
    int get_boundaries(It pmin, It pmax)
    {
        //TODO
      //  auto ppmin = func.x_min();
     //   std::copy(ppmin.begin(), ppmin.end(), pmin);
     //   auto ppmax = func.x_max();
     //   std::copy(ppmax.begin(), ppmax.end(), pmax);
     //   return 0; //have boundaries
        return -1; //no boundaries
    }
    template <typename It>
    int initial_guess(It begin)
    {
        if (0) {
            //read from a file
            //TODO
            return 0;
        } else if (1) {
            //or get from model
            //TODO
            return 0;
        } else {
            //or just not provide anything
            return -1;//initial parameters are not provided
        }
    }

    template <typename It>
    double genetic_algorithm_calls(It parameters_begin)
    {
        //unpack parameters
        //TODO

        double fitness = 0;/*
        for (auto & ap_period: ap_period_baselines) {
            std::vector<double> apRecord(ap_period.first.size());
            const double period = ap_period.second;
            //prepare model for a new CL
            //TODO
            int is_correct;
            solver.solve(model, std::vector<double> & y0, is_correct,
				double t0, double start_record, double tout, apRecord);
            //save something you want back to parameters
            //TODO
            
            double dist = 0;
            for (int i = 0; i < apRecord.size(); i++)
                dist += pow(apRecord[i] - ap_period.first[i], 2);
            fitness += dist / apRecord.size();
        }*/
        //make sure parameters were written back
        //TODO
        return fitness;
    }
    template <typename It>
    void genetic_algorithm_result(It parameters_begin)
    {
        //TODO
     //   for (int i = 0; i != func.get_xdim(); i++, parameters_begin++)
        //    result[i] = *parameters_begin;
    }

};




template <typename Func>
class MinimizeFunc
{
public:
    using num = typename Func::num;
    using value = typename Func::value;
    num operator()(const value & val) const
    {
        return val[0];
    }
};

template <typename Func>
class MaximizeFunc
{
public:
    using num = typename Func::num;
    using value = typename Func::value;
    num operator()(const value & val) const
    {
        return -val[0];
    }
};

template <typename Func, template<typename F> class Obj>
class FuncOptimization
{
    Func & func;
    Obj<Func> obj;
    using argument = typename Func::argument;
    argument result;
public:
    FuncOptimization(Func & func)
    : func(func)
    {}
    argument get_result()
    {
        return result;
    }
    int get_number_parameters() const
    {
        return func.get_xdim();
    }
    template <typename It>
    int get_boundaries(It pmin, It pmax)
    {
        auto ppmin = func.x_min();
        std::copy(ppmin.begin(), ppmin.end(), pmin);
        auto ppmax = func.x_max();
        std::copy(ppmax.begin(), ppmax.end(), pmax);
        return 0; //have boundaries
        //return -1; //no boundaries
    }
    template <typename It>
    int initial_guess(It begin)
    {
        //return 0; //if you provided some initial guess
        return -1; //no initial guess at all
    }
    template <typename It>
    double genetic_algorithm_calls(It parameters_begin)
    {
        argument params;
        for (int i = 0; i != func.get_xdim(); i++, parameters_begin++)
            params[i] = *parameters_begin;
        return obj(func(params));
    }
    template <typename It>
    void genetic_algorithm_result(It parameters_begin)
    {
        for (int i = 0; i != func.get_xdim(); i++, parameters_begin++)
            result[i] = *parameters_begin;
    }
};

#endif
