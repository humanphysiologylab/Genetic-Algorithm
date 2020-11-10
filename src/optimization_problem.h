#ifndef OPTIMIZATION_PROBLEM
#define OPTIMIZATION_PROBLEM

#include <string>
#include <vector>
/*
template <typename Model, typename Solver>
class ODEoptimization
{
    Model & model;
    Solver & solver;
    std::vector<std::pair<std::vector<double>, double>> ap_period_baselines;
    int number_parameters = 0;
public:
    ODEoptimization(Model & model, Solver & solver)
    : model(model), solver(solver)
    {
        
    }
    
    int get_number_parameters() const
    {
        return number_parameters;
    }
    void read_config()
    {
        number_parameters = 100;
    }
    int get_boundaries()
    {
        return -1; //no boundaries
    }
    int initial_guess(double * p)
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


    double genetic_algorithm_calls(double * parameters)
    {
        //unpack parameters
        //TODO

        double fitness = 0;
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
        }
        //make sure parameters were written back
        //TODO
        return fitness;
    }
};
*/



template <typename Func>
class MinimizeFunc
{
public:
    typename Func::num operator()(const typename Func::value & value) const
    {
        return value[0];
    }
};

template <typename Func>
class MaximizeFunc
{
public:
    typename Func::num operator()(const typename Func::value & value) const
    {
        return -value[0];
    }
};

template <typename Func, template<typename F> class Obj>
class FuncOptimization
{
    Func & func;
    Obj<Func> obj;
    typename Func::argument result;
public:
    FuncOptimization(Func & func)
    : func(func)
    {}
    typename Func::argument get_result()
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
        return -1; //no initial guess at all
    }
    template <typename It>
    double genetic_algorithm_calls(It parameters_begin)
    {
        typename Func::argument params;
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
