#ifndef OPTIMIZATION_PROBLEM
#define OPTIMIZATION_PROBLEM

#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>

#include <json.hpp>
using json = nlohmann::json;

template<int power = 2>
class MinimizeAPbaselines
{
public:
    using Baseline = std::vector<double>;
    using BaselineVec = std::vector<Baseline>;

    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        double res = 0;
        assert(a.size() == b.size());
        for (int i = 0; i < a.size(); i++)
            res += std::pow(std::abs(a[i] - b[i]), power);
        res /= a.size();
        return res;
    }

    double dist(const BaselineVec & a, const BaselineVec & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (int i = 0; i < a.size(); i++) {
            res += distBaselines(a[i], b[i]);
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
    using BaselineVec = typename Objective::BaselineVec;
    
    BaselineVec apbaselines, apmodel;
    
    //two types of variables
    struct Mutable
    {
        std::string name;
        double min_value;
        double max_value;
        int parameter_position;
        int model_position;
        double gamma;
        int is_mutation_applicable;
    };
    struct Value
    {
        std::string name;
        double value;
        int model_position;
    };

    struct Variables
    {
        std::vector<Mutable> mutableConstants;
        std::vector<Mutable> mutableStates;
        std::vector<Value>   valueConstants;
        std::vector<Value>   valueStates;
    };

    //global section of config
    Variables globalVariables;

    //baseline sections of config
    std::vector<Variables> BaselineVariables;

    //number of mutable variables passed to an optimization algorithm
    int number_parameters = 0;

    int mpi_rank;
    int mpi_size;
public:
    using ConstantsResult = std::unordered_map<std::string, double>;
    using StatesResult = std::unordered_map<std::string, double>;
    struct BaselineResult
    {
        ConstantsResult constantsResult;
        StatesResult statesResult;
    };
    using Results = std::vector<BaselineResult>;
private:
    Results results;
    std::vector<double> results_optimizer_format;
    std::unordered_map<std::string, int> constantsMapModel, statesMapModel;
public:

    Results get_results() const
    {
        return results;
    }
    std::vector<double> get_results_optimizer_format() const
    {
        return results_optimizer_format;
    }
    ODEoptimization(Model & model, Solver & solver, Objective & obj)
    : model(model), solver(solver), obj(obj)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    }

    int get_number_parameters() const
    {
        return number_parameters;
    }
    void read_config()
    {
        if (mpi_rank == 0) {
            std::ifstream configFile;
            configFile.open("config.json");
            if (!configFile.is_open())
                throw("Cannot open main config file");
            json config;

            configFile >> config;
            configFile.close();
            //TODO
        }
        //broadcast all this nonsense;
        //TODO
    }
    template <typename It>
    int get_boundaries(It pmin, It pmax)
    {
        for (const Mutable & gl: globalVariables.mutableConstants) {
            pmin[gl.model_position] = gl.min_value;
            pmax[gl.model_position] = gl.max_value;
        }
        for (const Mutable & gl: globalVariables.mutableStates) {
            pmin[gl.model_position] = gl.min_value;
            pmax[gl.model_position] = gl.max_value;
        }
        for (const Variables & vars: BaselineVariables) {
            for (const Mutable & gl: vars.mutableConstants) {
                pmin[gl.model_position] = gl.min_value;
                pmax[gl.model_position] = gl.max_value;
            }
            for (const Mutable & gl: vars.mutableStates) {
                pmin[gl.model_position] = gl.min_value;
                pmax[gl.model_position] = gl.max_value;
            }
        }
        return 0;
    }
    template <typename It>
    int initial_guess(It parameters_begin)
    {
        if (0) {
            //read from a file
            //TODO
            return 0;
        } else if (1) {
            //use model's default
            std::vector<double> y0(model.state_size());
            std::vector<double> vconstants(model.constants_size()];

            //first, find it
            model.initConsts(constants.data());
            model.initState(y0.data());
            
            //then, set mutable global parameters
            for (const Mutable & m: globalVariables.mutableConstants) {
                parameters_begin[m.parameter_position] = constants[m.model_position];
            }
            
            //finally, set mutable baseline parameters
            for (int i = 0; i < apbaselines.size(); i++) {
                for (const Mutable & m: BaselineVariables[i].mutableStates) {
                    parameters_begin[m.parameter_position] = y0[m.model_position];
                }
            }
            return 0;
        } else {
            //or just not provide anything
            return -1;//initial parameters are not provided
        }
    }

    template <typename It>
    void fill_constants_y0(It parameters_begin, double * constants, double * y0, int i)
    {
        //first, default
        model.initConsts(constants);
        model.initState(y0);
            
        //global section of config overwrites default
        for (const Mutable & m: globalVariables.mutableConstants) {
            constants[m.model_position] = parameters_begin[m.parameter_position];
        }
        for (const Mutable & m: globalVariables.mutableStates) {
            throw("Wait, that's illegal. Please contact us if you want it.");
        }
        for (const Value & v: globalVariables.valueConstants) {
            constants[v.model_position] = v.value;
        }
        for (const Value & v: globalVariables.valueStates) {
            throw("Wait, that's illegal. Please contact us if you want it.");
        }

        //baseline sections of config overwrites global and default
        for (const Mutable & m: BaselineVariables[i].mutableConstants) {
            throw("Wait, that's illegal. Please contact us if you want it.");
        }
        for (const Mutable & m: BaselineVariables[i].mutableStates) {
            y0[m.model_position] = parameters_begin[m.parameter_position];
        }
        for (const Value & v: BaselineVariables[i].valueConstants) {
            constants[v.model_position] = v.value;
        }
        for (const Value & v: BaselineVariables[i].valueStates) {
            throw("Wait, that's illegal. Please contact us if you want it.");
        }
    }

    template <typename It>
    double genetic_algorithm_calls(It parameters_begin)
    {
        for (int i = 0; i < apbaselines.size(); i++) {
            std::vector<double> y0(model.state_size());
            std::vector<double> vconstants(model.constants_size()];
            double * constants = vconstants.data();

            Baseline & apRecord = apmodel[i];
            model.set_constants(constants);

            fill_constants_y0(parameters_begin, constants, y0.data(), i);

            int is_correct;
            const int beats = 10;
            const double t0 = 0, start_record = period * (beats - 1),
                    tout = period * beats;

            solver.solve(model, std::vector<double> & y0, is_correct,
                            t0, start_record, tout, apRecord);

            //maybe if a solver fails rarely then we may consider it
            //as a really poor fitness
            if (!is_correct)
                throw("Solver failed");

            //save mutable variables from the state
            //since we let them to drift
            for (const Mutable & m: BaselineVariables[i].mutableStates) {
                parameters_begin[m.parameter_position] =  y0[m.model_position];
            }
        }
        return obj.dist(apbaselines, apmodel);
    }
    template <typename It>
    void genetic_algorithm_result(It parameters_begin)
    {
        //mirror the stage of initialization before solver.solve call
        //but just save the final result
        for (int i = 0; i < apbaselines.size(); i++) {
            BaselineResult res;
            std::vector<double> y0(model.state_size());
            std::vector<double> vconstants(model.constants_size()];
            fill_constants_y0(parameters_begin, vconstants.data(), y0.data(), i);

            for (auto & cit: constantsMapModel) {
                res.constantsResult[cit->first] = vconstants[cit->second];
            }
            for (auto & sit: statesMapModel) {
                res.statesResult[sit->first] = y0[sit->second];
            }
            results.append(res);
        }
        results_optimizer_format = std::vector<double>(number_parameters);
        std::copy(parameters_begin, parameters_begin + number_parameters, results_optimizer_format.begin());
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
