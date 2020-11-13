#ifndef OPTIMIZATION_PROBLEM
#define OPTIMIZATION_PROBLEM

#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>

#include <json.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/bimap/unconstrained_set_of.hpp>

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
        for (size_t i = 0; i < a.size(); i++)
            res += std::pow(std::abs(a[i] - b[i]), power);
        res /= a.size();
        return res;
    }

    double dist(const BaselineVec & a, const BaselineVec & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
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
    std::vector<Variables> baselineVariables;

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
    
    using BiMap = boost::bimap<boost::bimaps::unordered_set_of<std::string>,
            boost::bimaps::unordered_set_of<int>,
            boost::bimaps::list_of_relation>;
            
    BiMap constantsBiMapModel, statesBiMapModel;
      //no need of algebraicBiMapModel, ratesBiMapModel;
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

        std::unordered_map<int, std::string> statesMap(model.state_size()),
                constantsMap(model.constants_size()),
                algMap(model.get_alg_size()),
                ratesMap(model.state_size());

        model.get_maps(statesMap, constantsMap, algMap, ratesMap);

        //now save it to bimap
        for (auto & s: statesMap) {
            s.second.resize(s.second.find(" "));
            statesBiMapModel.left.insert(BiMap::left_value_type(s.second, s.first));
        }
        for (auto & s: constantsMap) {
            s.second.resize(s.second.find(" "));
            constantsBiMapModel.left.insert(BiMap::left_value_type(s.second, s.first));
        }
    }

    int get_number_parameters() const
    {
        return number_parameters;
    }
    
    void scan_baseline(Baseline & baseline, const std::string & filename)
    {
        std::ifstream file;
        file.open(filename);
        if (!file.is_open())
            throw(filename + " cannot be opened");
        double v;
        while (file >> v)
            baseline.push_back(v);
    }
    
    void scan_state()
    {
        //TODO
    }
    
    void read_config(const std::string & configFilename)
    {
        if (mpi_rank == 0) {
            std::ifstream configFile;
            configFile.open(configFilename);
            if (!configFile.is_open())
                throw("Cannot open main config file");
            json config;
            configFile >> config;
            configFile.close();
            
            number_parameters = 0;
            //read global variables
            for (auto variable: config["global"].items()) {
                auto v = variable.value();
                std::string name = v["name"].get<std::string>();
                
                int model_position;
                try {
                    //usually in global we have constants
                    model_position = constantsBiMapModel.left.at(name);
                } catch (...) {
                    //so name is not a constant
                    //1. it can be a state (which is quite reasonable but we dont need it now)
                    //2. or it is neither a constant nor a state variable
                    throw(name + " in global config is not a constant of a model");
                }
                bool is_value = (v.find("value") != v.end());
                if (is_value) {
                    globalVariables.valueConstants.push_back(
                    {.name = name,
                     .value = v["value"].get<double>(),
                     .model_position = model_position
                    });
                } else {
                    //parameter
                    globalVariables.mutableConstants.push_back(
                    {.name = name,
                     .min_value = v["bounds"][0].get<double>(),
                     .max_value = v["bounds"][1].get<double>(),
                     .parameter_position = number_parameters++,
                     .model_position = model_position,
                     .gamma = v["gamma"].get<double>(),
                     .is_mutation_applicable = 1
                    });
                }
            }
            //constants not listed in config will have default values


            for (auto baseline: config["baselines"].items()) {
                auto b = baseline.value();
                std::string apfilename = b["filename_phenotype"].get<std::string>();

                apbaselines.emplace_back();
                scan_baseline(apbaselines.back(), apfilename);

                bool initial_state = (b.find("filename_state") != b.end());
                if (initial_state) {
                    std::string statefilename = b["filename_state"].get<std::string>();

                    //TODO
                    //scan_state(..., statefilename);
                }


                baselineVariables.emplace_back();
                Variables & bVar = baselineVariables.back();
                BiMap statesBiMapDrifting = statesBiMapModel;
                for (auto variable: b["params"].items()) {
                    auto v = variable.value();
                    std::string name = v["name"].get<std::string>();
                    int model_position;
                    try {
                        //constant can be listed as a baseline value
                        model_position = constantsBiMapModel.left.at(name);
                        bool is_value = (v.find("value") != v.end());
                        if (!is_value)
                            throw(name + " in baseline config is required to be a value");
                        bVar.valueConstants.push_back(
                        {.name = name,
                         .value = v["value"].get<double>(),
                         .model_position = model_position
                        });
                    } catch (...) {
                        //name is not a constant
                        //so it is probably a state which can be mutated
                        try {
                            model_position = statesBiMapModel.left.at(name);
                        } catch (...) {
                            //or it is neither a constant nor a state variable
                            throw(name + " in baseline config is neither constant nor state variable");
                        }
                        bool is_value = (v.find("value") != v.end());
                        if (is_value)
                            throw(name + "cannot be a value");
                        //remove it from statesBiMapDrifting
                        statesBiMapDrifting.left.erase(name);

                        bVar.mutableStates.push_back(
                        {.name = name,
                         .min_value = v["bounds"][0].get<double>(),
                         .max_value = v["bounds"][1].get<double>(),
                         .parameter_position = number_parameters++,
                         .model_position = model_position,
                         .gamma = v["gamma"].get<double>(),
                         .is_mutation_applicable = 1
                        });
                    }
                }
                //other state mutable variables for each baseline
                //are drifting
                for (const auto & sit: statesBiMapDrifting) {
                    std::cout << sit.left << std::endl;
                    bVar.mutableStates.push_back(
                    {.name = sit.left,
                     .min_value = 0,
                     .max_value = 0,
                     .parameter_position = number_parameters++,
                     .model_position = sit.right,
                     .gamma = 0,
                     .is_mutation_applicable = 0
                    });
                }
            }
            std::cout << "Unknowns: " << number_parameters << std::endl;
            throw;
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
        for (const Variables & vars: baselineVariables) {
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
            std::vector<double> vconstants(model.constants_size());

            //first, find it
            model.initConsts(vconstants.data());
            model.initState(y0.data());
            
            //then, set mutable global parameters
            for (const Mutable & m: globalVariables.mutableConstants) {
                parameters_begin[m.parameter_position] = vconstants[m.model_position];
            }
            
            //finally, set mutable baseline parameters
            for (size_t i = 0; i < apbaselines.size(); i++) {
                for (const Mutable & m: baselineVariables[i].mutableStates) {
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
        for (const Mutable & m: baselineVariables[i].mutableConstants) {
            throw("Wait, that's illegal. Please contact us if you want it.");
        }
        for (const Mutable & m: baselineVariables[i].mutableStates) {
            y0[m.model_position] = parameters_begin[m.parameter_position];
        }
        for (const Value & v: baselineVariables[i].valueConstants) {
            constants[v.model_position] = v.value;
        }
        for (const Value & v: baselineVariables[i].valueStates) {
            throw("Wait, that's illegal. Please contact us if you want it.");
        }
    }

    template <typename It>
    double genetic_algorithm_calls(It parameters_begin)
    {
        for (size_t i = 0; i < apbaselines.size(); i++) {
            std::vector<double> y0(model.state_size());
            std::vector<double> vconstants(model.constants_size());
            double * constants = vconstants.data();

            Baseline & apRecord = apmodel[i];
            model.set_constants(constants);

            fill_constants_y0(parameters_begin, constants, y0.data(), i);

            const double period = vconstants[constantsBiMapModel.left.at("stim_period")];
            
            int is_correct;
            const int beats = 10;
            const double t0 = 0, start_record = period * (beats - 1),
                    tout = period * beats;

            solver.solve(model, y0, is_correct,
                            t0, start_record, tout, apRecord);

            //maybe if a solver fails rarely then we may consider it
            //as a really poor fitness
            if (!is_correct)
                throw("Solver failed");

            //save mutable variables from the state
            //since we let them to drift
            for (const Mutable & m: baselineVariables[i].mutableStates) {
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
        for (size_t i = 0; i < apbaselines.size(); i++) {
            BaselineResult res;
            std::vector<double> y0(model.state_size());
            std::vector<double> vconstants(model.constants_size());
            fill_constants_y0(parameters_begin, vconstants.data(), y0.data(), i);

            for (const auto & cit: constantsBiMapModel) {
                res.constantsResult[cit.left] = vconstants[cit.right];
            }
            for (const auto & sit: statesBiMapModel) {
                res.statesResult[sit.left] = y0[sit.right];
            }
            results.push_back(res);
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
