#ifndef OPTIMIZATION_PROBLEM
#define OPTIMIZATION_PROBLEM

#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <limits>

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
        res = std::pow(res, 1.0/power);
        if (std::isnan(res))
            res = std::numeric_limits<double>::max();
        return res;
    }
};

template <typename Model, typename Solver, typename Objective>
class ODEoptimization
{
    Model glob_model;//!!!!
    Solver solver;//!!!!
    Objective obj;//!!!!
    using Baseline = typename Objective::Baseline;
    using BaselineVec = typename Objective::BaselineVec;
    
    BaselineVec apbaselines;
    
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

    int beats;
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
    Results results, relative_results;//abs values and relative to default values
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
    Results get_relative_results() const
    {
        return relative_results;
    }
    std::vector<double> get_results_optimizer_format() const
    {
        return results_optimizer_format;
    }
    ODEoptimization(const Model & model, const Solver & solver, const Objective & obj)
    : glob_model(model), solver(solver), obj(obj)
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

    void write_baseline(const Baseline & apRecord, const std::string & filename) const
    {
        std::ofstream file;
        file.open(filename);
        if (!file.is_open())
            throw(filename + " cannot be opened");
        for (const auto & x: apRecord)
            file << x << std::endl;
    }
    
    void write_state(const std::vector<double> & state, const std::string & filename) const
    {
        std::ofstream file;
        file.open(filename);
        if (!file.is_open())
            throw(filename + " cannot be opened");
        for (const auto & x: state)
            file << x << std::endl;
    }
    void read_config(const std::string & configFilename)
    {
        //lets make it plain, no MPI-IO, TODO
       // if (mpi_rank == 0) {
            std::ifstream configFile;
            configFile.open(configFilename);
            if (!configFile.is_open())
                throw("Cannot open main config file");
            json config;
            configFile >> config;
            configFile.close();


            beats = config["n_beats"].get<int>();

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

            int baselines_num = 0;
            for (auto baseline: config["baselines"].items()) {
                baselines_num++;
                auto b = baseline.value();

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
            
            //we may need to generate test baselines first
            bool run_type = (config.find("mode") != config.end());
            apbaselines = BaselineVec(baselines_num); //sorry for this
            if (mpi_rank == 0 && run_type && config["mode"].get<std::string>() == "test") {
                
                std::cout << "Generating test baselines" << std::endl;
                std::cout << "ATTENTION! ALL BASELINES FROM CONFIG WILL BE OVERWRITTEN" << std::endl;
                const double t_sampling = config["t_sampling"].get<double>();
                const int beats_num = config["test_beats"].get<int>();//run model long enough to get a steady state

                int baseline_number = 0;

                for (auto baseline: config["baselines"].items()) {
                    std::vector<double> y0(glob_model.state_size());
                    Model model(glob_model);
                    std::vector<double> parameters(number_parameters);
                    initial_guess(parameters.begin());

                    std::vector<double> vconstants(model.constants_size());
                    double * constants = vconstants.data();
                    model.set_constants(constants);//!!!!!!!
                    fill_constants_y0(parameters.begin(), constants, y0.data(), baseline_number);
                    const double period = vconstants[constantsBiMapModel.left.at("stim_period")];

                    //now we need to resize apRecord according to period lenght and t_sampling
                    Baseline apRecord(1 + std::ceil( period / t_sampling ));
                    model_eval(y0, model, beats_num, period, apRecord);

                    if (std::isnan(apRecord[0])) {
                        std::cout << "nan" << std::endl;
                        throw;
                    }
                    //save state and baseline
                    auto b = baseline.value();
                    std::string apfilename = b["filename_phenotype"].get<std::string>();
                    write_baseline(apRecord, apfilename);
                    //initial state may not be provided
                    bool initial_state = (b.find("filename_state") != b.end());
                    if (initial_state) {
                        std::string statefilename = b["filename_state"].get<std::string>();
                        write_state(y0, statefilename);
                    }
                    baseline_number++;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);

            //read baselines
            for (auto baseline: config["baselines"].items()) {
                auto b = baseline.value();
                std::string apfilename = b["filename_phenotype"].get<std::string>();
                apbaselines = BaselineVec();
                apbaselines.emplace_back();
                scan_baseline(apbaselines.back(), apfilename);

                //initial state may not be provided
                bool initial_state = (b.find("filename_state") != b.end());
                if (initial_state) {
                    std::string statefilename = b["filename_state"].get<std::string>();

                    //TODO
                    //scan_state(..., statefilename);
                }
            }
       //     std::cout << "Unknowns: " << number_parameters << std::endl;
      //  }
        //broadcast all this nonsense;
        //TODO
    }
    template <typename It>
    int get_boundaries(It pmin, It pmax, int * is_mutation_applicable) const
    {
        for (const Mutable & gl: globalVariables.mutableConstants) {
            pmin[gl.parameter_position] = gl.min_value;
            pmax[gl.parameter_position] = gl.max_value;
            is_mutation_applicable[gl.parameter_position] = gl.is_mutation_applicable;
        }
        for (const Mutable & gl: globalVariables.mutableStates) {
            pmin[gl.parameter_position] = gl.min_value;
            pmax[gl.parameter_position] = gl.max_value;
            is_mutation_applicable[gl.parameter_position] = gl.is_mutation_applicable;
        }
        for (const Variables & vars: baselineVariables) {
            for (const Mutable & gl: vars.mutableConstants) {
                pmin[gl.parameter_position] = gl.min_value;
                pmax[gl.parameter_position] = gl.max_value;
                is_mutation_applicable[gl.parameter_position] = gl.is_mutation_applicable;
            }
            for (const Mutable & gl: vars.mutableStates) {
                pmin[gl.parameter_position] = gl.min_value;
                pmax[gl.parameter_position] = gl.max_value;
                is_mutation_applicable[gl.parameter_position] = gl.is_mutation_applicable;
            }
        }
        return 0;
    }
    template <typename It>
    int initial_guess(It parameters_begin) const
    {
        if (0) {
            //read from a file
            //TODO
            return 0;
        } else if (1) {
            //use model's default
            std::vector<double> y0(glob_model.state_size());
            std::vector<double> vconstants(glob_model.constants_size());

            //first, find it
            glob_model.initConsts(vconstants.data());
            glob_model.initState(y0.data());
            
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

    void get_default_values(double * constants, double * y0) const
    {
        glob_model.initConsts(constants);
        glob_model.initState(y0);
    }

    template <typename It>
    void fill_constants_y0(It parameters_begin, double * constants, double * y0, size_t i) const
    {
        //first, default
        get_default_values(constants, y0);

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

    void model_eval(std::vector<double> & y0, Model & model, int beats_number, double period, Baseline & apRecord) const
    {
        int is_correct;

        const double t0 = 0, start_record = period * (beats_number - 1),
                tout = period * beats_number;

        solver.solve(model, y0, is_correct,
                        t0, start_record, tout, apRecord);

        //maybe if a solver fails rarely then we may consider it
        //as a really poor fitness
        if (!is_correct)
            throw("Solver failed");
    }

    template <typename It>
    double genetic_algorithm_calls(It parameters_begin) const
    {
        BaselineVec apmodel;
        for (const auto & v: apbaselines)
            apmodel.push_back(Baseline(v.size()));
        
        Model model(glob_model);
        for (size_t i = 0; i < apbaselines.size(); i++) {
            std::vector<double> y0(model.state_size());
            Baseline & apRecord = apmodel[i];

            std::vector<double> vconstants(model.constants_size());
            double * constants = vconstants.data();
            model.set_constants(constants);//!!!!!!!
            fill_constants_y0(parameters_begin, constants, y0.data(), i);
            const double period = vconstants[constantsBiMapModel.left.at("stim_period")];
            

            model_eval(y0, model, beats, period, apRecord);
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
            BaselineResult res, relative_res;
            std::vector<double> y0(glob_model.state_size());
            std::vector<double> vconstants(glob_model.constants_size());
            fill_constants_y0(parameters_begin, vconstants.data(), y0.data(), i);

            std::vector<double> y0_default(glob_model.state_size());
            std::vector<double> vconstants_default(glob_model.constants_size());
            get_default_values(vconstants_default.data(), y0_default.data());
            //TODO y0_default should be specific for each baseline
            //not really from a default model

            for (const auto & cit: constantsBiMapModel) {
                res.constantsResult[cit.left] = vconstants[cit.right];
                relative_res.constantsResult[cit.left] = vconstants[cit.right] / vconstants_default[cit.right];
            }
            for (const auto & sit: statesBiMapModel) {
                res.statesResult[sit.left] = y0[sit.right];
                relative_res.statesResult[sit.left] = y0[sit.right] / y0_default[sit.right];
            }
            results.push_back(res);
            relative_results.push_back(relative_res);
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
    argument get_result() const
    {
        return result;
    }
    int get_number_parameters() const
    {
        return func.get_xdim();
    }
    template <typename It>
    int get_boundaries(It pmin, It pmax, int * is_mutation_applicable) const
    {
        auto ppmin = func.x_min();
        std::copy(ppmin.begin(), ppmin.end(), pmin);
        auto ppmax = func.x_max();
        std::copy(ppmax.begin(), ppmax.end(), pmax);
        std::fill(is_mutation_applicable, is_mutation_applicable + get_number_parameters(), 1);
        return 0; //have boundaries
        //return -1; //no boundaries
    }
    template <typename It>
    int initial_guess(It begin) const
    {
        //return 0; //if you provided some initial guess
        return -1; //no initial guess at all
    }
    template <typename It>
    double genetic_algorithm_calls(It parameters_begin) const
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
