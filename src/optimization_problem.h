#ifndef OPTIMIZATION_PROBLEM
#define OPTIMIZATION_PROBLEM

#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <limits>
#include <mpi.h>
#include <cmath>
#include <json.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/bimap/unconstrained_set_of.hpp>

#include <Eigen/Dense>

using json = nlohmann::json;

class BlockOfTable
: public Eigen::DenseBase<Eigen::MatrixXd >::ColsBlockXpr
{
    int * model_indices;
    using Base = Eigen::DenseBase<Eigen::MatrixXd >::ColsBlockXpr;
public:
    BlockOfTable(int * model_indices, Block xblock);
    int get_model_pos(int col) const;
};

class Table
: public Eigen::MatrixXd
{
    using Base = Eigen::MatrixXd;
    int states_cols, alg_cols;
    int * model_indices;
    std::vector<std::string> header;
public:
    Table(int m, int n, std::vector<int> states_model_indices,
            std::vector<int> alg_model_indices,
            const std::vector<std::string> & states_names,
            const std::vector<std::string> & alg_names);
    ~Table();
    BlockOfTable get_algebraic();
    BlockOfTable get_states();
    void export_csv(const std::string & filename);
};


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
        res /= (double) a.size();
        return res;
    }

    double dist(const BaselineVec & a, const BaselineVec & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += distBaselines(a[i], b[i]);
        }
        res = std::pow(res / (double) a.size(), 1.0/power);
        if (std::isnan(res))
            res = 1e50;
        return res;
    }
};


template<int power = 2>
class LeastSquaresMinimizeAPbaselines
{
public:
    using Baseline = std::vector<double>;
    using BaselineVec = std::vector<Baseline>;
    
    template<int p>
    double sum_of_elements(const Baseline & b) const
    {
        double s = 0;
        for (const auto & e: b)
            s += std::pow(e, p);
        return s;
    }
    double inner_product(const Baseline & a, const Baseline & b) const
    {
        double s = 0;
        for (unsigned i = 0; i < a.size(); i++)
            s += a[i] * b[i];
        return s;
    }
    
    double LSMdistBaselines(const Baseline & g, const Baseline & f) const
    {
        //watch for the correct order: g -- experimental optical mapping AP, f -- model
        //solve least square problem for f(t) = A * g(t) + B
        assert(g.size() == f.size());

        const double c = g.size();
        const double a = sum_of_elements<2>(g);
        const double b = sum_of_elements<1>(g);
        const double e = sum_of_elements<1>(f);
        const double p = inner_product(f, g);

        const double A = (c * p  - b * e) / (a * c - b * b);
        const double B = (-b * p + a * e) / (a * c - b * b);

        const double min_peak = 10;         // mV
        const double max_peak = 60;         // mV
        const double min_rest_potential = -100;   // mV
        const double max_rest_potential = -60;   // mV
        
        const double peak = std::min(std::max(A + B, min_peak), max_peak);
        const double rest_potential = std::min(std::max(B, min_rest_potential), max_rest_potential);

        const double amplitude = peak - rest_potential;
        /*
        std::cout << "     A+B: " << A + B << " A: " << amplitude << " B: " << rest_potential << " p: " << peak << std::endl;
        */
        if (std::isnan(amplitude)) {
            std::cout << "         a:" << a <<
            "         b:" << b <<
            "         c:" << c <<
            "         e:" << e <<
            "         p:" << p <<
            "         A:" << A <<
            "         B:" << B << std::endl;
            for (int i = 0; i < f.size(); i++) {
                if (std::isnan(f[i])) {
                    std::cout << "f is nan" << std::endl;
                    break;
                }
            }
        }
        double res = 0;
        for (size_t i = 0; i < f.size(); i++)
            res += std::pow(std::abs(f[i] - (amplitude * g[i] + rest_potential)), power);
            
       // for (size_t i = 40; i < 80; i++)
         //   res += std::pow(std::abs(f[i] - (amplitude * g[i] + rest_potential)), 4);

        res /= (double) f.size();
        return res;
    }

    double dist(const BaselineVec & a, const BaselineVec & b) const
    {
        //watch for the correct order: a -- experimental optical mapping AP, b -- model
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += LSMdistBaselines(a[i], b[i]);
        }
        res = std::pow(res / (double) a.size(), 1.0/power);
        if (std::isnan(res))
            res = 1e50;
        return res;
    }
};

template<int power = 2>
class MaximizeAPinnerProduct
{
public:
    using Baseline = std::vector<double>;
    using BaselineVec = std::vector<Baseline>;

    double innerProduct(const Baseline & a, const Baseline & b) const
    {
        double res = 0, lena = 0, lenb = 0;
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); i++) {
            res += a[i] * b[i];
            lena += a[i] * a[i];
            lenb += b[i] * b[i];
        }
        res /= sqrt(lena * lenb);
        return res;
    }

    double dist(const BaselineVec & a, const BaselineVec & b) const
    {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += innerProduct(a[i], b[i]);
        }
        res = res / (double) a.size(); //mean value
        if (res < -1 || res > 1)
            std::cout << "ATTENTION: inner product: " << res << std::endl;
        if (std::isnan(res))
            return 1e50;//std::numeric_limits<double>::max();
        return -res;
    }
};






template <typename Model, typename Solver, typename Objective>
class ODEoptimization
{
protected:
    Model glob_model;//!!!!//Maybe not to use it at all due to the multithreading?
    Solver solver;//!!!!
    Objective obj;//!!!!
    using Baseline = typename Objective::Baseline;
    using BaselineVec = typename Objective::BaselineVec;
    
    BaselineVec apbaselines;
    
    //two types of variables
    
    //TODO
    //maybe there should be a dict of parameters
    //required only by a method, and a method provides the keys
    //by itself
    struct Mutable
    {
        std::string name;
        double min_value;
        double max_value;
        int init_from_config = 0;
        double init_guess;
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

    int mpi_rank, mpi_size;

    int beats;
    bool is_AP_normalized = 0; //normalize every baseline vector independently!!!
public:
    using ConstantsResult = std::unordered_map<std::string, double>;
    using StatesResult = std::unordered_map<std::string, double>;
    struct BaselineResult
    {
        ConstantsResult constantsResult;
        StatesResult statesResult;
    };
    using Results = std::vector<BaselineResult>;
protected:
    Results results, relative_results;//abs values and relative to default values
    std::vector<double> results_optimizer_format;
    
    using BiMap = boost::bimap<boost::bimaps::unordered_set_of<std::string>,
            boost::bimaps::unordered_set_of<int>,
            boost::bimaps::list_of_relation>;
            
    BiMap constantsBiMapModel, statesBiMapModel, algebraicBiMapModel;
    //no need of ratesBiMapModel

    double log_scale(double x, double minOrig, double maxOrig) const
    {
        if (minOrig < 0) {
            throw("minOrig < 0");
        } else if (minOrig == 0) {
            return log(x + 1.0 / maxOrig) / log(maxOrig + 1.0 / maxOrig);
        } else {
            const double a = log(minOrig);
            const double b = log(maxOrig);
            return (log(x) - (a + b) / 2) / ((b - a) / 2);
        }
    }
    double log_scale_back(double x, double minOrig, double maxOrig) const
    {
        if (minOrig < 0) {
            throw("minOrig < 0");
        } else if (minOrig == 0) {
            return exp(x * log(maxOrig + 1.0 / maxOrig)) - 1.0 / maxOrig;
        } else {
            const double a = log(minOrig);
            const double b = log(maxOrig);
            return exp(x * ((b - a) / 2) + (a + b) / 2);
        }
    }
    double lin_scale(double x, double minOrig, double maxOrig) const
    {
        return (x - minOrig) / (maxOrig - minOrig);
    }
    double lin_scale_back(double x, double minOrig, double maxOrig) const
    {
        return (maxOrig - minOrig) * x + minOrig;
    }
    template <typename It, typename It2>
    void optimizer_model_scale(It optim_param_start, It2 model_param_start) const
    {
        //from optimizer to model
        for (const Mutable & m: globalVariables.mutableConstants) {
            const int pos = m.parameter_position;
            if (m.is_mutation_applicable == 1) {
                model_param_start[pos] = lin_scale_back(optim_param_start[pos], m.min_value, m.max_value);   
            } else if (m.is_mutation_applicable == 2) {
                model_param_start[pos] = log_scale_back(optim_param_start[pos], m.min_value, m.max_value); 
            } else {
                throw(m.name + ": unknown scale type");
            }
        }
        for (const auto & bv: baselineVariables) {
            for (const Mutable & m: bv.mutableStates) {
                const int pos = m.parameter_position;
                if (m.is_mutation_applicable == 1) {
                    model_param_start[pos] = lin_scale_back(optim_param_start[pos], m.min_value, m.max_value);   
                } else if (m.is_mutation_applicable == 2) {
                    model_param_start[pos] = log_scale_back(optim_param_start[pos], m.min_value, m.max_value); 
                } else {
                    throw(m.name + ": unknown scale type");
                }
            }
        }
    }

    template <typename It, typename It2>
    void model_optimizer_scale(It model_param_start, It2 optim_param_start) const
    {
        //from model to optimizer
        for (const Mutable & m: globalVariables.mutableConstants) {
            const int pos = m.parameter_position;
            if (m.is_mutation_applicable == 1) {
                optim_param_start[pos] = lin_scale(model_param_start[pos], m.min_value, m.max_value);   
            } else if (m.is_mutation_applicable == 2) {
                optim_param_start[pos] = log_scale(model_param_start[pos], m.min_value, m.max_value); 
            } else {
                throw(m.name + ": unknown scale type");
            }
        }
        for (const auto & bv: baselineVariables) {
            for (const Mutable & m: bv.mutableStates) {
                const int pos = m.parameter_position;
                if (m.is_mutation_applicable == 1) {
                    optim_param_start[pos] = lin_scale(model_param_start[pos], m.min_value, m.max_value);   
                } else if (m.is_mutation_applicable == 2) {
                    optim_param_start[pos] = log_scale(model_param_start[pos], m.min_value, m.max_value); 
                } else {
                    throw(m.name + ": unknown scale type");
                }
            }
        }
    }
public:
    void unfreeze_global_variable(const std::string & name, double min_value, double max_value, std::vector<double> & params)
    {//TODO more detailed! Use it at your own risk!
        globalVariables.valueConstants;//remove from valueConstants
        bool is_found_in_valueConstants = 0;
        Value v;
        for (auto wi = globalVariables.valueConstants.begin();
            wi != globalVariables.valueConstants.end(); wi++) {
            if (wi->name == name) {
                v = *wi;
                is_found_in_valueConstants = 1;
                globalVariables.valueConstants.erase(wi);
                break;
            }
        }
        if (!is_found_in_valueConstants)
            throw("tried to unfreeze global variable not from valueConstants list");

        globalVariables.mutableConstants.push_back(
            {.name = name,
             .min_value = min_value,
             .max_value = max_value,
             .init_from_config = 1,
             .init_guess = v.value,
             .parameter_position = number_parameters++, /////////////////////////////////////!!!!!!!!!!!!!!!!!!!!TODO
             .model_position = v.model_position,
             .gamma = 0, ////////////////////////////////////////////////////////////!!!!!!!!!!!TODO no gen algo for now
             .is_mutation_applicable = 1////////////////////////////////////////////////////////////!!!!!!!!!!!TODO no gen algo for now
        });
        //also add it to old params vector 
        //here we assume it is for restart of some method with unfrozen parameter
        params.push_back(v.value);
    }

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
        for (auto & s: algMap) {
            s.second.resize(s.second.find(" "));
            algebraicBiMapModel.left.insert(BiMap::left_value_type(s.second, s.first));
        }
    }

    int get_number_parameters() const
    {
        return number_parameters;
    }

    void normalize_baseline(Baseline & baseline) const
    {
        const auto [min, max] = std::minmax_element(baseline.begin(), baseline.end());
        const double vmin = *min, vmax = *max - vmin;
        for (auto & e: baseline)
            e = (e - vmin) / vmax;
    }

    void DEnormalize_baseline(Baseline & baseline) const
    {
        //TODO
        //call it from a script only
        const double vmin = -80, vmax = 34;
        for (auto & e: baseline)
            e =  e * (vmax - vmin) + vmin;
    }

    void read_baseline(Baseline & baseline, const std::string & filename)
    {
        std::ifstream file;
        file.open(filename);
        if (!file.is_open())
            throw(filename + " cannot be opened");
        double v;
        while (file >> v)
            baseline.push_back(v);
        if (is_AP_normalized)
            normalize_baseline(baseline);
    }

    void read_state()
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
        std::ifstream configFile;
        configFile.open(configFilename);
        if (!configFile.is_open())
            throw("Cannot open main config file");
        json config;
        configFile >> config;
        configFile.close();
        read_config(config);
    }
    void read_config(json & config)
    {
        //lets make it plain, no MPI-IO, TODO
       // if (mpi_rank == 0) {

            is_AP_normalized = config["is_AP_normalized"].get<int>();
            
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
                     .init_from_config = (v.find("init") != v.end()),
                     .init_guess = (v.find("init") != v.end()) ? v["init"].get<double>() : 0,
                     .parameter_position = number_parameters++,
                     .model_position = model_position,
                     .gamma = v["gamma"].get<double>(),
                     .is_mutation_applicable = (v["scale"].get<std::string>() == "linear"? 1: 2)
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
                    } catch (const char * const_isnt_a_value) {
                        throw(const_isnt_a_value);
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
                         .is_mutation_applicable = (v["scale"].get<std::string>() == "linear"? 1: 2)
                        });
                    }
                }
                //other state mutable variables for each baseline
                //are drifting
                if (config["RESET_STATES"].get<int>() == 0) {
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
            }
            
            //we may need to generate test baselines first
            bool run_type = (config.find("mode") != config.end());
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
                    initial_guess(parameters.begin(), 0);

                    std::vector<double> vconstants(model.constants_size());
                    double * constants = vconstants.data();
                    model.set_constants(constants);//!!!!!!!
                    fill_constants_y0(parameters.begin(), constants, y0.data(), baseline_number);
                    const double period = vconstants[constantsBiMapModel.left.at("stim_period")];

                    //now we need to resize apRecord according to period lenght and t_sampling
                    double num_rec = 1 + std::ceil( period / t_sampling );
                    assert(num_rec < std::numeric_limits<unsigned int>::max());
                    Baseline apRecord(static_cast<unsigned int> (num_rec));
                    model_eval(y0, model, beats_num, period, apRecord);

                    bool if_nan_in_AP = (apRecord.end() != std::find_if(apRecord.begin(), apRecord.end(), [](double v) { 
                                           return std::isnan(v); }));
                    if (if_nan_in_AP) {
                        std::cout << "NaN found in the generated baseline" << std::endl;
                        throw 1;
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
            
            //wait for a root node to complete baseline generation
            MPI_Barrier(MPI_COMM_WORLD);

            //read baselines
            apbaselines = BaselineVec();
            for (auto baseline: config["baselines"].items()) {
                auto b = baseline.value();
                std::string apfilename = b["filename_phenotype"].get<std::string>();
                apbaselines.emplace_back();
                read_baseline(apbaselines.back(), apfilename);

                //initial state may not be provided
                bool initial_state = (b.find("filename_state") != b.end());
                if (initial_state) {
                    std::string statefilename = b["filename_state"].get<std::string>();

                    //TODO
                    //read_state(..., statefilename);
                }
            }
       //     std::cout << "Unknowns: " << number_parameters << std::endl;
      //  }
        //broadcast all this nonsense;
        //TODO
    }
    std::vector<double> get_gamma_vector() const
    {
        /*
         * call it only after config read!
         * zero for a gene which does not mutate
         */
        std::vector<double> v_gamma(number_parameters);
        for (const Mutable & gl: globalVariables.mutableConstants) {
            v_gamma[gl.parameter_position] = gl.gamma;
        }
        for (const Mutable & gl: globalVariables.mutableStates) {
            v_gamma[gl.parameter_position] = gl.gamma;
        }
        for (const Variables & vars: baselineVariables) {
            for (const Mutable & gl: vars.mutableConstants) {
                v_gamma[gl.parameter_position] = gl.gamma;
            }
            for (const Mutable & gl: vars.mutableStates) {
                v_gamma[gl.parameter_position] = gl.gamma;
            }
        }
        return v_gamma;
    }
    template <typename T1, typename T2>
    int get_boundaries(T1 & Optpmin, T1 & Optpmax, T2 & is_mutation_applicable) const
    {
        std::vector<double> pmin(number_parameters), pmax(number_parameters),
                            optpmin(number_parameters), optpmax(number_parameters);

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

        //find borders for optimizer
        model_optimizer_scale(pmin.begin(), optpmin.begin());
        model_optimizer_scale(pmax.begin(), optpmax.begin());

        for (int i = 0; i < number_parameters; i++) {
            Optpmin[i] = optpmin[i];
            Optpmax[i] = optpmax[i];
        }
        return 0;
    }
    template <typename It>
    int initial_guess(It parameters_begin, const bool optimizer_called = 1) const
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
            
            
            std::vector<double> parameters(number_parameters);
            //then, set mutable global parameters
            for (const Mutable & m: globalVariables.mutableConstants) {
                parameters[m.parameter_position] = vconstants[m.model_position];
                if (m.init_from_config != 0) {
                    std::cout << m.init_guess << std::endl;
                    parameters[m.parameter_position] = m.init_guess;
                }
            }
            
            //finally, set mutable baseline parameters
            for (const Variables & vars: baselineVariables) {
                for (const Mutable & m: vars.mutableStates) {
                    parameters[m.parameter_position] = y0[m.model_position];
                }
            }
            
            if (optimizer_called) {
                model_optimizer_scale(parameters.begin(), parameters_begin);
            } else {
                std::copy(parameters.begin(), parameters.end(), parameters_begin);
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
protected:
    void direct_problem(std::vector<double> & y0, Model & model, double start_record, double tout, Table & table) const
    {
        int is_correct = 0;
        const double t0 = 0;

        try {
            solver.solve(model, y0, is_correct, t0, start_record, tout, table.get_states(), table.get_algebraic());
        } catch(...) {
            throw("Solver failed");
        }

        if (!is_correct)
            throw("Solver failed");
    }
public:
    void run_direct_and_dump(double start_record_time, double time, double dump_period, std::string filename,
        const std::vector<std::string> & dump_vars)
    {
        Model model(glob_model);

        std::vector<double> y0(model.state_size());
        std::vector<double> parameters(number_parameters);
        initial_guess(parameters.begin(), 0);
        std::vector<double> vconstants(model.constants_size());
        double * constants = vconstants.data();
        model.set_constants(constants);
        fill_constants_y0(parameters.begin(), constants, y0.data(), 0);

        std::vector<int> states_model_indices;
        std::vector<int> alg_model_indices;
        std::vector<std::string> states_names;
        std::vector<std::string> alg_names;
        //add dump_vars to a table
        for (const auto & svar: dump_vars) {
            int model_position;
            try {
                model_position = statesBiMapModel.left.at(svar);
            } catch (...) {
                continue;
            }
            std::cout << svar << std::endl;
            states_model_indices.push_back(model_position);
            states_names.push_back(svar);
        }
        for (const auto & svar: dump_vars) {
            int model_position;
            try {
                model_position = algebraicBiMapModel.left.at(svar);
            } catch (...) {
                continue;
            }
            std::cout << svar << std::endl;
            alg_model_indices.push_back(model_position);
            alg_names.push_back(svar);
        }
        
        

        Table table(1 + std::ceil((time - start_record_time) / dump_period),
                    states_model_indices.size() + alg_model_indices.size(),
                    states_model_indices, alg_model_indices,
                    states_names, alg_names);
        std::cout << "run model" << std::endl;
        direct_problem(y0, model, start_record_time, time, table);
        std::cout << "dump" << std::endl;
        table.export_csv(filename);
        std::cout << "dump complete" << std::endl;
    }
protected:
    void model_eval(std::vector<double> & y0, Model & model, int n_beats, double period, Baseline & apRecord) const
    {
        //TODO, this is temporary
        std::vector<int> states_model_indices = {0};
        std::vector<int> alg_model_indices;
        std::vector<std::string> states_names = {"V"};
        std::vector<std::string> alg_names;
        Table table(apRecord.size(), 1, states_model_indices, alg_model_indices,
                    states_names, alg_names);
        
        double start_record = period * (n_beats - 1);
        double tout = period * n_beats;
        direct_problem(y0, model, start_record, tout, table);

        //TODO
        //save to apRecord
        for (int i = 0; i < apRecord.size(); i++)
            apRecord[i] = table(i, 0);
    }
public:
    template <typename It>
    BaselineVec genetic_algorithm_calls_general(It optimizer_parameters_begin, int n_beats = -1) const
    {
        std::vector<double> model_scaled_parameters(number_parameters);
        optimizer_model_scale(optimizer_parameters_begin, model_scaled_parameters.begin());
        
        if (n_beats == -1) n_beats = beats;
        BaselineVec apmodel;
        for (const auto & v: apbaselines)
            apmodel.push_back(Baseline(v.size()));
        
        Model model(glob_model);
        bool is_solver_failed = 0;
        for (size_t i = 0; i < apbaselines.size(); i++) {
            std::vector<double> y0(model.state_size());
            Baseline & apRecord = apmodel[i];

            std::vector<double> vconstants(model.constants_size());
            double * constants = vconstants.data();
            model.set_constants(constants);//!!!!!!!
            fill_constants_y0(model_scaled_parameters.begin(), constants, y0.data(), i);
            const double period = vconstants[constantsBiMapModel.left.at("stim_period")];
            
            try {
                model_eval(y0, model, n_beats, period, apRecord);
            } catch (...) {
                //solver failed
                is_solver_failed = 1;
            }
            if (is_solver_failed)
                break;
            //save mutable variables from the state
            //since we let them to drift
            for (const Mutable & m: baselineVariables[i].mutableStates) {
                model_scaled_parameters[m.parameter_position] =  y0[m.model_position];
            }
        }
        if (is_solver_failed) {
            //i dk ugly design
            for (auto & b: apmodel)
                for (auto & e: b)
                    e = 1e6;
        } else {
            //hope its fine and no side effects
            model_optimizer_scale(model_scaled_parameters.begin(), optimizer_parameters_begin);
        }
        return apmodel;
    }

    template <typename It>
    double genetic_algorithm_calls(It parameters_begin) const
    {
        return obj.dist(apbaselines, genetic_algorithm_calls_general(parameters_begin));
    }
    template <typename V>
    void genetic_algorithm_result(const V & parameters)
    {
        std::vector<double> model_scaled_parameters(number_parameters);
        optimizer_model_scale(parameters.begin(), model_scaled_parameters.begin());
        //mirror the stage of initialization before solver.solve call
        //but just save the final result
        results = Results();
        relative_results = Results();
        for (size_t i = 0; i < apbaselines.size(); i++) {
            BaselineResult res, relative_res;
            std::vector<double> y0(glob_model.state_size());
            std::vector<double> vconstants(glob_model.constants_size());
            fill_constants_y0(model_scaled_parameters.begin(), vconstants.data(), y0.data(), i);

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
        std::copy(parameters.begin(), parameters.begin() + number_parameters, results_optimizer_format.begin());
    }
    template <typename It>
    void dump_ap(It parameters_begin, int i) const
    {
        const auto tmp_b = genetic_algorithm_calls_general(parameters_begin);
        const double dist = obj.dist(apbaselines, tmp_b);
        std::cout << "final dist: " << dist << std::endl;
        for (const auto &bs : tmp_b)
            write_baseline(bs, std::string("ap") + std::to_string(i++) + ".txt");
    }
};




template <typename Model, typename Solver, typename Objective>
class ODEoptimizationTrackVersion:
    public ODEoptimization<Model, Solver, Objective>
{
protected:
    using Base = ODEoptimization<Model, Solver, Objective>;
    using typename Base::Baseline;
    using typename Base::BaselineVec;
    using Base::apbaselines;
    using Base::obj;
    using Base::write_baseline;
    using Base::genetic_algorithm_calls_general;

    BaselineVec initial_guess_baseline;
    BaselineVec intermediate_baseline;
    
    double alpha;
public:
    using Base::is_AP_normalized;
    using Base::normalize_baseline;
    ODEoptimizationTrackVersion(const Model & model, const Solver & solver, const Objective & obj)
    : Base(model, solver, obj)
    {}
    
    template <typename It>
    void start_track(It parameters_begin)
    {
        initial_guess_baseline = genetic_algorithm_calls_general(parameters_begin, 300);
        write_baseline(initial_guess_baseline[0], "baseline_start.txt");        
        intermediate_baseline = initial_guess_baseline; //to set size and check correctness of lsm

        if (is_AP_normalized) {
            for (auto & b: initial_guess_baseline)
                normalize_baseline(b);
        }
        std::cout << "TEST: " << obj.dist(initial_guess_baseline, intermediate_baseline) << std::endl;
        set_alpha(0);
    }
    void set_alpha(double alpha_)
    {
        //alpha = 0 then initial_guess_baseline
        //alpha = 1 apbaselines
        alpha = std::pow(alpha_, 1);
        for (size_t i = 0; i < apbaselines.size(); i++) {
            for (size_t j = 0; j < apbaselines[i].size(); j++)
                intermediate_baseline[i][j] = alpha * apbaselines[i][j]
                            + (1 - alpha) * initial_guess_baseline[i][j];
            if (is_AP_normalized)
                normalize_baseline(intermediate_baseline[i]);
        }
    }
    template <typename It>
    double genetic_algorithm_calls(It parameters_begin) const
    {
     //   write_baseline(intermediate_baseline[0], "baseline_1.txt");///////////////////////////////
       // for (int i = 2; i < 10; i++) {///////////////////////////////////////////////////////////////////////
            const auto tmp_b = genetic_algorithm_calls_general(parameters_begin);
        //    write_baseline(tmp_b[0], std::string("baseline_") + std::to_string(i) + ".txt");/////////////////////////////////////////
       // }/////////////////////////////////////////////////////////
      //  auto tmp_b = intermediate_baseline;/////////////////////////////////////////////////////////////////////
        return obj.dist(intermediate_baseline, tmp_b);
    }
    
    /*
    template <typename It>
    void dump_ap(It parameters_begin, int i) const
    {
        const auto tmp_b = genetic_algorithm_calls_general(parameters_begin);
        write_baseline(tmp_b[0], std::string("ap") + std::to_string(i) + ".txt");
    }
    */
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
    template <typename T1, typename T2>
    int get_boundaries(T1 & pmin, T1 & pmax, T2 & is_mutation_applicable) const
    {
        auto ppmin = func.x_min();
        std::copy(ppmin.begin(), ppmin.end(), pmin.begin());
        auto ppmax = func.x_max();
        std::copy(ppmax.begin(), ppmax.end(), pmax.begin());
        std::fill(is_mutation_applicable.begin(), is_mutation_applicable.end(), 1);
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
    template <typename V>
    void genetic_algorithm_result(const V & parameters)
    {
        for (int i = 0; i != func.get_xdim(); i++)
            result[i] = parameters[i];
    }
    std::vector<double> get_gamma_vector() const
    {
        /*
         * call it only after config read!
         * zero for a gene which does not mutate
         */
        std::vector<double> v_gamma(get_number_parameters(), 1);
        return v_gamma;
    }
};

#endif
