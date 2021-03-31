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
#include <Eigen/StdVector>

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
protected:
    std::vector<std::string> header;
public:
    Table(int m, int n, const std::vector<std::string> & header);
    void export_csv(const std::string & filename);
};

class TableSplit
: public Table
{
    using Base = Table;
    int states_cols, alg_cols;
    int * model_indices;
public:
    TableSplit(int m, int n, std::vector<int> states_model_indices,
            std::vector<int> alg_model_indices,
            const std::vector<std::string> & states_names,
            const std::vector<std::string> & alg_names);
    ~TableSplit();
    BlockOfTable get_algebraic();
    BlockOfTable get_states();
};


int halfheight_index(const std::vector<double> & v);

template<int power = 2>
class MinimizeAPbaselines
{
public:
    using Baseline = std::vector<double>;
    using ListOfBaselines = std::vector<Baseline>;

    double distBaselines(const Baseline & a, const Baseline & b) const
    {
        double res = 0;
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); i++)
            res += std::pow(std::abs(a[i] - b[i]), power);
        res /= (double) a.size();
        return res;
    }

    double dist(const ListOfBaselines & a, const ListOfBaselines & b) const
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
    using ListOfBaselines = std::vector<Baseline>;
    
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
        const double max_rest_potential = -55;   // mV
        
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

    double dist(const ListOfBaselines & a, const ListOfBaselines & b) const
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
class ScaleMinimizeAPbaselines
{
public:
    using Baseline = std::vector<double>;
    using ListOfBaselines = std::vector<Baseline>;
    
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
    
    double LSdistBaselines(const Baseline & g, const Baseline & f) const
    {
        //watch for the correct order: g -- experimental data (not from (0,1) but wider) AP, f -- model
        //solve least square problem for f(t) = A * g(t)
        assert(g.size() == f.size());

        const double a = sum_of_elements<2>(g);
        const double p = inner_product(f, g);

        const double min_A = 1, max_A = 1.2;
        const double A = std::min(max_A, std::max(min_A, p / a));

        double res = 0;
        for (size_t i = 0; i < f.size(); i++)
            res += std::pow(std::abs(f[i] - A * g[i]), power);

        res /= (double) f.size();
        return res;
    }

    double dist(const ListOfBaselines & a, const ListOfBaselines & b) const
    {
        //watch for the correct order: a -- experimental optical mapping AP, b -- model
        assert(a.size() == b.size());
        double res = 0;
        for (size_t i = 0; i < a.size(); i++) {
            res += LSdistBaselines(a[i], b[i]);
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
    using ListOfBaselines = std::vector<Baseline>;

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

    double dist(const ListOfBaselines & a, const ListOfBaselines & b) const
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
    Solver solver;//!!!!
    Objective obj;//!!!!
    using Baseline = typename Objective::Baseline;
    using ListOfBaselines = typename Objective::ListOfBaselines;

    ListOfBaselines apbaselines;
    std::vector<int> apbaselines_halfheights;

    // There are two types of values:
    // 1. Unknowns
    // 2. Knowns

    // Each of the model parameters or model state variable
    // can be either unknown or known

    // Unknown parameters / variables
    struct Unknown
    {
        std::string name;
        std::string unique_name;
        double min_value;
        double max_value;
        int init_guess_available = 0;
        double init_guess;
        int optimizer_position;
        int model_position;
        double gamma;
        int is_mutation_applicable;
    };

    struct Known
    {
        std::string name;
        std::string unique_name;
        double value;
        int model_position;
    };
    
    // we group different values together
    struct Values
    {
        std::string groupName;
        std::vector<Unknown> unknownConstants;
        std::vector<Unknown> unknownStates;
        std::vector<Known>   knownConstants;
        std::vector<Known>   knownStates;
    };

    // global section of config
    Values globalValues;

    // baselines section of config
    // apbaselines[i] corresponds to baselineValues[i]
    std::vector<Values> baselineValues;

    //number of unknown values passed to an optimization algorithm
    int number_unknowns = 0;
    //pointers are evil!!!
    std::vector<Unknown *> pointers_unknowns;

    int mpi_rank, mpi_size;

    bool is_AP_normalized = 0; //normalize every baseline vector independently!!!
public:
int beats;
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
        for (auto pu : pointers_unknowns) {
            const auto & u = *pu;
            const int pos = u.optimizer_position;
            if (u.is_mutation_applicable == 0) {
                model_param_start[pos] = optim_param_start[pos];
            } else if (u.is_mutation_applicable == 1) {
                model_param_start[pos] = lin_scale_back(optim_param_start[pos], u.min_value, u.max_value);
            } else if (u.is_mutation_applicable == 2) {
                model_param_start[pos] = log_scale_back(optim_param_start[pos], u.min_value, u.max_value);
            } else {
                throw(u.name + ": unknown scale type");
            }
        }
    }

    template <typename It, typename It2>
    void model_optimizer_scale(It model_param_start, It2 optim_param_start) const
    {
        //from model to optimizer
        for (auto pu : pointers_unknowns) {
            const auto & u = *pu;
            const int pos = u.optimizer_position;
            if (u.is_mutation_applicable == 0) {
                optim_param_start[pos] = model_param_start[pos];
            } else if (u.is_mutation_applicable == 1) {
                optim_param_start[pos] = lin_scale(model_param_start[pos], u.min_value, u.max_value);
            } else if (u.is_mutation_applicable == 2) {
                optim_param_start[pos] = log_scale(model_param_start[pos], u.min_value, u.max_value);
            } else {
                throw(u.name + ": unknown scale type");
            }
        }
    }
public:
    void unfreeze_global_variable(const std::string & name, double min_value, double max_value, std::vector<double> & params)
    {//TODO a lot! Use it at your own risk!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        throw("unfreeze_global_variable not ready yet");
        globalValues.knownConstants;//remove from knownConstants
        bool is_found_in_knownConstants = 0;
        Known v;
        for (auto wi = globalValues.knownConstants.begin();
            wi != globalValues.knownConstants.end(); wi++) {
            if (wi->name == name) {
                v = *wi;
                is_found_in_knownConstants = 1;
                globalValues.knownConstants.erase(wi);
                break;
            }
        }
        if (!is_found_in_knownConstants)
            throw("tried to unfreeze global variable not from knownConstants list");

        globalValues.unknownConstants.push_back(
            {.name = name,
             .unique_name = name + "_global",
             .min_value = min_value,
             .max_value = max_value,
             .init_guess_available = 1,
             .init_guess = v.value,
             .optimizer_position = number_unknowns++, /////////////////////////////////////!!!!!!!!!!!!!!!!!!!!TODO
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
    : solver(solver), obj(obj)
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
        return number_unknowns;
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
            reg_alpha = config["regularization_alpha"].get<double>();
            number_unknowns = 0;
            //read global values
            globalValues.groupName = "global";
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
                    globalValues.knownConstants.push_back(
                    {.name = name,
                     .unique_name = name + "_global",
                     .value = v["value"].get<double>(),
                     .model_position = model_position
                    });
                } else {
                    //unknown constant
                    globalValues.unknownConstants.push_back(
                    {.name = name,
                     .unique_name = name + "_global",
                     .min_value = v["bounds"][0].get<double>(),
                     .max_value = v["bounds"][1].get<double>(),
                     .init_guess_available = (v.find("init") != v.end()),
                     .init_guess = (v.find("init") != v.end()) ? v["init"].get<double>() : 0,
                     .optimizer_position = number_unknowns++,
                     .model_position = model_position,
                     .gamma = v["gamma"].get<double>(),
                     .is_mutation_applicable = (v["scale"].get<std::string>() == "linear"? 1: 2)
                    });
                }
            }
            //constants not listed in config will have default values


            //now lets read baseline values
            int baselines_num = 0;
            for (auto baseline: config["baselines"].items()) {
                baselines_num++;
                auto b = baseline.value();

                baselineValues.emplace_back();
                Values & bVar = baselineValues.back();
                BiMap statesBiMapDrifting = statesBiMapModel;
                bVar.groupName = b["name"].get<std::string>();
                for (auto variable: b["params"].items()) {
                    auto v = variable.value();
                    std::string name = v["name"].get<std::string>();
                    int model_position;
                    try {
                        //constant can be listed as a baseline value
                        model_position = constantsBiMapModel.left.at(name);
                        bool is_value = (v.find("value") != v.end());
                        if (is_value) {
                            bVar.knownConstants.push_back(
                            {.name = name,
                             .unique_name = name + "_" + bVar.groupName,
                             .value = v["value"].get<double>(),
                             .model_position = model_position
                            });
                        } else {
                            //unknown constant specified for this baseline
                            //i.e. stimulation shift
                            bVar.unknownConstants.push_back(
                            {.name = name,
                             .unique_name = name + "_" + bVar.groupName,
                             .min_value = v["bounds"][0].get<double>(),
                             .max_value = v["bounds"][1].get<double>(),
                             .init_guess_available = (v.find("init") != v.end()),
                             .init_guess = (v.find("init") != v.end()) ? v["init"].get<double>() : 0,
                             .optimizer_position = number_unknowns++,
                             .model_position = model_position,
                             .gamma = v["gamma"].get<double>(),
                             .is_mutation_applicable = (v["scale"].get<std::string>() == "linear"? 1: 2)
                            });
                        }
                    } catch (...) {
                        //name is not a constant
                        try {
                            model_position = statesBiMapModel.left.at(name);
                        } catch (...) {
                            //or it is neither a constant nor a state variable
                            throw(name + " in baseline config is neither constant nor state variable");
                        }
                        //remove it from statesBiMapDrifting
                        statesBiMapDrifting.left.erase(name);
                        
                        bool is_value = (v.find("value") != v.end());
                        if (is_value) {
                            bVar.knownStates.push_back(
                            {.name = name,
                             .unique_name = name + "_" + bVar.groupName,
                             .value = v["value"].get<double>(),
                             .model_position = model_position
                            });
                        } else {
                            bVar.unknownStates.push_back(
                            {.name = name,
                             .unique_name = name + "_" + bVar.groupName,
                             .min_value = v["bounds"][0].get<double>(),
                             .max_value = v["bounds"][1].get<double>(),
                             .optimizer_position = number_unknowns++,
                             .model_position = model_position,
                             .gamma = v["gamma"].get<double>(),
                             .is_mutation_applicable = (v["scale"].get<std::string>() == "linear"? 1: 2)
                            });
                        }
                    }
                }
                // statesBiMapDrifting may not be empty
                // if RESET_STATES == 1 then statesBiMapDrifting states
                // are reset to default before model runs
                //
                // if RESET_STATES == 0 then statesBiMapDrifting states
                // drift (values are overwritten by model)
                if (config["RESET_STATES"].get<int>() == 0) {
                    for (const auto & sit: statesBiMapDrifting) {
                        bVar.unknownStates.push_back(
                        {.name = sit.left,
                         .unique_name = sit.left + "_" + bVar.groupName,
                         .min_value = 0,
                         .max_value = 0,
                         .optimizer_position = number_unknowns++,
                         .model_position = sit.right,
                         .gamma = 0,
                         .is_mutation_applicable = 0
                        });
                    }
                }
            }
            
            //fill pointers_unknowns!!!!!!!!!!!!!!!!!!!!!
            pointers_unknowns.resize(number_unknowns);
            for (Unknown & gl: globalValues.unknownConstants) {
                pointers_unknowns[gl.optimizer_position] = &gl;
            }
            for (Unknown & gl: globalValues.unknownStates) {
                pointers_unknowns[gl.optimizer_position] = &gl;
            }

            for (Values & vars: baselineValues) {
                for (Unknown & gl: vars.unknownConstants) {
                    pointers_unknowns[gl.optimizer_position] = &gl;
                }
                for (Unknown & gl: vars.unknownStates) {
                    pointers_unknowns[gl.optimizer_position] = &gl;
                }
            }
            if (mpi_rank == 0) {
                for (auto pu : pointers_unknowns) {
                    const auto & u = *pu;
                    std::cout << u.optimizer_position << " " << u.unique_name << std::endl;
                }
            }
/*obsolete since we can call direct problem directly
            //we may need to generate test baselines first
            bool run_type = (config.find("mode") != config.end());
            if (mpi_rank == 0 && run_type && config["mode"].get<std::string>() == "Generate synthetic baselines") {

                std::cout << "Generating synthetic baselines" << std::endl;
                std::cout << "ATTENTION! ALL BASELINES FROM CONFIG WILL BE OVERWRITTEN" << std::endl;
                const double t_sampling = config["t_sampling"].get<double>();
                const int beats_num = config["generator_beats"].get<int>();//run model long enough to get a steady state

                int baseline_number = 0;

                for (auto baseline: config["baselines"].items()) {
                    Model model;
                    std::vector<double> y0(model.state_size());
                    
                    std::vector<double> parameters(number_unknowns);
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
*/

      //  }
        //broadcast all this nonsense;
        //TODO
    }
    void read_baselines(json & config)
    {
        apbaselines = ListOfBaselines();
        for (auto baseline: config["baselines"].items()) {
            auto b = baseline.value();
            std::string apfilename = b["filename_phenotype"].get<std::string>();
            apbaselines.emplace_back();
            read_baseline(apbaselines.back(), apfilename);

            apbaselines_halfheights.emplace_back();
            const int hh_index = halfheight_index(apbaselines.back());
            if (hh_index == -1)
                throw("Fix baselines");
            apbaselines_halfheights.back() = hh_index;
            std::cout << "hh_index: " << hh_index << std::endl;
            //initial state may not be provided
            bool initial_state = (b.find("filename_state") != b.end());
            if (initial_state) {
                std::string statefilename = b["filename_state"].get<std::string>();

                //TODO
                //read_state(..., statefilename);
            }
        }
    }
    std::vector<double> get_gamma_vector() const
    {
        /*
         * call it only after config read!
         * zero for a gene which does not mutate
         */
        std::vector<double> v_gamma(number_unknowns);
        for (auto pu : pointers_unknowns) {
            const auto & u = *pu;
            v_gamma[u.optimizer_position] = u.gamma;
        }
        return v_gamma;
    }


    template <typename T1, typename T2>
    int get_boundaries_updated(T1 & Optpmin, T1 & Optpmax, T2 & is_mutation_applicable) const
    {
        const double eps = 0.1;
        for (const Unknown & gl: globalValues.unknownConstants) {
            Optpmin[gl.optimizer_position] = (1 - eps) * results_optimizer_format[gl.optimizer_position];
            Optpmax[gl.optimizer_position] = (1 + eps) * results_optimizer_format[gl.optimizer_position];
            is_mutation_applicable[gl.optimizer_position] = 1;//gl.is_mutation_applicable; //no need to know that it was log scaled
        }
        for (const Unknown & gl: globalValues.unknownStates) {
            Optpmin[gl.optimizer_position] = (1 - eps) * results_optimizer_format[gl.optimizer_position];
            Optpmax[gl.optimizer_position] = (1 + eps) * results_optimizer_format[gl.optimizer_position];
            is_mutation_applicable[gl.optimizer_position] = 1;//gl.is_mutation_applicable;
        }
        for (const Values & vars: baselineValues) {
            for (const Unknown & gl: vars.unknownConstants) {
                Optpmin[gl.optimizer_position] = (1 - eps) * results_optimizer_format[gl.optimizer_position];
                Optpmax[gl.optimizer_position] = (1 + eps) * results_optimizer_format[gl.optimizer_position];
                is_mutation_applicable[gl.optimizer_position] = 1;//gl.is_mutation_applicable;
            }
            for (const Unknown & gl: vars.unknownStates) {
                Optpmin[gl.optimizer_position] = (1 - eps) * results_optimizer_format[gl.optimizer_position];
                Optpmax[gl.optimizer_position] = (1 + eps) * results_optimizer_format[gl.optimizer_position];
                is_mutation_applicable[gl.optimizer_position] = 1;//gl.is_mutation_applicable;
                //TODO
                //tricky place
                //we have two types of unknownStates
                //for drifting states we probably don't have max and min
            }
        }
        return 0;
    }

    template <typename T1, typename T2>
    int get_boundaries(T1 & Optpmin, T1 & Optpmax, T2 & is_mutation_applicable) const
    {
        //why we pass is_mutation_applicable = 1 ?
        //because some optimization methods still assume they
        //need to scale unknowns on their own
        std::vector<double> pmin(number_unknowns), pmax(number_unknowns);

        for (auto pu : pointers_unknowns) {
            const auto & u = *pu;
            pmin[u.optimizer_position] = u.min_value;
            pmax[u.optimizer_position] = u.max_value;
            if (u.is_mutation_applicable == 0)
                is_mutation_applicable[u.optimizer_position] = 0;
            else //no need to know that it was log scaled
                is_mutation_applicable[u.optimizer_position] = 1;
            
        }

        //find borders for optimizer
        model_optimizer_scale(pmin.begin(), Optpmin.begin());
        model_optimizer_scale(pmax.begin(), Optpmax.begin());

        return 0;
    }

protected:   
    void get_default_values(double * constants, double * y0) const
    {
        Model::initConsts(constants);
        Model::initState(y0);
    }


    template <typename It>
    void initial_guess(It parameters_begin) const
    {
        std::vector<double> y0(Model::state_size());
        std::vector<double> vconstants(Model::constants_size());
        
        //first, set model's default values
        get_default_values(vconstants.data(), y0.data()); 

        for (const Unknown & m: globalValues.unknownConstants) {
            parameters_begin[m.optimizer_position] = vconstants[m.model_position];
        }
        for (const Values & vars: baselineValues) {
            for (const Unknown & m: vars.unknownStates) {
                parameters_begin[m.optimizer_position] = y0[m.model_position];
            }
            for (const Unknown & m: vars.unknownConstants) {
                parameters_begin[m.optimizer_position] = vconstants[m.model_position];
            }
        }

        //then, set unknowns initial guesses if available
        for (auto pu : pointers_unknowns) {
            const auto & u = *pu;
            if (u.init_guess_available) {
                parameters_begin[u.optimizer_position] = u.init_guess;
            }
        }
    }
public:
    template <typename It>
    int initial_guess_for_optimizer(It parameters_begin) const
    {
        if (1) {
            std::vector<double> parameters(number_unknowns);
            initial_guess(parameters.begin());
            model_optimizer_scale(parameters.begin(), parameters_begin);
            return 0;
        } else {
            //or just not provide anything
            return -1;
        }
    }
protected:
    template <typename It>
    void fill_constants_y0(It parameters_begin, double * constants, double * y0, size_t baseline_index) const
    {
        //first, default
        get_default_values(constants, y0);

        //global section of config overwrites default
        for (const Unknown & m: globalValues.unknownConstants) {
            constants[m.model_position] = parameters_begin[m.optimizer_position];
        }
        for (const Unknown & m: globalValues.unknownStates) {
            throw("Wait, that's illegal. Please contact us if you want it.");
        }
        for (const Known & v: globalValues.knownConstants) {
            constants[v.model_position] = v.value;
        }
        for (const Known & v: globalValues.knownStates) {
            throw("Wait, that's illegal. Please contact us if you want it.");
        }

        //baseline sections of config overwrites global and default
        for (const Unknown & m: baselineValues[baseline_index].unknownConstants) {
            constants[m.model_position] = parameters_begin[m.optimizer_position];
        }
        for (const Unknown & m: baselineValues[baseline_index].unknownStates) {
            y0[m.model_position] = parameters_begin[m.optimizer_position];
        }
        for (const Known & v: baselineValues[baseline_index].knownConstants) {
            constants[v.model_position] = v.value;
        }
        for (const Known & v: baselineValues[baseline_index].knownStates) {
            y0[v.model_position] = v.value;
        }
    }

    void direct_problem(std::vector<double> & y0, Model & model, double start_record, double tout, TableSplit & table) const
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
        for (int baseline_index = 0; baseline_index < baselineValues.size(); baseline_index++) {
            
            const auto & baseline = baselineValues[baseline_index];
            Model model;

            std::vector<double> y0(model.state_size());
            
            std::vector<double> parameters(number_unknowns);
            //Initial guess from config file will be respected!
            //Maybe user would like to verify that initial guess is a correct model
            initial_guess(parameters.begin());
            
            std::vector<double> vconstants(model.constants_size());
            double * constants = vconstants.data();
            model.set_constants(constants);
            
            fill_constants_y0(parameters.begin(), constants, y0.data(), baseline_index);

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


            TableSplit table(1 + std::ceil((time - start_record_time) / dump_period),
                        states_model_indices.size() + alg_model_indices.size(),
                        states_model_indices, alg_model_indices,
                        states_names, alg_names);
            direct_problem(y0, model, start_record_time, time, table);
            table.export_csv(filename + "_" + baseline.groupName);
        }
    }
protected:
    void model_eval(std::vector<double> & y0, Model & model, int n_beats, double period, Baseline & apRecord) const
    {
        // called while solving optimization problem
        // when we try to fit AP curves
        std::vector<int> states_model_indices = {0};
        std::vector<int> alg_model_indices;
        std::vector<std::string> states_names = {"V"};
        std::vector<std::string> alg_names;
        TableSplit table(apRecord.size(), 1, states_model_indices, alg_model_indices,
                    states_names, alg_names);
        
        double start_record = period * (n_beats - 1);
        double tout = period * n_beats;
        direct_problem(y0, model, start_record, tout, table);

        //TODO
        //save to apRecord
        for (int i = 0; i < apRecord.size(); i++)
            apRecord[i] = table(i, 0);
    }
    double reg_alpha;
public:

    double parameter_penalty(std::vector<double> & parameters) const
    {
        double penalty = 0;
        for (auto pu : pointers_unknowns) {
            const auto & u = *pu;
            const int pos = u.optimizer_position;
            const double val = parameters[pos];
            if (!u.is_mutation_applicable)
                continue;
            if (u.min_value > val)
                penalty += std::pow(20 * (val - u.min_value) / u.min_value, 2);///////////// is div by zero possible? TODO
            if (u.max_value < val)
                penalty += std::pow(20 * (val - u.max_value) / u.max_value, 2);///////////// is div by zero possible? TODO
        }
        return penalty;
    }
    template <typename It>
    ListOfBaselines genetic_algorithm_calls_general(It optimizer_parameters_begin, double & extra_penalty, int n_beats = -1) const
    {
        std::vector<double> model_scaled_parameters(number_unknowns);
        optimizer_model_scale(optimizer_parameters_begin, model_scaled_parameters.begin());

        if (n_beats == -1) n_beats = beats;
        ListOfBaselines apmodel;
        for (const auto & v: apbaselines)
            apmodel.push_back(Baseline(v.size()));

        Model model;
        bool solver_failed = 0;
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
                solver_failed = 1;
            }
            if (solver_failed)
                break;
            //save mutable variables from the state
            //since we let them to drift
            for (const Unknown & m: baselineValues[i].unknownStates) {
                model_scaled_parameters[m.optimizer_position] = y0[m.model_position];
            }
            extra_penalty = parameter_penalty(model_scaled_parameters);
            
            
            // rotate apRecord so halfheights are same
            const int indexApRecord = halfheight_index(apRecord);
            if (indexApRecord == -1) {
                // assume it is a problem with this exact guess made by optimizer
                // and do nothing
            } else {
                const int diff = indexApRecord - apbaselines_halfheights[i];
                if (diff >= 0)
                    std::rotate(apRecord.begin(), apRecord.begin() + diff, apRecord.end());
                else
                    std::rotate(apRecord.begin(), apRecord.end() + diff, apRecord.end());
            }
        }
        if (solver_failed) {
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
    double regularization(It parameters_begin) const
    {
        double reg = 0;
        std::vector<double> prior(pointers_unknowns.size());
        initial_guess_for_optimizer(prior.begin());
        for (auto pu : pointers_unknowns) {
            const auto & u = *pu;
            const int pos = u.optimizer_position;
            reg += std::pow(parameters_begin[pos] - prior[pos], 2);
        }
        return reg_alpha * reg;
    }

    template <typename It>
    double genetic_algorithm_calls(It parameters_begin) const
    {
        double boundaries_penalty;
        double main_penalty = obj.dist(apbaselines, genetic_algorithm_calls_general(parameters_begin, boundaries_penalty));
        return main_penalty + boundaries_penalty + regularization(parameters_begin);
    }
    template <typename V>
    void genetic_algorithm_result(const V & parameters)
    {
        std::vector<double> model_scaled_parameters(number_unknowns);
        optimizer_model_scale(parameters.begin(), model_scaled_parameters.begin());
        //mirror the stage of initialization before solver.solve call
        //but just save the final result
        results = Results();
        relative_results = Results();
        for (size_t i = 0; i < apbaselines.size(); i++) {
            BaselineResult res, relative_res;
            std::vector<double> y0(Model::state_size());
            std::vector<double> vconstants(Model::constants_size());
            fill_constants_y0(model_scaled_parameters.begin(), vconstants.data(), y0.data(), i);

            std::vector<double> y0_default(Model::state_size());
            std::vector<double> vconstants_default(Model::constants_size());
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
        results_optimizer_format = std::vector<double>(number_unknowns);
        std::copy(parameters.begin(), parameters.begin() + number_unknowns, results_optimizer_format.begin());
    }

    template <typename It>
    void dump_ap(It parameters_begin, int i, int num_beats = -1) const
    {
        double extra_penalty;
        const auto tmp_b = genetic_algorithm_calls_general(parameters_begin, extra_penalty, num_beats);
        const double dist = obj.dist(apbaselines, tmp_b) + extra_penalty;
        std::cout << "final dist: " << dist << std::endl;
        for (const auto &bs : tmp_b)
            write_baseline(bs, std::string("ap") + std::to_string(i++) + ".txt");
    }

    double loglikelihood(std::vector<double> pars)
    {
        double ll = 0;
        double extra_penalty;
        ListOfBaselines res = genetic_algorithm_calls_general(pars.begin(), extra_penalty); //extra_penalty is not used further
        assert(apbaselines.size() == res.size());
        const double sigma = 100;//it is just an assumption, need to find sigma from experimental data TODO
        const double add = -log(sqrt(2*M_PI) * sigma);//with normalization
        for (size_t i = 0; i < apbaselines.size(); i++) {
            const Baseline & a = apbaselines[i];
            const Baseline & r = res[i];
            double p = 0;
            //maybe i should use sizes of these vectors...
            for (size_t j = 0; j < a.size(); j++) {
                const double arg = (a[j] - r[j]) / sigma;
                p += (- 0.5 * arg * arg + add); 
            }
            ll += p ;// a.size();
        }
        return ll;
    }
    std::vector<std::string> get_unique_parameter_names()
    {
        std::vector<std::string> names(number_unknowns);
        for (int i = 0; i < number_unknowns; i++) {
            names[i] = pointers_unknowns[i]->unique_name;
        }
        /*
        for (const Unknown & gl: globalValues.unknownConstants) {
            names[gl.optimizer_position] = gl.name;
        }
        for (const Unknown & gl: globalValues.unknownStates) {
            names[gl.optimizer_position] = gl.name;
        }

        for (const Values & vars: baselineValues) {
            std::string unique_suffix = vars.groupName;
            for (const Unknown & gl: vars.unknownConstants) {
                names[gl.optimizer_position] = gl.name + unique_suffix;
            }
            for (const Unknown & gl: vars.unknownStates) {
                names[gl.optimizer_position] = gl.name + unique_suffix;
            }
        }*/
        return names;
    }

    std::vector<Table,Eigen::aligned_allocator<Table>> gen_algo_tables;

    void gen_algo_stats(const std::vector<std::pair<double, int>> & sd_n_index, const std::vector<double> & all_genes, int gen, int total_gen)
    {
        const int percentile_step = 5;
        if (gen == 0) {
            for (int i = 0; i < number_unknowns; i++) {
                std::vector<std::string> header;
                header.push_back("best_organism");
                for (int p = percentile_step; p <= 100; p += percentile_step) {
                    header.push_back(std::string("mean_") + std::to_string(p));
                    header.push_back(std::string("var_") + std::to_string(p));
                }
                gen_algo_tables.push_back(Table(total_gen, header.size(), header));
            }
        }

        const int number_organisms = all_genes.size() / number_unknowns;

        // first save best organism's parameters
        for (int i = 0; i < number_unknowns; i++) {
            std::vector<double> model_params(number_unknowns);
            optimizer_model_scale(all_genes.data() + sd_n_index[0].second * number_unknowns, model_params.begin());
            gen_algo_tables[i](gen, 0) = model_params[i];
        }

        // now write statistics by percentiles
        std::vector<double> sums_x2(number_unknowns), sums_x(number_unknowns);
        const double sz = (double) percentile_step * number_organisms / 100;

        for (int i = 0, percentile = 0; i < number_organisms; i++) {
            std::vector<double> model_params(number_unknowns);
            optimizer_model_scale(all_genes.data() + sd_n_index[i].second * number_unknowns, model_params.begin());
            for (int j = 0; j < number_unknowns; j++) {
                const double v = model_params[j];
                sums_x[j] += v;
                sums_x2[j] += std::pow(v, 2);
            }
            if ((int) std::floor((i + 1) / sz) > percentile) {
                percentile += 1;
                for (int j = 0; j < number_unknowns; j++) {
                    const double i_th_median = sums_x[j] / (i + 1);
                    const double i_th_var = sums_x2[j] / (i + 1) - std::pow(i_th_median, 2);
                    gen_algo_tables[j](gen, 1 + 2 * (percentile - 1)     ) = i_th_median;
                    gen_algo_tables[j](gen, 1 + 2 * (percentile - 1)  + 1) = i_th_var;
                }
            }
        }
    }
    void export_gen_algo_tables()
    {
        for (int i = 0; i < number_unknowns; i++) {
            gen_algo_tables[i].export_csv(std::string("table_") + pointers_unknowns[i]->unique_name);
        }
    }
};




template <typename Model, typename Solver, typename Objective>
class ODEoptimizationTrackVersion:
    public ODEoptimization<Model, Solver, Objective>
{
protected:
    using Base = ODEoptimization<Model, Solver, Objective>;
    using typename Base::Baseline;
    using typename Base::ListOfBaselines;
    using Base::apbaselines;
    using Base::obj;
    using Base::write_baseline;
    using Base::genetic_algorithm_calls_general;

    ListOfBaselines initial_guess_baseline;
    ListOfBaselines intermediate_baseline;
    
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
        double extra_penalty;
        initial_guess_baseline = genetic_algorithm_calls_general(parameters_begin, extra_penalty, 300);
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
       double extra_penalty;
            const auto tmp_b = genetic_algorithm_calls_general(parameters_begin, extra_penalty);
        //    write_baseline(tmp_b[0], std::string("baseline_") + std::to_string(i) + ".txt");/////////////////////////////////////////
       // }/////////////////////////////////////////////////////////
      //  auto tmp_b = intermediate_baseline;/////////////////////////////////////////////////////////////////////
        return obj.dist(intermediate_baseline, tmp_b) + extra_penalty;
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
