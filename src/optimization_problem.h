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

#include "stimulation.h"
#include "objective.h"

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


/**
 * @brief Finds depolarization halfheight position in @p v
 *
 * It is assumed that action potential is completely inside the single period
 * represented by @p v.
 *
 * @param[in] v Action potential recording
 *
 * @return halfheight index
 */
int halfheight_index(const std::vector<double> & v);


template <typename Model, typename Solver>
class ODEoptimization
{
protected:
    using Baseline = std::vector<double>;
    using VectorOfBaselines = std::vector<Baseline>;

    Solver solver;
    std::unique_ptr<BaseObjective<Baseline, VectorOfBaselines>> obj;

    //set_constants is not thread-safe
    Model g_model; //use this model only as a reference!!!! Do not call it directly but make a copy of it!

    VectorOfBaselines apbaselines;

    bool ignore_before_halfheight;
    std::vector<int> apbaselines_halfheights;

    std::vector<StimulationBase*> stimulation_protocols;

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

    // special section for spontaneously beating checks
    Values spontBeatValues;
    enum ExpectedSpontBeat {NoChecks, NonSpont, Spont};
    ExpectedSpontBeat expectedSpontBeat = NoChecks;
    double spontThreshold = 0;
    int spontBeatCheckSize = 0;
    double spontBeatCheckFailScaler = 1;

    template<typename It>
    bool spontBeatCheck(It optimizer_parameters_begin) const
    {
        if (expectedSpontBeat == NoChecks)
            return 1;
        std::vector<double> model_scaled_parameters(number_unknowns);
        optimizer_model_scale(optimizer_parameters_begin, model_scaled_parameters.begin());

        Baseline baseline = generate_one_baseline(model_scaled_parameters, spontBeatCheckSize, spontBeatValues, 0, 1);
        model_optimizer_scale(model_scaled_parameters.begin(), optimizer_parameters_begin);

        double maxV = *std::max_element(baseline.begin(), baseline.end());
        bool noActivation = (maxV < spontThreshold);
        return (expectedSpontBeat == NonSpont) == noActivation;
    }


    //number of unknown values passed to an optimization algorithm
    int number_unknowns = 0;
    std::vector<Unknown *> pointers_unknowns;

    int mpi_rank, mpi_size;

    bool is_AP_normalized = 0; //normalize every baseline vector independently!!!
public:
	~ODEoptimization()
	{
		for (auto s: stimulation_protocols)
			delete [] s;
	}

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


    /**
     * @brief Log scale to unit interval
     *
     * The function log scales x from [minOrig, maxOrig] into [0, 1]
     */
    double log_scale(double x, double minOrig, double maxOrig) const
    {
        if (minOrig <= 0) {
            throw("minOrig <= 0");
        } else {
            return log(x / minOrig) / log(maxOrig / minOrig);
        }
    }

    /**
     * @brief Log scale to [minOrig, maxOrig]
     *
     * The function log scales x from [0, 1] into [minOrig, maxOrig]
     *
     */
    double log_scale_back(double x, double minOrig, double maxOrig) const
    {
        if (minOrig <= 0) {
            throw("minOrig <= 0");
        } else {
            return exp(x * log(maxOrig / minOrig)) * minOrig;
        }
    }

    /**
     * @brief Linear scale to unit interval
     *
     * The function linearly scales x from [minOrig, maxOrig] into [0, 1]
     */
    double lin_scale(double x, double minOrig, double maxOrig) const
    {
        return (x - minOrig) / (maxOrig - minOrig);
    }

    /**
     * @brief Linear scale to [minOrig, maxOrig]
     *
     * The function linearly scales x from [0, 1] into [minOrig, maxOrig]
     *
     */
    double lin_scale_back(double x, double minOrig, double maxOrig) const
    {
        return (maxOrig - minOrig) * x + minOrig;
    }

    /**
     * @brief Optimizer to model parameters transform
     *
     *
     */
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

    /**
     * @brief Model to optimizer parameters transform
     */
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
    {/// @todo not ready, use it at your own risk!
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
    ODEoptimization(const Model & model, const Solver & solver)
    : solver(solver), g_model(model)
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

    void denormalize_baseline(Baseline & baseline) const
    {
        /// @todo
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
        /// @todo
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

    void read_baseline_config(json & config, json & b, Values & bVar)
    {
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

    void read_config(json & config)
    {
        const std::string sname = config["script"].get<std::string>();
        if (sname != "Direct Problem") {
            num_repeated_obj_runs = config["repeated_obj_runs"].get<int>();
            ignore_before_halfheight = config["ignore_before_halfheight"].get<int>();
            obj = new_objective<Baseline, VectorOfBaselines> (config["Objective"].get<std::string>());
            reg_alpha = config["regularization_alpha"].get<double>();
            is_AP_normalized = config["is_AP_normalized"].get<int>();
        }
        beats = config["n_beats"].get<int>();
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
        for (auto baseline: config["baselines"].items()) {
            auto b = baseline.value();
            std::string expectedSpontBeatString;
            try {
                expectedSpontBeatString = b["expectedSpontBeat"].get<std::string>();
            } catch(...) {
                expectedSpontBeatString = "";
            }
            if (expectedSpontBeatString == "") {
                // ordinary baseline
                baselineValues.emplace_back();
                Values & bVar = baselineValues.back();
                read_baseline_config(config, b, bVar);
                try {
                    if (b["stimProtocol"].get<std::string>() == "Biphasic")
                        stimulation_protocols.push_back(new BiphasicStim(b["stimAmplitude"].get<double>(), b["pcl"].get<double>(), b["stim_shift"].get<double>(), b["pulseDuration"].get<double>()));
                    else if (b["stimProtocol"].get<std::string>() == "BiphasicStim_CaSR_Protocol")
                        stimulation_protocols.push_back(new BiphasicStim_CaSR_Protocol(b["stimAmplitude"].get<double>(),
                                        b["pcl_start"].get<double>(), b["pcl_end"].get<double>(), b["growth_time"].get<double>(),
                                        b["pcl_end_duration"].get<double>(),  b["stim_shift"].get<double>(), b["pulseDuration"].get<double>()));

                }
                catch (...) {
                    stimulation_protocols.push_back(new StimulationNone());
                }
            } else {
                // spontBeatChecker
                if (expectedSpontBeatString == "NoChecks")
                    expectedSpontBeat = NoChecks;
                else if (expectedSpontBeatString == "NonSpont")
                    expectedSpontBeat = NonSpont;
                else if (expectedSpontBeatString == "Spont")
                    expectedSpontBeat = Spont;
                else
                    throw("unknown expectedSpontBeat type in config");

                spontBeatCheckFailScaler = b["FailScaler"].get<double>();
                spontThreshold = b["spontThreshold"].get<double>();
                spontBeatCheckSize = b["spontBeatCheckSize"].get<int>();
                read_baseline_config(config, b, spontBeatValues);
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

        for (Unknown & gl: spontBeatValues.unknownConstants) {
            pointers_unknowns[gl.optimizer_position] = &gl;
        }
        for (Unknown & gl: spontBeatValues.unknownStates) {
            pointers_unknowns[gl.optimizer_position] = &gl;
        }


        if (mpi_rank == 0) {
            for (auto pu : pointers_unknowns) {
                const auto & u = *pu;
                std::cout << u.optimizer_position << " " << u.unique_name << std::endl;
            }
        }
    }
    void read_baselines(json & config)
    {
        apbaselines = VectorOfBaselines();
        for (auto baseline: config["baselines"].items()) {
            auto b = baseline.value();
            std::string apfilename;
            try {
                apfilename = b["filename_phenotype"].get<std::string>();
            } catch(...) {
                continue;
            }
            apbaselines.emplace_back();
            read_baseline(apbaselines.back(), apfilename);

            apbaselines_halfheights.emplace_back();
            const int hh_index = halfheight_index(apbaselines.back());
            if (hh_index == -1)
                throw("Fix baselines");
            apbaselines_halfheights.back() = hh_index;
            //initial state may not be provided
            // do we really need it?
            bool initial_state = (b.find("filename_state") != b.end());
            if (initial_state) {
                std::string statefilename = b["filename_state"].get<std::string>();

                /// @todo
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

/*
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
*/
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
        g_model.initConsts(constants);
        g_model.initState(y0);
    }

    template <typename It>
    void initial_guess(It parameters_begin) const
    {
        std::vector<double> y0(g_model.state_size());
        std::vector<double> vconstants(g_model.constants_size());

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
        for (const Unknown & m: spontBeatValues.unknownStates) {
            parameters_begin[m.optimizer_position] = y0[m.model_position];
        }
        for (const Unknown & m: spontBeatValues.unknownConstants) {
            parameters_begin[m.optimizer_position] = vconstants[m.model_position];
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
    void fill_constants_y0(It parameters_begin, double * constants, double * y0, const Values & baselineValue) const
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
        for (const Unknown & m: baselineValue.unknownConstants) {
            constants[m.model_position] = parameters_begin[m.optimizer_position];
        }
        for (const Unknown & m: baselineValue.unknownStates) {
            y0[m.model_position] = parameters_begin[m.optimizer_position];
        }
        for (const Known & v: baselineValue.knownConstants) {
            constants[v.model_position] = v.value;
        }
        for (const Known & v: baselineValue.knownStates) {
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
        const std::vector<std::string> & dump_vars) const
    {
        for (unsigned baseline_index = 0; baseline_index < baselineValues.size(); baseline_index++) {

            const auto & baseline = baselineValues[baseline_index];
            Model model(g_model);

            std::vector<double> y0(model.state_size());

            std::vector<double> parameters(number_unknowns);
            //Initial guess from config file will be respected!
            //Maybe user would like to verify that initial guess is a correct model
            initial_guess(parameters.begin());

            std::vector<double> vconstants(model.constants_size());
            double * constants = vconstants.data();
            model.set_constants(constants);
            model.set_stimulation(stimulation_protocols[baseline_index]);
            fill_constants_y0(parameters.begin(), constants, y0.data(), baselineValues[baseline_index]);

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
    void get_voltage_trace(std::vector<double> & y0, Model & model, int n_beats, double period, Baseline & apRecord) const
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

        //save to apRecord
        for (size_t i = 0; i < apRecord.size(); i++)
            apRecord[i] = table(i, 0);
    }
    double reg_alpha;
public:

    /**
     * @brief Boundary restriction penalty in model scale
     *
     * The function provides boundary restriction penalty for parameters in model scale.
     * It is meaningful only for unconstained optimizators such as nelder-mead, gradient descent etc.
     *
     */
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
    double parameter_penalty_optimizer(It optimizer_parameters_begin) const
    {
        std::vector<double> model_scaled_parameters(number_unknowns);
        optimizer_model_scale(optimizer_parameters_begin, model_scaled_parameters.begin());
        return parameter_penalty(model_scaled_parameters);
    }


    ///@todo create generate_one_baseline function and call it from generate_baselines
    Baseline generate_one_baseline(std::vector<double> & model_scaled_parameters, int size,
            const Values & baselineValue, StimulationBase * stimulation_protocol = nullptr, int n_beats = -1) const
    {
        if (n_beats == -1) n_beats = beats;
        Model model(g_model);

        Baseline baseline(size);

        std::vector<double> y0(model.state_size());
        std::vector<double> vconstants(model.constants_size());
        double * constants = vconstants.data();
        model.set_constants(constants);

        auto sp_unique = std::make_unique<StimulationNone>();
        if (stimulation_protocol == nullptr) {
            model.set_stimulation(sp_unique.get());
        } else {
            model.set_stimulation(stimulation_protocol);
        }
        fill_constants_y0(model_scaled_parameters.begin(), constants, y0.data(), baselineValue);
        const double period = vconstants[constantsBiMapModel.left.at("stim_period")];

        try {
            get_voltage_trace(y0, model, n_beats, period, baseline);
            //save mutable variables from the state
            //since we let them drift
            for (const Unknown & m: baselineValue.unknownStates) {
                model_scaled_parameters[m.optimizer_position] = y0[m.model_position];
            }
        } catch (...) {
            //solver failed
            for(auto &e: baseline)
                e = 1e6;
        }
        return baseline;
    }
    template <typename It>
    VectorOfBaselines generate_baselines(It optimizer_parameters_begin, int n_beats = -1) const
    {
        std::vector<double> model_scaled_parameters(number_unknowns);
        optimizer_model_scale(optimizer_parameters_begin, model_scaled_parameters.begin());

        VectorOfBaselines gen_baselines;

        for (size_t i = 0; i < apbaselines.size(); i++) {
            gen_baselines.push_back(generate_one_baseline(model_scaled_parameters, apbaselines[i].size(),
                        baselineValues[i], stimulation_protocols[i], n_beats));

            Baseline & apRecord = gen_baselines.back();
            // roll apRecord so halfheights are aligned
            const int indexApRecord = halfheight_index(apRecord);
            if (indexApRecord == -1) {
                // assume it is a problem with this exact guess made by optimizer
                // and do nothing
            } else {
                const int diff = indexApRecord - apbaselines_halfheights[i];
                if (diff >= 0)
                    std::rotate(apRecord.begin(), apRecord.begin() + diff, apRecord.end());
                else // diff < 0
                    std::rotate(apRecord.begin(), apRecord.end() + diff, apRecord.end());
            }
        }

        model_optimizer_scale(model_scaled_parameters.begin(), optimizer_parameters_begin);
        return gen_baselines;
    }

    /**
     * @brief L2-regularization wrt initial guess in optimizer scale
     *
     */
    template <typename It>
    double regularization_optimizer_scale(It parameters_begin) const
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

    VectorOfBaselines get_trimmed_baselines(const VectorOfBaselines b) const
    {
        VectorOfBaselines res;

        for (unsigned i = 0; i < b.size(); i++) {
            res.push_back({b[i].begin() + apbaselines_halfheights[i], b[i].end()});
        }

        return res;
    }

    int num_repeated_obj_runs = 1;

    template <typename It>
    double get_objective_value(It parameters_begin) const
    {
        double main_penalty = 0;
        assert(num_repeated_obj_runs > 0);
        for (int i = 0; i < num_repeated_obj_runs; i++) {
            if (!ignore_before_halfheight) {
                main_penalty += obj->dist(apbaselines, generate_baselines(parameters_begin));
            } else {
                //trim 0:halfheight_index
                main_penalty += obj->dist(get_trimmed_baselines(apbaselines), get_trimmed_baselines(generate_baselines(parameters_begin)));
            }
        }
        main_penalty /= num_repeated_obj_runs;
        double param_penalty_value = parameter_penalty_optimizer(parameters_begin);

        double res = main_penalty + param_penalty_value + regularization_optimizer_scale(parameters_begin);

        bool spontBeatCheckPassed = spontBeatCheck(parameters_begin);
        if (!spontBeatCheckPassed) {
            if (res < 0)
                throw("res should be nonnegative");
            res *= spontBeatCheckFailScaler;
        }
        if (std::isnan(res))
            res = 1e50;
        return res;
    }

    /**
     * @brief The function returns objective function value for a given unknown parameter vector
     *
     */
    /*
    template <typename It>
    double get_objective_value(It parameters_begin) const
    {
        return get_objective_value(parameters_begin, generate_baselines(parameters_begin));
    }
    */

    /**
     * @brief Method to submit result made by optimizer
     *
     */
    template <typename V>
    void submit_result(const V & parameters)
    {
        results_optimizer_format = std::vector<double>(number_unknowns);
        std::copy(parameters.begin(), parameters.begin() + number_unknowns, results_optimizer_format.begin());

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
            fill_constants_y0(model_scaled_parameters.begin(), vconstants.data(), y0.data(), baselineValues[i]);

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
    }

    template <typename It>
    void dump_ap(It parameters_begin, int i, int num_beats = -1) const
    {
        const auto tmp_b = generate_baselines(parameters_begin, num_beats);
        const double error = obj->dist(apbaselines, tmp_b);
        std::cout << "Final error: " << error << std::endl;
        std::cout << "SPONTANEOUS BEATING REQUIREMENT SATISFIED?: " << spontBeatCheck(parameters_begin) << std::endl;
        for (const auto &bs : tmp_b)
            write_baseline(bs, std::string("ap") + std::to_string(i++) + ".txt");
    }

    double loglikelihood(std::vector<double> pars) const
    {
        double ll = 0;
        VectorOfBaselines res = generate_baselines(pars.begin());
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
    std::vector<std::string> get_unique_parameter_names() const
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

    void gen_algo_stats(const std::vector<std::pair<double, int>> & sd_n_index,
            const std::vector<double> & all_genes, int gen, int total_gen)
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


#endif
