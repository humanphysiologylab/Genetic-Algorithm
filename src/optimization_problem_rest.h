#ifndef OP_REST
#define OP_REST

template <typename Model, typename Solver>
class ODEoptimizationTrackVersion:
    public ODEoptimization<Model, Solver>
{
protected:
    using Base = ODEoptimization<Model, Solver>;
    using typename Base::Baseline;
    using typename Base::VectorOfBaselines;
    using Base::apbaselines;
    using Base::obj;
    using Base::write_baseline;
    using Base::generate_baselines;

    VectorOfBaselines initial_guess_baseline;
    VectorOfBaselines intermediate_baseline;

    double alpha;
public:
    using Base::is_AP_normalized;
    using Base::normalize_baseline;
    ODEoptimizationTrackVersion(const Model & model, const Solver & solver)
    : Base(model, solver)
    {}

    template <typename It>
    void start_track(It parameters_begin)
    {
        double extra_penalty;
        initial_guess_baseline = generate_baselines(parameters_begin, extra_penalty, 300);
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
    double get_objective_value(It parameters_begin) const
    {
     //   write_baseline(intermediate_baseline[0], "baseline_1.txt");///////////////////////////////
       // for (int i = 2; i < 10; i++) {///////////////////////////////////////////////////////////////////////
       double extra_penalty;
            const auto tmp_b = generate_baselines(parameters_begin, extra_penalty);
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
    void submit_result(const V & parameters)
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
