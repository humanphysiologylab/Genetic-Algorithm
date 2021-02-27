#ifndef CELLML_ODE_SOLVER
#define CELLML_ODE_SOLVER

#include "LSODA.h"
#include <iostream>


template <typename Model, typename Table>
void make_a_record(double t, const Model & model, Table recorded_states, Table recorded_algs, const std::vector<double> & y0, int i)
{
    for (int j = 0; j < recorded_states.cols(); j++) {
        const double r = y0[recorded_states.get_model_pos(j)];
        if (isnan(r)) {
            throw -100;
        }
        recorded_states(i, j) = r;
    }

    if (recorded_algs.cols() == 0)
        return;

    //if so then compute algebraics
    std::vector<double> algebraic(model.get_alg_size());
    model.compute_algebraic(t, y0.data(), algebraic.data());

    for (int j = 0; j < recorded_algs.cols(); j++) {
        const double r = algebraic[recorded_algs.get_model_pos(j)];
        if (isnan(r)) {
            throw -100;
        }
        recorded_algs(i, j) = r;
    }
}

class ODESolver
{
public:
	template <typename Model, typename Table>
	void solve(const Model & model, std::vector<double> & y0, int & is_correct,
                double t0, double start_record, double tout, Table recorded_states, Table recorded_algs) const
	{
        //TODO for now Table is expected to be Eigen::Block, that's why it is passed by value
        //TODO throw reasonable exceptions
        //TODO record from start_record to end_record, not just to tout
		assert(recorded_states.rows() > 0);
		assert(t0 <= start_record);
		assert(start_record <= tout);
		assert(y0.size() > 0);


		LSODA lsoda;
		int istate = 1;//state to check if everything is fine

		std::vector<double> state_out;//anyway it will be resized lol

		const double record_step = (tout - start_record) / (recorded_states.rows() - 1);
		const double max_step = std::min(model.max_step(), record_step);
		double t = t0;
        int itask = 4;

		if (start_record > t)
			lsoda.lsoda_update(model, model.state_size(), y0,
				state_out, &t, start_record, &istate, nullptr, model.get_r_a_tol(), max_step, itask);
		else
			state_out = y0;

		y0 = state_out;

        make_a_record(t, model, recorded_states, recorded_algs, y0, 0);
		for (size_t i = 1; i < recorded_states.rows(); i++) {
            if (istate <= 0) {
                is_correct = 0;
                throw istate;
            }

			lsoda.lsoda_update(model, model.state_size(), y0,
				state_out, &t, t + record_step, &istate, nullptr, model.get_r_a_tol(), max_step, itask);

			y0 = state_out;
            make_a_record(t, model, recorded_states, recorded_algs, y0, i);
		}

		if (istate <= 0) {
			is_correct = 0;
			throw istate;
		} else {
			is_correct = 1;
		}
	}
};


class ODESolverCheckLSODA
{
public:
	template <typename Model, typename Table>
	void solve(const Model & model, std::vector<double> & y0, int & is_correct,
				double t0, double start_record, double tout, Table recorded_states, Table recorded_algs) const
	{
        //this is a version recording every step of lsoda solver
        //TODO throw reasonable exceptions
        //TODO record from start_record to end_record, not just to tout
		assert(recorded_states.rows() > 0);
		assert(t0 <= start_record);
		assert(start_record <= tout);
		assert(y0.size() > 0);


		LSODA lsoda;
		int istate = 1;//state to check if everything is fine

		std::vector<double> state_out;//anyway it will be resized lol

		const double max_step = 1;
		double t = t0;
        int itask = 4;

		if (start_record > t)
			lsoda.lsoda_update(model, model.state_size(), y0,
				state_out, &t, start_record, &istate, nullptr, model.get_r_a_tol(), max_step, itask);
		else
			state_out = y0;

        itask = 5; //now record every step of the solver

		y0 = state_out;

        make_a_record(t, model, recorded_states, recorded_algs, y0, 0);

        int i = 0;
        while (t < tout) {
            double tdone = t + 1;
            while (1) {
                if (istate <= 0) {
                    is_correct = 0;
                    throw istate;
                }
                if (t == tdone)
                    break;
                lsoda.lsoda_update(model, model.state_size(), y0,
                    state_out, &t, tdone, &istate, nullptr, model.get_r_a_tol(), max_step, itask);
                y0 = state_out;
                make_a_record(t, model, recorded_states, recorded_algs, y0, i);
                i++;
                if (recorded_states.rows() == i) {
                    std::cout << "not enough space to record" << std::endl;
                    throw 0;
                }
            }
        }
		if (istate <= 0) {
			is_correct = 0;
			throw istate;
		} else {
			is_correct = 1;
		}
	}
};




class MOCKODESolver
{
public:
	template <typename Model>
	void solve(const Model & model, std::vector<double> & y0, int & is_correct,
				double t0, double start_record, double tout, std::vector<double> & ap) const
	{
		assert(ap.size() > 0);
		assert(t0 <= start_record);
		assert(start_record <= tout);
		assert(y0.size() > 0);

		const double record_step = (tout - start_record) / (ap.size() - 1);
		double t = t0;

        t = start_record;
        model(y0, t);
		ap[0] = y0[0];
		for (size_t i = 1; i < ap.size(); i++) {
            model(y0, t + record_step);
            t += record_step;
			ap[i] = y0[0];
		}
        is_correct = 1;
	}
};

#endif
