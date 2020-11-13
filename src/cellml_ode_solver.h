#ifndef CELLML_ODE_SOLVER
#define CELLML_ODE_SOLVER

#include "LSODA.h"
#include <iostream>



class ODESolver
{
	const double rtol = 1e-9, atol = 1e-6;

public:
	template <typename Model>
	void solve(const Model & model, std::vector<double> & y0, int & is_correct,
				double t0, double start_record, double tout, std::vector<double> & ap) const
	{
		assert(ap.size() > 0);
		assert(t0 <= start_record);
		assert(start_record <= tout);
		assert(y0.size() > 0);

		LSODA lsoda;
		int istate = 1;//state to check if everything is fine

		std::vector<double> state_out;//anyway it will be resized lol

		const double record_step = (tout - start_record) / (ap.size() - 1);
		const double max_step = std::min(model.max_step(), record_step);
		double t = t0;
		int itask = 4;

		if (start_record > t)
			lsoda.lsoda_update(model, model.state_size(), y0,
				state_out, &t, start_record, &istate, nullptr, rtol, atol, max_step, itask);
		else
			state_out = y0;

		y0 = state_out;
		ap[0] = state_out[0];
		for (size_t i = 1; i < ap.size(); i++) {

			lsoda.lsoda_update(model, model.state_size(), y0,
				state_out, &t, t + record_step, &istate, nullptr, rtol, atol, max_step, itask);

			y0 = state_out;
			ap[i] = state_out[0];
		}

		if (istate <= 0) {
			is_correct = 0;
			throw;
		} else {
			is_correct = 1;
		}
	}
};
#endif
