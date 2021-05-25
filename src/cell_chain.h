#ifndef CELL_CHAIN_H
#define CELL_CHAIN_H

#include "stimulation.h"

template <typename CellModel>
class CellChainModel
{
    int n;
    int mainCellIndex;
    CellModel cellModel;
public:
    //chain of n homogeneous cells
    CellChainModel(int n, int mainCellIndex)
    : n(n), mainCellIndex(mainCellIndex)
    {
        assert(n > 0);
        assert(mainCellIndex >= 0);
        assert(mainCellIndex < n);
    }

    std::vector<std::pair<double, double>> get_r_a_tol() const
    {
		std::vector<std::pair<double, double>> r_a_tol(state_size());
		auto s = CellModel::get_r_a_tol();
		for (int i = 0; i < n; i++) {
			std::copy(s.begin(), s.end(), r_a_tol.begin() + i * cellModel.state_size());
		}
		return r_a_tol;
	}

    void set_constants(double *c)
    {
        cellModel.set_constants(c);
    }
    double max_step() const
    {
        return cellModel.max_step();
    }
    int state_size() const
    {
        return n * cellModel.state_size();
    }
    int constants_size() const
    {
        return cellModel.constants_size();
    }
    int get_alg_size() const
    {
        return cellModel.get_alg_size(); // only mainCellIndex algebraics are dumped
        //without gap junction current
    }
    void initConsts(double * constants) const
    {
        cellModel.initConsts(constants);
    }
    void initState(double * state) const
    {
        const int state_size = cellModel.state_size();
        for (int i = 0; i < n; i++) {
            cellModel.initState(state);
            state += state_size;
        }
    }


	double Igj(double g, double v_from, double v_to) const
	{
		return g * (v_from - v_to);
	}

    void operator()(double t, double * __restrict x, double * __restrict dxdt, void * __restrict data) const
    {
        // we assume V and dV/dt are at 0 indices of x and dxdt corr.
        const int v_index = 0;
        
        // cell 0 is stimulated
        const double Cm = 60; // pF // from Kernik-Clancy
		const double g_gj = 5; // gap junction conductance, nS 
        //hardcode for Kernik-Clancy
        const int stim_flag_index = 82; // stim_flag = 0: no stim
        double * constants = cellModel.get_constants();
        const double orig_stim_flag = constants[stim_flag_index];
    
        if (n == 1) {
			// single stimulated cell
			cellModel(t, x, dxdt, data);
		} else {
			for (int i = 0; i < n; i++) {
				if (i == 0) {
					//first, stimulated cell
					constants[stim_flag_index] = orig_stim_flag;
					cellModel(t, x, dxdt, data);
				} else if (i == n - 1) {
					constants[stim_flag_index] = 0;
					cellModel(t, x + i * cellModel.state_size(), dxdt + i * cellModel.state_size(), data);
				} else {
					constants[stim_flag_index] = 0;
					cellModel(t, x + i * cellModel.state_size(), dxdt + i * cellModel.state_size(), data);
				}
			}
			// now add gap junction to the right hand side
			for (int i = 0; i < n; i++) {
				if (i == 0) {
					dxdt[v_index] -= Igj(g_gj, x[v_index], x[v_index + cellModel.state_size()]) / Cm; 
				} else if (i == n - 1) {
					dxdt[v_index + i * cellModel.state_size()] -=
						Igj(g_gj, x[v_index + i * cellModel.state_size()], x[v_index + (i - 1) * cellModel.state_size()]) / Cm;
				} else {
					dxdt[v_index + i * cellModel.state_size()] -=
						(
						Igj(g_gj, x[v_index + i * cellModel.state_size()], x[v_index + (i - 1) * cellModel.state_size()])
						+
						Igj(g_gj, x[v_index + i * cellModel.state_size()], x[v_index + (i + 1) * cellModel.state_size()])
						) / Cm;
				}
			}
		}
		//restore stimulation after all
		constants[stim_flag_index] = orig_stim_flag;
    }

    void compute_algebraic(double t, const double *  __restrict states, double * __restrict algebraic) const
    {
        cellModel.compute_algebraic(t, states + cellModel.state_size() * mainCellIndex, algebraic);
    }

    template<typename Map1, typename Map2, typename Map3, typename Map4>
    void get_maps(Map1 & legend_states, Map2 & legend_constants, Map3 & legend_algebraic, Map4 & legend_rates) const
    {
		Map1 legend_states_local;
        // we are interested only in mainCellIndex cell
        cellModel.get_maps(legend_states_local, legend_constants, legend_algebraic, legend_rates);
        // guess we need to modify legend_states and legend_rates somehow
        
        // this is fine for direct problem
        // BUT it may cause issues for inverse problem!!!
        for (auto & s: legend_states_local) {
			legend_states[s.first + cellModel.state_size() * mainCellIndex] = s.second;
		}
    }
    void set_stimulation(const StimulationBase *s)
    {
		cellModel.set_stimulation(s);
	}

};

#endif
