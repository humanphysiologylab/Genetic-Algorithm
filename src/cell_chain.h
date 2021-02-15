#ifndef CELL_CHAIN_H
#define CELL_CHAIN_H

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
        assert(mainCellIndex > 0);
        assert(mainCellIndex <= n);
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
        return n * cellModel.get_alg_size();
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

    void operator()(double t, double * __restrict x, double * __restrict dxdt, void * __restrict data) const
    {
        
    }

    void compute_algebraic(double t, const double *  __restrict states, double * __restrict algebraic) const
    {
        //how to calculate gap junctions?

    }

    template<typename Map1, typename Map2, typename Map3, typename Map4>
    void get_maps(Map1 & legend_states, Map2 & legend_constants, Map3 & legend_algebraic, Map4 & legend_rates) const
    {
        
    }
};

#endif
