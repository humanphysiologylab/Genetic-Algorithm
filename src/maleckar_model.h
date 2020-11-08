#ifndef MALECKAR_MODEL
#define MALECKAR_MODEL

#include <vector>
#include <cassert>

class MaleckarModel
{
    const int states_size = 30, alg_size = 70, const_size = 51;
    double * constants;
    std::vector<double> algebraic;
    void computerates(double VOI, double*  __restrict constants, double*  __restrict rates, double*  __restrict states, double* __restrict  algebraic);

public:
    void set_constants(double *c)
    {
        constants = c;
    }
   
    MaleckarModel()
    : constants(0), algebraic(alg_size)
    {}
    
    double max_step() const
    {
        return 1e-3;
    }
    int state_size() const
    {
        return states_size;
    }
    int constants_size() const
    {
        return const_size;
    }
    
    void operator()(double t, double * __restrict x, double * __restrict dxdt, void * __restrict data)
    {
        //the last parameter data was passed to lsoda_update (consider it null_ptr)
        //basically it was poor man's functor
        //here for real functor we do not need it
        assert(constants != 0);
        computerates(t, constants, dxdt, x, algebraic.data()); 
    }
    void initConsts(double * constants) const;
    void initState(double * state) const;
};
#endif
