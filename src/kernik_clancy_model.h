#ifndef KERNIK_CLANCY
#define KERNIK_CLANCY

#include <vector>
#include <cassert>
#include <cmath>
#include <string>
#include <array>

#include "stimulation.h"

class KernikClancyModel
{
    // C++ class of Kernik-Clancy iPSC-CM model
    // original description:
    // Kernik-Clancy iPSC-CM model
    //**********************************************
    //Kernik DC, Morotti S, Wu H, Garg P, Duff HJ, Kurokawa J, Jalife J, Wu JC, Grandi E, Clancy CE.
    //A computational model of induced pluripotent stem-cell derived cardiomyocytes
    //incorporating experimental variability from multiple data sources"
    //J Physiol. 2019 Jul 6. doi: 10.1113/JP277724
    //**********************************************
    //
    //Converted to C-code by Mao-Tsuen Jeng
    //
    //Colleen Clancy Lab @ UC davis
    //

    static const int states_size = 25, alg_size = 22, const_size = 99;
    static constexpr double max_step_v = 1; //(ms)
    double * constants;
    
    const StimulationBase * stimulation;
    
    void computerates(const double VOI,
                      const double *  __restrict constants,
                      double *  __restrict rates,
                      const double *  __restrict states,
                      double * __restrict algebraic = nullptr) const;

public:
    static std::vector<std::pair<double, double>> get_r_a_tol();
    void compute_algebraic(double t, const double *  __restrict states, double * __restrict algebraic) const;
    void set_constants(double *c);
	double * get_constants() const;
	void set_stimulation(const StimulationBase *s);

    KernikClancyModel();

    static double max_step();
    static int state_size();
    static int constants_size();
    static int get_alg_size();

    void operator()(double t, double * __restrict x, double * __restrict dxdt, void * __restrict data) const;
    static void initConsts(double * constants);
    static void initState(double * state);


    template<typename Map1, typename Map2, typename Map3, typename Map4>
    static void get_maps(Map1 & legend_states, Map2 & legend_constants, Map3 & legend_algebraic, Map4 & legend_rates)
    {
        //simply copy it from python's version of cellml code
        //unfortunately, there is no cellml code for Kernik-Clancy

        legend_states[0] = "V in component membrane (millivolt)";

        //Ionic Flux
        legend_states[1] = "Ca_SR (millimolar)";
        legend_states[2] = "Cai (millimolar)";
        legend_states[3] = "Nai (millimolar)";
        legend_states[4] = "Ki (millimolar)";
        legend_states[5] = "Ca_ligand (millimolar)"; //see .cpp ????????????????????????????????????????????????????????????????????

        //Current Gating (dimensionless)
        legend_states[6] = "d (activation in i_CaL)";
        legend_states[7] = "f1 (inactivation in i_CaL)";
        legend_states[8] = "fCa (calcium-dependent inactivation in i_CaL)";
        legend_states[9] = "Xr1 (activation in i_Kr)";
        legend_states[10] = "Xr2 (inactivation in i_Kr)";
        legend_states[11] = "Xs (activation in i_Ks)";
        legend_states[12] = "h (inactivation in i_Na)";
        legend_states[13] = "j (slow inactivation in i_Na)";
        legend_states[14] = "m (activation in i_Na)";
        legend_states[15] = "Xf (inactivation in i_f)";
        legend_states[16] = "s (inactivation in i_to)";
        legend_states[17] = "r (activation in i_to)";
        legend_states[18] = "dCaT (activation in i_CaT)";
        legend_states[19] = "fCaT (inactivation in i_CaT)";
        legend_states[20] = "R (in Irel)";
        legend_states[21] = "O (in Irel)";
        legend_states[22] = "I (in Irel)";
        legend_states[23] = "a_ur (ultra_rapid_K_current_aur_gate)";
        legend_states[24] = "i_ur (ultra_rapid_K_current_iur_gate)";

        legend_rates[0] = "d/dt V in component membrane (millivolt)";
        //Ionic Flux
        legend_rates[1] = "d/dt Ca_SR (millimolar)";
        legend_rates[2] = "d/dt Cai (millimolar)";
        legend_rates[3] = "d/dt Nai (millimolar)";
        legend_rates[4] = "d/dt Ki (millimolar)";
        legend_rates[5] = "d/dt Ca_ligand (millimolar)"; //see .cpp ????????????????????????????????????????????????????????????????????

        //Current Gating (dimensionless)
        legend_rates[6] = "d/dt d (activation in i_CaL)";
        legend_rates[7] = "d/dt f1 (inactivation in i_CaL)";
        legend_rates[8] = "d/dt fCa (calcium-dependent inactivation in i_CaL)";
        legend_rates[9] = "d/dt Xr1 (activation in i_Kr)";
        legend_rates[10] = "d/dt Xr2 (inactivation in i_Kr)";
        legend_rates[11] = "d/dt Xs (activation in i_Ks)";
        legend_rates[12] = "d/dt h (inactivation in i_Na)";
        legend_rates[13] = "d/dt j (slow inactivation in i_Na)";
        legend_rates[14] = "d/dt m (activation in i_Na)";
        legend_rates[15] = "d/dt Xf (inactivation in i_f)";
        legend_rates[16] = "d/dt s (inactivation in i_to)";
        legend_rates[17] = "d/dt r (activation in i_to)";
        legend_rates[18] = "d/dt dCaT (activation in i_CaT)";
        legend_rates[19] = "d/dt fCaT (inactivation in i_CaT)";
        legend_rates[20] = "d/dt R (in Irel)";
        legend_rates[21] = "d/dt O (in Irel)";
        legend_rates[22] = "d/dt I (in Irel)";
        legend_rates[23] = "d/dt a_ur (ultra_rapid_K_current_aur_gate)";
        legend_rates[24] = "d/dt i_ur (ultra_rapid_K_current_iur_gate)";

        legend_constants[0] = "g_K1_scaler (dimensionless)";
        legend_constants[1] = "g_Kr_scaler (dimensionless)";
        legend_constants[2] = "g_Ks_scaler (dimensionless)";
        legend_constants[3] = "g_to_scaler (dimensionless)";
        legend_constants[4] = "g_CaL_scaler (dimensionless)";
        legend_constants[5] = "g_CaT_scaler (dimensionless)";
        legend_constants[6] = "g_Na_scaler (dimensionless)";
        legend_constants[7] = "g_f_scaler (dimensionless)";
        legend_constants[8] = "kNaCa_scaler (dimensionless)";
        legend_constants[9] = "VmaxUp_scaler (dimensionless)";
        legend_constants[10] = "ks_scaler (dimensionless)";
        legend_constants[11] = "V_leak_scaler (dimensionless)";
        legend_constants[12] = "PNaK_scaler (dimensionless)";
        legend_constants[13] = "g_b_Na_scaler (dimensionless)";
        legend_constants[14] = "g_b_Ca_scaler (dimensionless)";
        legend_constants[15] = "g_PCa_scaler (dimensionless)";

        legend_constants[16] = "g_K1 (nS_per_pF)";
        for (int i = 17; i < 22; i++)
            legend_constants[i] = std::string("x_K1_") + std::to_string(i - 17) + " ";

        for (int i = 22; i < 33; i++)
            legend_constants[i] = std::string("x_KR_") + std::to_string(i - 22) + " ";

        for (int i = 33; i < 39; i++)
            legend_constants[i] = std::string("x_IKS_") + std::to_string(i - 33) + " ";

        for (int i = 39; i < 50; i++)
            legend_constants[i] = std::string("xTO_") + std::to_string(i - 39) + " ";

        for (int i = 50; i < 61; i++)
            legend_constants[i] = std::string("x_cal_") + std::to_string(i - 50) + " ";

        legend_constants[61] = "x_cat ";

        legend_constants[62] = "g_Na (nS_per_pF)";
        for (int i = 63; i < 76; i++)
            legend_constants[i] = std::string("x_NA_") + std::to_string(i - 63) + " ";

        for (int i = 76; i < 82; i++)
            legend_constants[i] = std::string("x_F_") + std::to_string(i - 76) + " ";

        legend_constants[82] = "stim_flag (int)"; //0 -- no stim, 1 -- rectangular pulse, 2 -- biphasic pulse
        legend_constants[83] = "stim_shift (ms)";
        legend_constants[84] = "voltageclamp ";

        legend_constants[85] = "g_kur_scaler (dimensionless)";
        legend_constants[86] = "g_kur (nS_per_pF)";
        legend_constants[87] = "stim_period (ms)";
        legend_constants[88] = "i_stim_Amplitude (pA/pF)";
        legend_constants[89] = "Na_fixed (boolean)";
        legend_constants[90] = "K_fixed (boolean)";
        legend_constants[91] = "i_stim_PulseDuration (ms)";
		legend_constants[92] = "g_b_K_scaler (dimensionless)";
        legend_constants[93] = "g_seal (nS_per_pF)";
        legend_constants[94] = "g_seal_scaler (dimensionless)";
        legend_constants[95] = "Cm_scaler (dimensionless, think really about it as a scaler for Cm/Vc)";
        legend_constants[96] = "Ma2011_IK1 (boolean)";
        legend_constants[97] = "Ma2011_INa (boolean)";
        legend_constants[98] = "abstract_seal_leak_scaler (dimensionless)";

        legend_algebraic[0] = "i_K1 (pA/pF)";
        legend_algebraic[1] = "i_to (pA/pF)";
        legend_algebraic[2] = "i_Kr (pA/pF)";
        legend_algebraic[3] = "i_Ks (pA/pF)";
        legend_algebraic[4] = "i_CaL (pA/pF)";
        legend_algebraic[5] = "i_NaK (pA/pF)";
        legend_algebraic[6] = "i_Na (pA/pF)";
        legend_algebraic[7] = "i_NaCa (pA/pF)";
        legend_algebraic[8] = "i_PCa (pA/pF)";
        legend_algebraic[9] = "i_f (pA/pF)";
        legend_algebraic[10] = "i_b_Na (pA/pF)";
        legend_algebraic[11] = "i_b_Ca (pA/pF)";
        legend_algebraic[12] = "i_rel (pA/pF)";
        legend_algebraic[13] = "i_up (pA/pF)";
        legend_algebraic[14] = "i_leak (pA/pF)";
        legend_algebraic[15] = "i_stim (pA/pF)";
        legend_algebraic[16] = "i_CaT (pA/pF)";
        legend_algebraic[17] = "i_Kur (pA/pF)";
        legend_algebraic[18] = "i_voltageclamp (pA/pF)";
        legend_algebraic[19] = "time (ms)";
        legend_algebraic[20] = "i_b_K (pA/pF)";
        legend_algebraic[21] = "i_abstract_seal_leak (pA/pF)";
    }
};

#endif
