#ifndef MALECKAR_MODEL
#define MALECKAR_MODEL

#include <vector>
#include <cassert>

class MaleckarModel
{
    static const int states_size = 30, alg_size = 70, const_size = 51;
    double * constants;
    void computerates(double VOI, const double*  __restrict constants, double*  __restrict rates, double*  __restrict states) const;

public:
    void set_constants(double *c)
    {
        constants = c;
    }
   
    MaleckarModel()
    : constants(0)
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
    int get_alg_size() const
    {
        return alg_size;
    }

    void operator()(double t, double * __restrict x, double * __restrict dxdt, void * __restrict data) const
    {
        //the last parameter data was passed to lsoda_update (consider it null_ptr)
        //basically it was poor man's functor
        //here for real functor we do not need it
        assert(constants != 0);
        computerates(t, constants, dxdt, x); 
    }
    void initConsts(double * constants) const;
    void initState(double * state) const;


    template<typename Map1, typename Map2, typename Map3, typename Map4>
    void get_maps(Map1 & legend_states, Map2 & legend_constants, Map3 & legend_algebraic, Map4 & legend_rates) const
    {
        //simply copy it from python's version of cellml code
        legend_states[0] = "V in component membrane (millivolt)";
        legend_constants[0] = "R in component membrane (millijoule_per_mole_kelvin)";
        legend_constants[1] = "T in component membrane (kelvin)";
        legend_constants[2] = "F in component membrane (coulomb_per_mole)";
        legend_constants[3] = "Cm in component membrane (nanoF)";
        legend_algebraic[0] = "Q_tot in component membrane (millivolt)";
        legend_algebraic[37] = "i_Na in component sodium_current (picoA)";
        legend_algebraic[41] = "i_Ca_L in component L_type_Ca_channel (picoA)";
        legend_algebraic[44] = "i_t in component Ca_independent_transient_outward_K_current (picoA)";
        legend_algebraic[45] = "i_Kur in component ultra_rapid_K_current (picoA)";
        legend_algebraic[46] = "i_K1 in component inward_rectifier (picoA)";
        legend_algebraic[49] = "i_Kr in component delayed_rectifier_K_currents (picoA)";
        legend_algebraic[47] = "i_Ks in component delayed_rectifier_K_currents (picoA)";
        legend_algebraic[50] = "i_B_Na in component background_currents (picoA)";
        legend_algebraic[52] = "i_B_Ca in component background_currents (picoA)";
        legend_algebraic[54] = "i_NaK in component sodium_potassium_pump (picoA)";
        legend_algebraic[55] = "i_CaP in component sarcolemmal_calcium_pump_current (picoA)";
        legend_algebraic[56] = "i_NaCa in component Na_Ca_ion_exchanger_current (picoA)";
        legend_algebraic[57] = "i_KACh in component ACh_dependent_K_current (picoA)";
        legend_algebraic[59] = "I in component membrane (pA_per_nF)";
        legend_algebraic[24] = "i_Stim in component membrane (pA_per_nF)";
        legend_constants[4] = "stim_offset in component membrane (second)";
        legend_constants[5] = "stim_period in component membrane (second)";
        legend_constants[6] = "stim_duration in component membrane (second)";
        legend_constants[7] = "stim_amplitude in component membrane (pA_per_nF)";
        legend_algebraic[1] = "past in component membrane (second)";
        legend_algebraic[35] = "E_Na in component sodium_current (millivolt)";
        legend_constants[8] = "P_Na in component sodium_current (nanolitre_per_second)";
        legend_states[1] = "Na_c in component cleft_space_ion_concentrations (millimolar)";
        legend_states[2] = "Na_i in component intracellular_ion_concentrations (millimolar)";
        legend_states[3] = "m in component sodium_current_m_gate (dimensionless)";
        legend_states[4] = "h1 in component sodium_current_h1_gate (dimensionless)";
        legend_states[5] = "h2 in component sodium_current_h2_gate (dimensionless)";
        legend_algebraic[14] = "m_infinity in component sodium_current_m_gate (dimensionless)";
        legend_algebraic[2] = "m_factor in component sodium_current_m_gate (dimensionless)";
        legend_algebraic[26] = "tau_m in component sodium_current_m_gate (second)";
        legend_algebraic[3] = "h_infinity in component sodium_current_h1_gate (dimensionless)";
        legend_algebraic[15] = "h_factor in component sodium_current_h1_gate (dimensionless)";
        legend_algebraic[27] = "tau_h1 in component sodium_current_h1_gate (second)";
        legend_algebraic[28] = "tau_h2 in component sodium_current_h2_gate (second)";
        legend_constants[9] = "g_Ca_L in component L_type_Ca_channel (nanoS)";
        legend_constants[10] = "E_Ca_app(not_fixed!) in component L_type_Ca_channel (millivolt)";
        legend_algebraic[39] = "f_Ca in component L_type_Ca_channel (dimensionless)";
        legend_constants[11] = "k_Ca in component L_type_Ca_channel (millimolar)";
        legend_states[6] = "Ca_d in component intracellular_ion_concentrations (millimolar)";
        legend_states[7] = "d_L in component L_type_Ca_channel_d_L_gate (dimensionless)";
        legend_states[8] = "f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless)";
        legend_states[9] = "f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless)";
        legend_algebraic[4] = "d_L_infinity in component L_type_Ca_channel_d_L_gate (dimensionless)";
        legend_algebraic[16] = "d_L_factor in component L_type_Ca_channel_d_L_gate (dimensionless)";
        legend_algebraic[29] = "tau_d_L in component L_type_Ca_channel_d_L_gate (second)";
        legend_algebraic[5] = "f_L_infinity in component L_type_Ca_channel_f_L1_gate (dimensionless)";
        legend_algebraic[17] = "f_L_factor in component L_type_Ca_channel_f_L1_gate (millivolt)";
        legend_algebraic[30] = "tau_f_L1 in component L_type_Ca_channel_f_L1_gate (second)";
        legend_algebraic[31] = "tau_f_L2 in component L_type_Ca_channel_f_L2_gate (second)";
        legend_algebraic[43] = "E_K in component Ca_independent_transient_outward_K_current (millivolt)";
        legend_constants[12] = "g_t in component Ca_independent_transient_outward_K_current (nanoS)";
        legend_states[10] = "K_c in component cleft_space_ion_concentrations (millimolar)";
        legend_states[11] = "K_i in component intracellular_ion_concentrations (millimolar)";
        legend_states[12] = "r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless)";
        legend_states[13] = "s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)";
        legend_algebraic[18] = "tau_r in component Ca_independent_transient_outward_K_current_r_gate (second)";
        legend_algebraic[6] = "r_infinity in component Ca_independent_transient_outward_K_current_r_gate (dimensionless)";
        legend_algebraic[32] = "tau_s in component Ca_independent_transient_outward_K_current_s_gate (second)";
        legend_algebraic[7] = "s_infinity in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)";
        legend_algebraic[19] = "s_factor in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)";
        legend_constants[13] = "g_kur in component ultra_rapid_K_current (nanoS)";
        legend_states[14] = "a_ur in component ultra_rapid_K_current_aur_gate (dimensionless)";
        legend_states[15] = "i_ur in component ultra_rapid_K_current_iur_gate (dimensionless)";
        legend_algebraic[8] = "a_ur_infinity in component ultra_rapid_K_current_aur_gate (dimensionless)";
        legend_algebraic[20] = "tau_a_ur in component ultra_rapid_K_current_aur_gate (second)";
        legend_algebraic[9] = "i_ur_infinity in component ultra_rapid_K_current_iur_gate (dimensionless)";
        legend_algebraic[21] = "tau_i_ur in component ultra_rapid_K_current_iur_gate (second)";
        legend_constants[14] = "g_K1 in component inward_rectifier (nanoS)";
        legend_constants[15] = "g_Ks in component delayed_rectifier_K_currents (nanoS)";
        legend_constants[16] = "g_Kr in component delayed_rectifier_K_currents (nanoS)";
        legend_states[16] = "n in component delayed_rectifier_K_currents_n_gate (dimensionless)";
        legend_states[17] = "pa in component delayed_rectifier_K_currents_pa_gate (dimensionless)";
        legend_algebraic[48] = "pip in component delayed_rectifier_K_currents_pi_gate (dimensionless)";
        legend_algebraic[33] = "tau_n in component delayed_rectifier_K_currents_n_gate (second)";
        legend_algebraic[10] = "n_infinity in component delayed_rectifier_K_currents_n_gate (dimensionless)";
        legend_algebraic[22] = "n_factor in component delayed_rectifier_K_currents_n_gate (dimensionless)";
        legend_algebraic[34] = "tau_pa in component delayed_rectifier_K_currents_pa_gate (second)";
        legend_algebraic[23] = "pa_factor in component delayed_rectifier_K_currents_pa_gate (dimensionless)";
        legend_algebraic[11] = "p_a_infinity in component delayed_rectifier_K_currents_pa_gate (dimensionless)";
        legend_constants[17] = "g_B_Na in component background_currents (nanoS)";
        legend_constants[18] = "g_B_Ca in component background_currents (nanoS)";
        legend_algebraic[51] = "E_Ca in component background_currents (millivolt)";
        legend_states[18] = "Ca_c in component cleft_space_ion_concentrations (millimolar)";
        legend_states[19] = "Ca_i in component intracellular_ion_concentrations (millimolar)";
        legend_constants[19] = "K_NaK_K in component sodium_potassium_pump (millimolar)";
        legend_constants[20] = "i_NaK_max in component sodium_potassium_pump (picoA)";
        legend_constants[21] = "pow_K_NaK_Na_15 in component sodium_potassium_pump (millimolar15)";
        legend_algebraic[53] = "pow_Na_i_15 in component sodium_potassium_pump (millimolar15)";
        legend_constants[22] = "i_CaP_max in component sarcolemmal_calcium_pump_current (picoA)";
        legend_constants[23] = "k_CaP in component sarcolemmal_calcium_pump_current (millimolar)";
        legend_constants[24] = "K_NaCa in component Na_Ca_ion_exchanger_current (picoA_per_millimolar_4)";
        legend_constants[25] = "d_NaCa in component Na_Ca_ion_exchanger_current (per_millimolar_4)";
        legend_constants[26] = "gamma_Na in component Na_Ca_ion_exchanger_current (dimensionless)";
        legend_constants[27] = "ACh in component ACh_dependent_K_current (millimolar)";
        legend_constants[28] = "phi_Na_en in component intracellular_ion_concentrations (picoA)";
        legend_constants[29] = "Vol_i in component intracellular_ion_concentrations (nanolitre)";
        legend_constants[30] = "Vol_d in component intracellular_ion_concentrations (nanolitre)";
        legend_algebraic[58] = "i_di in component intracellular_ion_concentrations (picoA)";
        legend_constants[31] = "tau_di in component intracellular_ion_concentrations (second)";
        legend_algebraic[67] = "i_up in component Ca_handling_by_the_SR (picoA)";
        legend_algebraic[66] = "i_rel in component Ca_handling_by_the_SR (picoA)";
        legend_algebraic[63] = "J_O in component intracellular_Ca_buffering (per_second)";
        legend_states[20] = "O_C in component intracellular_Ca_buffering (dimensionless)";
        legend_states[21] = "O_TC in component intracellular_Ca_buffering (dimensionless)";
        legend_states[22] = "O_TMgC in component intracellular_Ca_buffering (dimensionless)";
        legend_states[23] = "O_TMgMg in component intracellular_Ca_buffering (dimensionless)";
        legend_states[24] = "O in component intracellular_Ca_buffering (dimensionless)";
        legend_algebraic[60] = "J_O_C in component intracellular_Ca_buffering (per_second)";
        legend_algebraic[61] = "J_O_TC in component intracellular_Ca_buffering (per_second)";
        legend_algebraic[62] = "J_O_TMgC in component intracellular_Ca_buffering (per_second)";
        legend_algebraic[12] = "J_O_TMgMg in component intracellular_Ca_buffering (per_second)";
        legend_constants[32] = "Mg_i in component intracellular_Ca_buffering (millimolar)";
        legend_constants[33] = "Vol_c in component cleft_space_ion_concentrations (nanolitre)";
        legend_constants[34] = "tau_Na in component cleft_space_ion_concentrations (second)";
        legend_constants[35] = "tau_K in component cleft_space_ion_concentrations (second)";
        legend_constants[36] = "tau_Ca in component cleft_space_ion_concentrations (second)";
        legend_constants[37] = "Na_b in component cleft_space_ion_concentrations (millimolar)";
        legend_constants[38] = "Ca_b in component cleft_space_ion_concentrations (millimolar)";
        legend_constants[39] = "K_b in component cleft_space_ion_concentrations (millimolar)";
        legend_algebraic[68] = "i_tr in component Ca_handling_by_the_SR (picoA)";
        legend_constants[40] = "I_up_max in component Ca_handling_by_the_SR (picoA)";
        legend_constants[41] = "k_cyca in component Ca_handling_by_the_SR (millimolar)";
        legend_constants[42] = "k_srca in component Ca_handling_by_the_SR (millimolar)";
        legend_constants[43] = "k_xcs in component Ca_handling_by_the_SR (dimensionless)";
        legend_constants[44] = "alpha_rel in component Ca_handling_by_the_SR (picoA_per_millimolar)";
        legend_states[25] = "Ca_rel in component Ca_handling_by_the_SR (millimolar)";
        legend_states[26] = "Ca_up in component Ca_handling_by_the_SR (millimolar)";
        legend_constants[45] = "Vol_up in component Ca_handling_by_the_SR (nanolitre)";
        legend_constants[46] = "Vol_rel in component Ca_handling_by_the_SR (nanolitre)";
        legend_algebraic[40] = "r_act in component Ca_handling_by_the_SR (per_second)";
        legend_algebraic[42] = "r_inact in component Ca_handling_by_the_SR (per_second)";
        legend_constants[47] = "r_recov in component Ca_handling_by_the_SR (per_second)";
        legend_algebraic[13] = "r_Ca_d_term in component Ca_handling_by_the_SR (dimensionless)";
        legend_algebraic[25] = "r_Ca_i_term in component Ca_handling_by_the_SR (dimensionless)";
        legend_algebraic[36] = "r_Ca_d_factor in component Ca_handling_by_the_SR (dimensionless)";
        legend_algebraic[38] = "r_Ca_i_factor in component Ca_handling_by_the_SR (dimensionless)";
        legend_algebraic[64] = "i_rel_f2 in component Ca_handling_by_the_SR (dimensionless)";
        legend_algebraic[65] = "i_rel_factor in component Ca_handling_by_the_SR (dimensionless)";
        legend_states[27] = "O_Calse in component Ca_handling_by_the_SR (dimensionless)";
        legend_algebraic[69] = "J_O_Calse in component Ca_handling_by_the_SR (per_second)";
        legend_states[28] = "F1 in component Ca_handling_by_the_SR (dimensionless)";
        legend_states[29] = "F2 in component Ca_handling_by_the_SR (dimensionless)";
        legend_constants[48] = "tau_tr in component Ca_handling_by_the_SR (second)";
        legend_constants[49] = "k_rel_i in component Ca_handling_by_the_SR (millimolar)";
        legend_constants[50] = "k_rel_d in component Ca_handling_by_the_SR (millimolar)";
        legend_rates[0] = "d/dt V in component membrane (millivolt)";
        legend_rates[3] = "d/dt m in component sodium_current_m_gate (dimensionless)";
        legend_rates[4] = "d/dt h1 in component sodium_current_h1_gate (dimensionless)";
        legend_rates[5] = "d/dt h2 in component sodium_current_h2_gate (dimensionless)";
        legend_rates[7] = "d/dt d_L in component L_type_Ca_channel_d_L_gate (dimensionless)";
        legend_rates[8] = "d/dt f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless)";
        legend_rates[9] = "d/dt f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless)";
        legend_rates[12] = "d/dt r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless)";
        legend_rates[13] = "d/dt s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)";
        legend_rates[14] = "d/dt a_ur in component ultra_rapid_K_current_aur_gate (dimensionless)";
        legend_rates[15] = "d/dt i_ur in component ultra_rapid_K_current_iur_gate (dimensionless)";
        legend_rates[16] = "d/dt n in component delayed_rectifier_K_currents_n_gate (dimensionless)";
        legend_rates[17] = "d/dt pa in component delayed_rectifier_K_currents_pa_gate (dimensionless)";
        legend_rates[11] = "d/dt K_i in component intracellular_ion_concentrations (millimolar)";
        legend_rates[2] = "d/dt Na_i in component intracellular_ion_concentrations (millimolar)";
        legend_rates[19] = "d/dt Ca_i in component intracellular_ion_concentrations (millimolar)";
        legend_rates[6] = "d/dt Ca_d in component intracellular_ion_concentrations (millimolar)";
        legend_rates[20] = "d/dt O_C in component intracellular_Ca_buffering (dimensionless)";
        legend_rates[21] = "d/dt O_TC in component intracellular_Ca_buffering (dimensionless)";
        legend_rates[22] = "d/dt O_TMgC in component intracellular_Ca_buffering (dimensionless)";
        legend_rates[23] = "d/dt O_TMgMg in component intracellular_Ca_buffering (dimensionless)";
        legend_rates[24] = "d/dt O in component intracellular_Ca_buffering (dimensionless)";
        legend_rates[18] = "d/dt Ca_c in component cleft_space_ion_concentrations (millimolar)";
        legend_rates[10] = "d/dt K_c in component cleft_space_ion_concentrations (millimolar)";
        legend_rates[1] = "d/dt Na_c in component cleft_space_ion_concentrations (millimolar)";
        legend_rates[28] = "d/dt F1 in component Ca_handling_by_the_SR (dimensionless)";
        legend_rates[29] = "d/dt F2 in component Ca_handling_by_the_SR (dimensionless)";
        legend_rates[27] = "d/dt O_Calse in component Ca_handling_by_the_SR (dimensionless)";
        legend_rates[26] = "d/dt Ca_up in component Ca_handling_by_the_SR (millimolar)";
        legend_rates[25] = "d/dt Ca_rel in component Ca_handling_by_the_SR (millimolar)";
    }
};
#endif
