#include "maleckar_model.h"
#include <cmath>
/*
   There are a total of 70 entries in the algebraic variable array.
   There are a total of 30 entries in each of the rate and state variable arrays.
   There are a total of 51 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (second).
 * states[0] is V in component membrane (millivolt).
 * constants[0] is R in component membrane (millijoule_per_mole_kelvin).
 * constants[1] is T in component membrane (kelvin).
 * constants[2] is F in component membrane (coulomb_per_mole).
 * constants[3] is Cm in component membrane (nanoF).
 * algebraic[0] is Q_tot in component membrane (millivolt).
 * algebraic[37] is i_Na in component sodium_current (picoA).
 * algebraic[41] is i_Ca_L in component L_type_Ca_channel (picoA).
 * algebraic[44] is i_t in component Ca_independent_transient_outward_K_current (picoA).
 * algebraic[45] is i_Kur in component ultra_rapid_K_current (picoA).
 * algebraic[46] is i_K1 in component inward_rectifier (picoA).
 * algebraic[49] is i_Kr in component delayed_rectifier_K_currents (picoA).
 * algebraic[47] is i_Ks in component delayed_rectifier_K_currents (picoA).
 * algebraic[50] is i_B_Na in component background_currents (picoA).
 * algebraic[52] is i_B_Ca in component background_currents (picoA).
 * algebraic[54] is i_NaK in component sodium_potassium_pump (picoA).
 * algebraic[55] is i_CaP in component sarcolemmal_calcium_pump_current (picoA).
 * algebraic[56] is i_NaCa in component Na_Ca_ion_exchanger_current (picoA).
 * algebraic[57] is i_KACh in component ACh_dependent_K_current (picoA).
 * algebraic[59] is I in component membrane (pA_per_nF).
 * algebraic[24] is i_Stim in component membrane (pA_per_nF).
 * constants[4] is stim_offset in component membrane (second).
 * constants[5] is stim_period in component membrane (second).
 * constants[6] is stim_duration in component membrane (second).
 * constants[7] is stim_amplitude in component membrane (pA_per_nF).
 * algebraic[1] is past in component membrane (second).
 * algebraic[35] is E_Na in component sodium_current (millivolt).
 * constants[8] is P_Na in component sodium_current (nanolitre_per_second).
 * states[1] is Na_c in component cleft_space_ion_concentrations (millimolar).
 * states[2] is Na_i in component intracellular_ion_concentrations (millimolar).
 * states[3] is m in component sodium_current_m_gate (dimensionless).
 * states[4] is h1 in component sodium_current_h1_gate (dimensionless).
 * states[5] is h2 in component sodium_current_h2_gate (dimensionless).
 * algebraic[14] is m_infinity in component sodium_current_m_gate (dimensionless).
 * algebraic[2] is m_factor in component sodium_current_m_gate (dimensionless).
 * algebraic[26] is tau_m in component sodium_current_m_gate (second).
 * algebraic[3] is h_infinity in component sodium_current_h1_gate (dimensionless).
 * algebraic[15] is h_factor in component sodium_current_h1_gate (dimensionless).
 * algebraic[27] is tau_h1 in component sodium_current_h1_gate (second).
 * algebraic[28] is tau_h2 in component sodium_current_h2_gate (second).
 * constants[9] is g_Ca_L in component L_type_Ca_channel (nanoS).
 * constants[10] is E_Ca_app in component L_type_Ca_channel (millivolt).     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 * algebraic[39] is f_Ca in component L_type_Ca_channel (dimensionless).
 * constants[11] is k_Ca in component L_type_Ca_channel (millimolar).
 * states[6] is Ca_d in component intracellular_ion_concentrations (millimolar).
 * states[7] is d_L in component L_type_Ca_channel_d_L_gate (dimensionless).
 * states[8] is f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless).
 * states[9] is f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless).
 * algebraic[4] is d_L_infinity in component L_type_Ca_channel_d_L_gate (dimensionless).
 * algebraic[16] is d_L_factor in component L_type_Ca_channel_d_L_gate (dimensionless).
 * algebraic[29] is tau_d_L in component L_type_Ca_channel_d_L_gate (second).
 * algebraic[5] is f_L_infinity in component L_type_Ca_channel_f_L1_gate (dimensionless).
 * algebraic[17] is f_L_factor in component L_type_Ca_channel_f_L1_gate (millivolt).
 * algebraic[30] is tau_f_L1 in component L_type_Ca_channel_f_L1_gate (second).
 * algebraic[31] is tau_f_L2 in component L_type_Ca_channel_f_L2_gate (second).
 * algebraic[43] is E_K in component Ca_independent_transient_outward_K_current (millivolt).
 * constants[12] is g_t in component Ca_independent_transient_outward_K_current (nanoS).
 * states[10] is K_c in component cleft_space_ion_concentrations (millimolar).
 * states[11] is K_i in component intracellular_ion_concentrations (millimolar).
 * states[12] is r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
 * states[13] is s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
 * algebraic[18] is tau_r in component Ca_independent_transient_outward_K_current_r_gate (second).
 * algebraic[6] is r_infinity in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
 * algebraic[32] is tau_s in component Ca_independent_transient_outward_K_current_s_gate (second).
 * algebraic[7] is s_infinity in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
 * algebraic[19] is s_factor in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
 * constants[13] is g_kur in component ultra_rapid_K_current (nanoS).
 * states[14] is a_ur in component ultra_rapid_K_current_aur_gate (dimensionless).
 * states[15] is i_ur in component ultra_rapid_K_current_iur_gate (dimensionless).
 * algebraic[8] is a_ur_infinity in component ultra_rapid_K_current_aur_gate (dimensionless).
 * algebraic[20] is tau_a_ur in component ultra_rapid_K_current_aur_gate (second).
 * algebraic[9] is i_ur_infinity in component ultra_rapid_K_current_iur_gate (dimensionless).
 * algebraic[21] is tau_i_ur in component ultra_rapid_K_current_iur_gate (second).
 * constants[14] is g_K1 in component inward_rectifier (nanoS).
 * constants[15] is g_Ks in component delayed_rectifier_K_currents (nanoS).
 * constants[16] is g_Kr in component delayed_rectifier_K_currents (nanoS).
 * states[16] is n in component delayed_rectifier_K_currents_n_gate (dimensionless).
 * states[17] is pa in component delayed_rectifier_K_currents_pa_gate (dimensionless).
 * algebraic[48] is pip in component delayed_rectifier_K_currents_pi_gate (dimensionless).
 * algebraic[33] is tau_n in component delayed_rectifier_K_currents_n_gate (second).
 * algebraic[10] is n_infinity in component delayed_rectifier_K_currents_n_gate (dimensionless).
 * algebraic[22] is n_factor in component delayed_rectifier_K_currents_n_gate (dimensionless).
 * algebraic[34] is tau_pa in component delayed_rectifier_K_currents_pa_gate (second).
 * algebraic[23] is pa_factor in component delayed_rectifier_K_currents_pa_gate (dimensionless).
 * algebraic[11] is p_a_infinity in component delayed_rectifier_K_currents_pa_gate (dimensionless).
 * constants[17] is g_B_Na in component background_currents (nanoS).
 * constants[18] is g_B_Ca in component background_currents (nanoS).
 * algebraic[51] is E_Ca in component background_currents (millivolt).
 * states[18] is Ca_c in component cleft_space_ion_concentrations (millimolar).
 * states[19] is Ca_i in component intracellular_ion_concentrations (millimolar).
 * constants[19] is K_NaK_K in component sodium_potassium_pump (millimolar).
 * constants[20] is i_NaK_max in component sodium_potassium_pump (picoA).
 * constants[21] is pow_K_NaK_Na_15 in component sodium_potassium_pump (millimolar15).
 * algebraic[53] is pow_Na_i_15 in component sodium_potassium_pump (millimolar15).
 * constants[22] is i_CaP_max in component sarcolemmal_calcium_pump_current (picoA).
 * constants[23] is k_CaP in component sarcolemmal_calcium_pump_current (millimolar).
 * constants[24] is K_NaCa in component Na_Ca_ion_exchanger_current (picoA_per_millimolar_4).
 * constants[25] is d_NaCa in component Na_Ca_ion_exchanger_current (per_millimolar_4).
 * constants[26] is gamma_Na in component Na_Ca_ion_exchanger_current (dimensionless).
 * constants[27] is ACh in component ACh_dependent_K_current (millimolar).
 * constants[28] is phi_Na_en in component intracellular_ion_concentrations (picoA).
 * constants[29] is Vol_i in component intracellular_ion_concentrations (nanolitre).
 * constants[30] is Vol_d in component intracellular_ion_concentrations (nanolitre).
 * algebraic[58] is i_di in component intracellular_ion_concentrations (picoA).
 * constants[31] is tau_di in component intracellular_ion_concentrations (second).
 * algebraic[67] is i_up in component Ca_handling_by_the_SR (picoA).
 * algebraic[66] is i_rel in component Ca_handling_by_the_SR (picoA).
 * algebraic[63] is J_O in component intracellular_Ca_buffering (per_second).
 * states[20] is O_C in component intracellular_Ca_buffering (dimensionless).
 * states[21] is O_TC in component intracellular_Ca_buffering (dimensionless).
 * states[22] is O_TMgC in component intracellular_Ca_buffering (dimensionless).
 * states[23] is O_TMgMg in component intracellular_Ca_buffering (dimensionless).
 * states[24] is O in component intracellular_Ca_buffering (dimensionless).
 * algebraic[60] is J_O_C in component intracellular_Ca_buffering (per_second).
 * algebraic[61] is J_O_TC in component intracellular_Ca_buffering (per_second).
 * algebraic[62] is J_O_TMgC in component intracellular_Ca_buffering (per_second).
 * algebraic[12] is J_O_TMgMg in component intracellular_Ca_buffering (per_second).
 * constants[32] is Mg_i in component intracellular_Ca_buffering (millimolar).
 * constants[33] is Vol_c in component cleft_space_ion_concentrations (nanolitre).
 * constants[34] is tau_Na in component cleft_space_ion_concentrations (second).
 * constants[35] is tau_K in component cleft_space_ion_concentrations (second).
 * constants[36] is tau_Ca in component cleft_space_ion_concentrations (second).
 * constants[37] is Na_b in component cleft_space_ion_concentrations (millimolar).
 * constants[38] is Ca_b in component cleft_space_ion_concentrations (millimolar).
 * constants[39] is K_b in component cleft_space_ion_concentrations (millimolar).
 * algebraic[68] is i_tr in component Ca_handling_by_the_SR (picoA).
 * constants[40] is I_up_max in component Ca_handling_by_the_SR (picoA).
 * constants[41] is k_cyca in component Ca_handling_by_the_SR (millimolar).
 * constants[42] is k_srca in component Ca_handling_by_the_SR (millimolar).
 * constants[43] is k_xcs in component Ca_handling_by_the_SR (dimensionless).
 * constants[44] is alpha_rel in component Ca_handling_by_the_SR (picoA_per_millimolar).
 * states[25] is Ca_rel in component Ca_handling_by_the_SR (millimolar).
 * states[26] is Ca_up in component Ca_handling_by_the_SR (millimolar).
 * constants[45] is Vol_up in component Ca_handling_by_the_SR (nanolitre).
 * constants[46] is Vol_rel in component Ca_handling_by_the_SR (nanolitre).
 * algebraic[40] is r_act in component Ca_handling_by_the_SR (per_second).
 * algebraic[42] is r_inact in component Ca_handling_by_the_SR (per_second).
 * constants[47] is r_recov in component Ca_handling_by_the_SR (per_second).
 * algebraic[13] is r_Ca_d_term in component Ca_handling_by_the_SR (dimensionless).
 * algebraic[25] is r_Ca_i_term in component Ca_handling_by_the_SR (dimensionless).
 * algebraic[36] is r_Ca_d_factor in component Ca_handling_by_the_SR (dimensionless).
 * algebraic[38] is r_Ca_i_factor in component Ca_handling_by_the_SR (dimensionless).
 * algebraic[64] is i_rel_f2 in component Ca_handling_by_the_SR (dimensionless).
 * algebraic[65] is i_rel_factor in component Ca_handling_by_the_SR (dimensionless).
 * states[27] is O_Calse in component Ca_handling_by_the_SR (dimensionless).
 * algebraic[69] is J_O_Calse in component Ca_handling_by_the_SR (per_second).
 * states[28] is F1 in component Ca_handling_by_the_SR (dimensionless).
 * states[29] is F2 in component Ca_handling_by_the_SR (dimensionless).
 * constants[48] is tau_tr in component Ca_handling_by_the_SR (second).
 * constants[49] is k_rel_i in component Ca_handling_by_the_SR (millimolar).
 * constants[50] is k_rel_d in component Ca_handling_by_the_SR (millimolar).
 * rates[0] is d/dt V in component membrane (millivolt).
 * rates[3] is d/dt m in component sodium_current_m_gate (dimensionless).
 * rates[4] is d/dt h1 in component sodium_current_h1_gate (dimensionless).
 * rates[5] is d/dt h2 in component sodium_current_h2_gate (dimensionless).
 * rates[7] is d/dt d_L in component L_type_Ca_channel_d_L_gate (dimensionless).
 * rates[8] is d/dt f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless).
 * rates[9] is d/dt f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless).
 * rates[12] is d/dt r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
 * rates[13] is d/dt s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
 * rates[14] is d/dt a_ur in component ultra_rapid_K_current_aur_gate (dimensionless).
 * rates[15] is d/dt i_ur in component ultra_rapid_K_current_iur_gate (dimensionless).
 * rates[16] is d/dt n in component delayed_rectifier_K_currents_n_gate (dimensionless).
 * rates[17] is d/dt pa in component delayed_rectifier_K_currents_pa_gate (dimensionless).
 * rates[11] is d/dt K_i in component intracellular_ion_concentrations (millimolar).
 * rates[2] is d/dt Na_i in component intracellular_ion_concentrations (millimolar).
 * rates[19] is d/dt Ca_i in component intracellular_ion_concentrations (millimolar).
 * rates[6] is d/dt Ca_d in component intracellular_ion_concentrations (millimolar).
 * rates[20] is d/dt O_C in component intracellular_Ca_buffering (dimensionless).
 * rates[21] is d/dt O_TC in component intracellular_Ca_buffering (dimensionless).
 * rates[22] is d/dt O_TMgC in component intracellular_Ca_buffering (dimensionless).
 * rates[23] is d/dt O_TMgMg in component intracellular_Ca_buffering (dimensionless).
 * rates[24] is d/dt O in component intracellular_Ca_buffering (dimensionless).
 * rates[18] is d/dt Ca_c in component cleft_space_ion_concentrations (millimolar).
 * rates[10] is d/dt K_c in component cleft_space_ion_concentrations (millimolar).
 * rates[1] is d/dt Na_c in component cleft_space_ion_concentrations (millimolar).
 * rates[28] is d/dt F1 in component Ca_handling_by_the_SR (dimensionless).
 * rates[29] is d/dt F2 in component Ca_handling_by_the_SR (dimensionless).
 * rates[27] is d/dt O_Calse in component Ca_handling_by_the_SR (dimensionless).
 * rates[26] is d/dt Ca_up in component Ca_handling_by_the_SR (millimolar).
 * rates[25] is d/dt Ca_rel in component Ca_handling_by_the_SR (millimolar).
 */

void MaleckarModel::initState(double * states) const
{
    states[0] = -74.031982;
    states[1] = 130.022096;
    states[2] = 8.516766;
    states[3] = 0.003289;
    states[4] = 0.877202;
    states[5] = 0.873881;
    states[6] = 7.1e-5;
    states[7] = 0.000014;
    states[8] = 0.998597;
    states[9] = 0.998586;
    states[10] = 5.560224;
    states[11] = 129.485991;
    states[12] = 0.001089;
    states[13] = 0.948597;
    states[14] = 0.000367;
    states[15] = 0.96729;
    states[16] = 0.004374;
    states[17] = 0.000053;
    states[18] = 1.815768;
    states[19] = 6.5e-5;
    states[20] = 0.026766;
    states[21] = 0.012922;
    states[22] = 0.190369;
    states[23] = 0.714463;
    states[24] = 1.38222;
    states[25] = 0.632613;
    states[26] = 0.649195;
    states[27] = 0.431547;
    states[28] = 0.470055;
    states[29] = 0.002814;
}

void MaleckarModel::initConsts(double * constants) const
{
    constants[0] = 8314;
    constants[1] = 306.15;
    constants[2] = 96487;
    constants[3] = 50;
    constants[4] = 0;
    constants[5] = 1;
    constants[6] = 0.006;
    constants[7] = -15;
    constants[8] = 0.0018;
    constants[9] = 6.75;
    constants[10] = 60;
    constants[11] = 0.025;
    constants[12] = 8.25;
    constants[13] = 2.25;
    constants[14] = 3.1;
    constants[15] = 1;
    constants[16] = 0.5;
    constants[17] = 0.060599;
    constants[18] = 0.078681;
    constants[19] = 1;
    constants[20] = 68.55;
    constants[21] = 36.4829;
    constants[22] = 4;
    constants[23] = 0.0002;
    constants[24] = 0.0374842;
    constants[25] = 0.0003;
    constants[26] = 0.45;
    constants[27] = 1e-24;
    constants[28] = 0;
    constants[29] = 0.005884;
    constants[30] = 0.00011768;
    constants[31] = 0.01;
    constants[32] = 2.5;
    constants[33] = 0.000800224;
    constants[34] = 14.3;
    constants[35] = 10;
    constants[36] = 24.7;
    constants[37] = 130;
    constants[38] = 1.8;
    constants[39] = 5.4;
    constants[40] = 2800;
    constants[41] = 0.0003;
    constants[42] = 0.5;
    constants[43] = 0.4;
    constants[44] = 200000;
    constants[45] = 0.0003969;
    constants[46] = 0.0000441;
    constants[47] = 0.815;
    constants[48] = 0.01;
    constants[49] = 0.0003;
    constants[50] = 0.003;
}

void MaleckarModel::computerates(double VOI, double*  __restrict constants, double*  __restrict rates, double*  __restrict states, double*  __restrict algebraic)
{
algebraic[12] =  2000.00*constants[32]*((1.00000 - states[22]) - states[23]) -  666.000*states[23];
rates[23] = algebraic[12];
algebraic[18] =  0.00350000*exp((( - states[0]*states[0])/30.0000)/30.0000)+0.00150000;
algebraic[6] = 1.00000/(1.00000+exp((states[0] - 1.00000)/- 11.0000));
rates[12] = (algebraic[6] - states[12])/algebraic[18];
algebraic[8] = 1.00000/(1.00000+exp(- (states[0]+6.00000)/8.60000));
algebraic[20] = 0.00900000/(1.00000+exp((states[0]+5.00000)/12.0000))+0.000500000;
rates[14] = (algebraic[8] - states[14])/algebraic[20];
algebraic[9] = 1.00000/(1.00000+exp((states[0]+7.50000)/10.0000));
algebraic[21] = 0.590000/(1.00000+exp((states[0]+60.0000)/10.0000))+3.05000;
rates[15] = (algebraic[9] - states[15])/algebraic[21];
algebraic[14] = 1.00000/(1.00000+exp((states[0]+27.1200)/- 8.21000));
algebraic[2] = (states[0]+25.5700)/28.8000;
algebraic[26] =  4.20000e-05*exp( - algebraic[2]*algebraic[2])+2.40000e-05;
rates[3] = (algebraic[14] - states[3])/algebraic[26];
algebraic[3] = 1.00000/(1.00000+exp((states[0]+63.6000)/5.30000));
algebraic[15] = 1.00000/(1.00000+exp((states[0]+35.1000)/3.20000));
algebraic[27] =  0.0300000*algebraic[15]+0.000300000;
rates[4] = (algebraic[3] - states[4])/algebraic[27];
algebraic[28] =  0.120000*algebraic[15]+0.00300000;
rates[5] = (algebraic[3] - states[5])/algebraic[28];
algebraic[4] = 1.00000/(1.00000+exp((states[0]+9.00000)/- 5.80000));
algebraic[16] = (states[0]+35.0000)/30.0000;
algebraic[29] =  0.00270000*exp( - algebraic[16]*algebraic[16])+0.00200000;
rates[7] = (algebraic[4] - states[7])/algebraic[29];
algebraic[5] = 1.00000/(1.00000+exp((states[0]+27.4000)/7.10000));
algebraic[17] = states[0]+40.0000;
algebraic[30] =  0.161000*exp((( - algebraic[17]*algebraic[17])/14.4000)/14.4000)+0.0100000;
rates[8] = (algebraic[5] - states[8])/algebraic[30];
algebraic[31] =  1.33230*exp((( - algebraic[17]*algebraic[17])/14.2000)/14.2000)+0.0626000;
rates[9] = (algebraic[5] - states[9])/algebraic[31];
algebraic[19] = (states[0]+52.4500)/15.8827;
algebraic[32] =  0.0256350*exp( - algebraic[19]*algebraic[19])+0.0141400;
algebraic[7] = 1.00000/(1.00000+exp((states[0]+40.5000)/11.5000));
rates[13] = (algebraic[7] - states[13])/algebraic[32];
algebraic[22] = (states[0] - 20.0000)/20.0000;
algebraic[33] = 0.700000+ 0.400000*exp( - algebraic[22]*algebraic[22]);
algebraic[10] = 1.00000/(1.00000+exp((states[0] - 19.9000)/- 12.7000));
rates[16] = (algebraic[10] - states[16])/algebraic[33];
algebraic[23] = (states[0]+20.1376)/22.1996;
algebraic[34] = 0.0311800+ 0.217180*exp( - algebraic[23]*algebraic[23]);
algebraic[11] = 1.00000/(1.00000+exp((states[0]+15.0000)/- 6.00000));
rates[17] = (algebraic[11] - states[17])/algebraic[34];
algebraic[13] = states[6]/(states[6]+constants[50]);
algebraic[36] =  algebraic[13]*algebraic[13]*algebraic[13]*algebraic[13];
algebraic[25] = states[19]/(states[19]+constants[49]);
algebraic[38] =  algebraic[25]*algebraic[25]*algebraic[25]*algebraic[25];
algebraic[40] =  203.800*(algebraic[38]+algebraic[36]);
rates[28] =  constants[47]*((1.00000 - states[28]) - states[29]) -  algebraic[40]*states[28];
algebraic[42] = 33.9600+ 339.600*algebraic[38];
rates[29] =  algebraic[40]*states[28] -  algebraic[42]*states[29];
algebraic[43] =  (( constants[0]*constants[1])/constants[2])*log(states[10]/states[11]);
algebraic[44] =  constants[12]*states[12]*states[13]*(states[0] - algebraic[43]);
algebraic[45] =  constants[13]*states[14]*states[15]*(states[0] - algebraic[43]);
algebraic[46] = ( constants[14]*pow(states[10]/1.00000, 0.445700)*(states[0] - algebraic[43]))/(1.00000+exp(( 1.50000*((states[0] - algebraic[43])+3.60000)*constants[2])/( constants[0]*constants[1])));
algebraic[48] = 1.00000/(1.00000+exp((states[0]+55.0000)/24.0000));
algebraic[49] =  constants[16]*states[17]*algebraic[48]*(states[0] - algebraic[43]);
algebraic[47] =  constants[15]*states[16]*(states[0] - algebraic[43]);
algebraic[53] = pow(states[2], 1.50000);
algebraic[54] = ( (( (( constants[20]*states[10])/(states[10]+constants[19]))*algebraic[53])/(algebraic[53]+constants[21]))*(states[0]+150.000))/(states[0]+200.000);
algebraic[1] =  floor(VOI/constants[5])*constants[5];
algebraic[24] = (VOI - algebraic[1]>=constants[4]&&VOI - algebraic[1]<=constants[4]+constants[6] ? constants[7] : 0.00000);
rates[11] = - (((algebraic[44]+algebraic[45]+algebraic[46]+algebraic[47]+algebraic[49]) -  2.00000*algebraic[54])+ algebraic[24]*constants[3])/( constants[29]*constants[2]);
rates[10] = (constants[39] - states[10])/constants[35]+((algebraic[44]+algebraic[45]+algebraic[46]+algebraic[47]+algebraic[49]) -  2.00000*algebraic[54])/( constants[33]*constants[2]);
algebraic[35] =  (( constants[0]*constants[1])/constants[2])*log(states[1]/states[2]);
algebraic[37] = ( (( constants[8]*states[3]*states[3]*states[3]*( 0.900000*states[4]+ 0.100000*states[5])*states[1]*states[0]*constants[2]*constants[2])/( constants[0]*constants[1]))*(exp(( (states[0] - algebraic[35])*constants[2])/( constants[0]*constants[1])) - 1.00000))/(exp(( states[0]*constants[2])/( constants[0]*constants[1])) - 1.00000);
algebraic[50] =  constants[17]*(states[0] - algebraic[35]);
algebraic[56] = ( constants[24]*( states[2]*states[2]*states[2]*states[18]*exp(( constants[2]*states[0]*constants[26])/( constants[0]*constants[1])) -  states[1]*states[1]*states[1]*states[19]*exp(( (constants[26] - 1.00000)*states[0]*constants[2])/( constants[0]*constants[1]))))/(1.00000+ constants[25]*( states[1]*states[1]*states[1]*states[19]+ states[2]*states[2]*states[2]*states[18]));
rates[2] = - (algebraic[37]+algebraic[50]+ 3.00000*algebraic[56]+ 3.00000*algebraic[54]+constants[28])/( constants[29]*constants[2]);
algebraic[39] = states[6]/(states[6]+constants[11]);


//const double E_Ca_app = constants[10];
const double E_Ca_app = ((constants[0] * constants[1]) / (2. * constants[2])) * log(states[18] / states[6]); // C->R * C->T / (2 * C->F) * log(S->Ca_c / S->Ca_d)
algebraic[41] =  constants[9]*states[7]*( algebraic[39]*states[8]+ (1.00000 - algebraic[39])*states[9])*(states[0] - E_Ca_app);



algebraic[51] =  (( constants[0]*constants[1])/( 2.00000*constants[2]))*log(states[18]/states[19]);
algebraic[52] =  constants[18]*(states[0] - algebraic[51]);
algebraic[55] = ( constants[22]*states[19])/(states[19]+constants[23]);
rates[18] = (constants[38] - states[18])/constants[36]+((algebraic[41]+algebraic[52]+algebraic[55]) -  2.00000*algebraic[56])/( 2.00000*constants[33]*constants[2]);
rates[1] = (constants[37] - states[1])/constants[34]+(algebraic[37]+algebraic[50]+ 3.00000*algebraic[56]+ 3.00000*algebraic[54]+constants[28])/( constants[33]*constants[2]);
algebraic[58] = ( (states[6] - states[19])*2.00000*constants[30]*constants[2])/constants[31];
rates[6] = - (algebraic[41]+algebraic[58])/( 2.00000*constants[30]*constants[2]);
algebraic[57] =  (10.0000/(1.00000+( 9.13652*pow(1.00000, 0.477811))/pow(constants[27], 0.477811)))*(0.0517000+0.451600/(1.00000+exp((states[0]+59.5300)/17.1800)))*(states[0] - algebraic[43])*constants[3];
algebraic[59] = (algebraic[37]+algebraic[41]+algebraic[44]+algebraic[45]+algebraic[46]+algebraic[49]+algebraic[47]+algebraic[50]+algebraic[52]+algebraic[54]+algebraic[55]+algebraic[56]+algebraic[57])/constants[3]+algebraic[24];
rates[0] =  - algebraic[59]*1000.00;
algebraic[60] =  200000.*states[19]*(1.00000 - states[20]) -  476.000*states[20];
rates[20] = algebraic[60];
algebraic[61] =  78400.0*states[19]*(1.00000 - states[21]) -  392.000*states[21];
rates[21] = algebraic[61];
algebraic[62] =  200000.*states[19]*((1.00000 - states[22]) - states[23]) -  6.60000*states[22];
rates[22] = algebraic[62];
algebraic[63] =  0.0800000*algebraic[61]+ 0.160000*algebraic[62]+ 0.0450000*algebraic[60];
rates[24] = algebraic[63];
algebraic[67] = ( constants[40]*(states[19]/constants[41] - ( constants[43]*constants[43]*states[26])/constants[42]))/((states[19]+constants[41])/constants[41]+( constants[43]*(states[26]+constants[42]))/constants[42]);
algebraic[64] = states[29]/(states[29]+0.250000);
algebraic[65] =  algebraic[64]*algebraic[64];
algebraic[66] =  constants[44]*algebraic[65]*(states[25] - states[19]);
rates[19] = - ((algebraic[52]+algebraic[55]+algebraic[67]) - (algebraic[58]+algebraic[66]+ 2.00000*algebraic[56]))/( 2.00000*constants[29]*constants[2]) -  1.00000*algebraic[63];
algebraic[68] = ( (states[26] - states[25])*2.00000*constants[46]*constants[2])/constants[48];
rates[26] = (algebraic[67] - algebraic[68])/( 2.00000*constants[45]*constants[2]);
algebraic[69] =  480.000*states[25]*(1.00000 - states[27]) -  400.000*states[27];
rates[27] = algebraic[69];
rates[25] = (algebraic[68] - algebraic[66])/( 2.00000*constants[46]*constants[2]) -  31.0000*algebraic[69];
}


/*not maintained at all
void
computeVariables(double VOI, double* constants, double* rates, double* states, double* algebraic)
{
algebraic[12] =  2000.00*constants[32]*((1.00000 - states[22]) - states[23]) -  666.000*states[23];
algebraic[18] =  0.00350000*exp((( - states[0]*states[0])/30.0000)/30.0000)+0.00150000;
algebraic[6] = 1.00000/(1.00000+exp((states[0] - 1.00000)/- 11.0000));
algebraic[8] = 1.00000/(1.00000+exp(- (states[0]+6.00000)/8.60000));
algebraic[20] = 0.00900000/(1.00000+exp((states[0]+5.00000)/12.0000))+0.000500000;
algebraic[9] = 1.00000/(1.00000+exp((states[0]+7.50000)/10.0000));
algebraic[21] = 0.590000/(1.00000+exp((states[0]+60.0000)/10.0000))+3.05000;
algebraic[14] = 1.00000/(1.00000+exp((states[0]+27.1200)/- 8.21000));
algebraic[2] = (states[0]+25.5700)/28.8000;
algebraic[26] =  4.20000e-05*exp( - algebraic[2]*algebraic[2])+2.40000e-05;
algebraic[3] = 1.00000/(1.00000+exp((states[0]+63.6000)/5.30000));
algebraic[15] = 1.00000/(1.00000+exp((states[0]+35.1000)/3.20000));
algebraic[27] =  0.0300000*algebraic[15]+0.000300000;
algebraic[28] =  0.120000*algebraic[15]+0.00300000;
algebraic[4] = 1.00000/(1.00000+exp((states[0]+9.00000)/- 5.80000));
algebraic[16] = (states[0]+35.0000)/30.0000;
algebraic[29] =  0.00270000*exp( - algebraic[16]*algebraic[16])+0.00200000;
algebraic[5] = 1.00000/(1.00000+exp((states[0]+27.4000)/7.10000));
algebraic[17] = states[0]+40.0000;
algebraic[30] =  0.161000*exp((( - algebraic[17]*algebraic[17])/14.4000)/14.4000)+0.0100000;
algebraic[31] =  1.33230*exp((( - algebraic[17]*algebraic[17])/14.2000)/14.2000)+0.0626000;
algebraic[19] = (states[0]+52.4500)/15.8827;
algebraic[32] =  0.0256350*exp( - algebraic[19]*algebraic[19])+0.0141400;
algebraic[7] = 1.00000/(1.00000+exp((states[0]+40.5000)/11.5000));
algebraic[22] = (states[0] - 20.0000)/20.0000;
algebraic[33] = 0.700000+ 0.400000*exp( - algebraic[22]*algebraic[22]);
algebraic[10] = 1.00000/(1.00000+exp((states[0] - 19.9000)/- 12.7000));
algebraic[23] = (states[0]+20.1376)/22.1996;
algebraic[34] = 0.0311800+ 0.217180*exp( - algebraic[23]*algebraic[23]);
algebraic[11] = 1.00000/(1.00000+exp((states[0]+15.0000)/- 6.00000));
algebraic[13] = states[6]/(states[6]+constants[50]);
algebraic[36] =  algebraic[13]*algebraic[13]*algebraic[13]*algebraic[13];
algebraic[25] = states[19]/(states[19]+constants[49]);
algebraic[38] =  algebraic[25]*algebraic[25]*algebraic[25]*algebraic[25];
algebraic[40] =  203.800*(algebraic[38]+algebraic[36]);
algebraic[42] = 33.9600+ 339.600*algebraic[38];
algebraic[43] =  (( constants[0]*constants[1])/constants[2])*log(states[10]/states[11]);
algebraic[44] =  constants[12]*states[12]*states[13]*(states[0] - algebraic[43]);
algebraic[45] =  constants[13]*states[14]*states[15]*(states[0] - algebraic[43]);
algebraic[46] = ( constants[14]*pow(states[10]/1.00000, 0.445700)*(states[0] - algebraic[43]))/(1.00000+exp(( 1.50000*((states[0] - algebraic[43])+3.60000)*constants[2])/( constants[0]*constants[1])));
algebraic[48] = 1.00000/(1.00000+exp((states[0]+55.0000)/24.0000));
algebraic[49] =  constants[16]*states[17]*algebraic[48]*(states[0] - algebraic[43]);
algebraic[47] =  constants[15]*states[16]*(states[0] - algebraic[43]);
algebraic[53] = pow(states[2], 1.50000);
algebraic[54] = ( (( (( constants[20]*states[10])/(states[10]+constants[19]))*algebraic[53])/(algebraic[53]+constants[21]))*(states[0]+150.000))/(states[0]+200.000);
algebraic[1] =  floor(VOI/constants[5])*constants[5];
algebraic[24] = (VOI - algebraic[1]>=constants[4]&&VOI - algebraic[1]<=constants[4]+constants[6] ? constants[7] : 0.00000);
algebraic[35] =  (( constants[0]*constants[1])/constants[2])*log(states[1]/states[2]);
algebraic[37] = ( (( constants[8]*states[3]*states[3]*states[3]*( 0.900000*states[4]+ 0.100000*states[5])*states[1]*states[0]*constants[2]*constants[2])/( constants[0]*constants[1]))*(exp(( (states[0] - algebraic[35])*constants[2])/( constants[0]*constants[1])) - 1.00000))/(exp(( states[0]*constants[2])/( constants[0]*constants[1])) - 1.00000);
algebraic[50] =  constants[17]*(states[0] - algebraic[35]);
algebraic[56] = ( constants[24]*( states[2]*states[2]*states[2]*states[18]*exp(( constants[2]*states[0]*constants[26])/( constants[0]*constants[1])) -  states[1]*states[1]*states[1]*states[19]*exp(( (constants[26] - 1.00000)*states[0]*constants[2])/( constants[0]*constants[1]))))/(1.00000+ constants[25]*( states[1]*states[1]*states[1]*states[19]+ states[2]*states[2]*states[2]*states[18]));
algebraic[39] = states[6]/(states[6]+constants[11]);
algebraic[41] =  constants[9]*states[7]*( algebraic[39]*states[8]+ (1.00000 - algebraic[39])*states[9])*(states[0] - constants[10]);
algebraic[51] =  (( constants[0]*constants[1])/( 2.00000*constants[2]))*log(states[18]/states[19]);
algebraic[52] =  constants[18]*(states[0] - algebraic[51]);
algebraic[55] = ( constants[22]*states[19])/(states[19]+constants[23]);
algebraic[58] = ( (states[6] - states[19])*2.00000*constants[30]*constants[2])/constants[31];
algebraic[57] =  (10.0000/(1.00000+( 9.13652*pow(1.00000, 0.477811))/pow(constants[27], 0.477811)))*(0.0517000+0.451600/(1.00000+exp((states[0]+59.5300)/17.1800)))*(states[0] - algebraic[43])*constants[3];
algebraic[59] = (algebraic[37]+algebraic[41]+algebraic[44]+algebraic[45]+algebraic[46]+algebraic[49]+algebraic[47]+algebraic[50]+algebraic[52]+algebraic[54]+algebraic[55]+algebraic[56]+algebraic[57])/constants[3]+algebraic[24];
algebraic[60] =  200000.*states[19]*(1.00000 - states[20]) -  476.000*states[20];
algebraic[61] =  78400.0*states[19]*(1.00000 - states[21]) -  392.000*states[21];
algebraic[62] =  200000.*states[19]*((1.00000 - states[22]) - states[23]) -  6.60000*states[22];
algebraic[63] =  0.0800000*algebraic[61]+ 0.160000*algebraic[62]+ 0.0450000*algebraic[60];
algebraic[67] = ( constants[40]*(states[19]/constants[41] - ( constants[43]*constants[43]*states[26])/constants[42]))/((states[19]+constants[41])/constants[41]+( constants[43]*(states[26]+constants[42]))/constants[42]);
algebraic[64] = states[29]/(states[29]+0.250000);
algebraic[65] =  algebraic[64]*algebraic[64];
algebraic[66] =  constants[44]*algebraic[65]*(states[25] - states[19]);
algebraic[68] = ( (states[26] - states[25])*2.00000*constants[46]*constants[2])/constants[48];
algebraic[69] =  480.000*states[25]*(1.00000 - states[27]) -  400.000*states[27];
algebraic[0] =  0.0500000*states[0];
}
*/
