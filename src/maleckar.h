//
// Created by andrey on 16/03/20.
//

#ifndef GA_MALECKAR_H

#include <cmath>
#include <iostream>

const int STATE_ARRAY_SIZE = 39;

struct State {
	double V; // in component membrane (millivolt).
	double Na_c; // in component cleft_space_ion_concentrations (millimolar).
	double Na_i; // in component intracellular_ion_concentrations (millimolar).
	double m; // in component sodium_current_m_gate (dimensionless).
	double h1; // in component sodium_current_h1_gate (dimensionless).
	double h2; // in component sodium_current_h2_gate (dimensionless).
	double Ca_d; // in component intracellular_ion_concentrations (millimolar).
	double d_L; // in component L_type_Ca_channel_d_L_gate (dimensionless).
	double f_L1; // in component L_type_Ca_channel_f_L1_gate (dimensionless).
	double f_L2; // in component L_type_Ca_channel_f_L2_gate (dimensionless).
	double K_c; // in component cleft_space_ion_concentrations (millimolar).
	double K_i; // in component intracellular_ion_concentrations (millimolar).
	double r; // in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
	double s; // in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
	double a_ur; // in component ultra_rapid_K_current_aur_gate (dimensionless).
	double i_ur; // in component ultra_rapid_K_current_iur_gate (dimensionless).
	double n; // in component delayed_rectifier_K_currents_n_gate (dimensionless).
	double pa; // in component delayed_rectifier_K_currents_pa_gate (dimensionless).
	double Ca_c; // in component cleft_space_ion_concentrations (millimolar).
	double Ca_i; // in component intracellular_ion_concentrations (millimolar).
	double O_C; // in component intracellular_Ca_buffering (dimensionless).
	double O_TC; // in component intracellular_Ca_buffering (dimensionless).
	double O_TMgC; // in component intracellular_Ca_buffering (dimensionless).
	double O_TMgMg; // in component intracellular_Ca_buffering (dimensionless).
	double O; // in component intracellular_Ca_buffering (dimensionless).
	double Ca_rel; // in component Ca_handling_by_the_SR (millimolar).
	double Ca_up; // in component Ca_handling_by_the_SR (millimolar).
	double O_Calse; // in component Ca_handling_by_the_SR (dimensionless).
	double F1; // in component Ca_handling_by_the_SR (dimensionless).
	double F2; // in component Ca_handling_by_the_SR (dimensionless).

    double d_ord; // ORd
    double ff;
    double fs;
    double fcaf;
    double fcas;
    double jca;
    double ffp;
    double fcafp;
    double nca;
};

struct Algebraic {
	double Q_tot; // in component membrane (millivolt).
	double past; // in component membrane (second).
	double m_factor; // in component sodium_current_m_gate (dimensionless).
	double h_infinity; // in component sodium_current_h1_gate (dimensionless).
	double d_L_infinity; // in component L_type_Ca_channel_d_L_gate (dimensionless).
	double f_L_infinity; // in component L_type_Ca_channel_f_L1_gate (dimensionless).
	double r_infinity; // in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
	double s_infinity; // in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
	double a_ur_infinity; // in component ultra_rapid_K_current_aur_gate (dimensionless).
	double i_ur_infinity; // in component ultra_rapid_K_current_iur_gate (dimensionless).
	double n_infinity; // in component delayed_rectifier_K_currents_n_gate (dimensionless).
	double p_a_infinity; // in component delayed_rectifier_K_currents_pa_gate (dimensionless).
	double J_O_TMgMg; // in component intracellular_Ca_buffering (per_second).
	double r_Ca_d_term; // in component Ca_handling_by_the_SR (dimensionless).
	double m_infinity; // in component sodium_current_m_gate (dimensionless).
	double h_factor; // in component sodium_current_h1_gate (dimensionless).
	double d_L_factor; // in component L_type_Ca_channel_d_L_gate (dimensionless).
	double f_L_factor; // in component L_type_Ca_channel_f_L1_gate (millivolt).
	double tau_r; // in component Ca_independent_transient_outward_K_current_r_gate (second).
	double s_factor; // in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
	double tau_a_ur; // in component ultra_rapid_K_current_aur_gate (second).
	double tau_i_ur; // in component ultra_rapid_K_current_iur_gate (second).
	double n_factor; // in component delayed_rectifier_K_currents_n_gate (dimensionless).
	double pa_factor; // in component delayed_rectifier_K_currents_pa_gate (dimensionless).
	double i_Stim; // in component membrane (pA_per_nF).
	double r_Ca_i_term; // in component Ca_handling_by_the_SR (dimensionless).
	double tau_m; // in component sodium_current_m_gate (second).
	double tau_h1; // in component sodium_current_h1_gate (second).
	double tau_h2; // in component sodium_current_h2_gate (second).
	double tau_d_L; // in component L_type_Ca_channel_d_L_gate (second).
	double tau_f_L1; // in component L_type_Ca_channel_f_L1_gate (second).
	double tau_f_L2; // in component L_type_Ca_channel_f_L2_gate (second).
	double tau_s; // in component Ca_independent_transient_outward_K_current_s_gate (second).
	double tau_n; // in component delayed_rectifier_K_currents_n_gate (second).
	double tau_pa; // in component delayed_rectifier_K_currents_pa_gate (second).
	double E_Na; // in component sodium_current (millivolt).
	double r_Ca_d_factor; // in component Ca_handling_by_the_SR (dimensionless).
	double i_Na; // in component sodium_current (picoA).
	double r_Ca_i_factor; // in component Ca_handling_by_the_SR (dimensionless).
	double f_Ca; // in component L_type_Ca_channel (dimensionless).
	double r_act; // in component Ca_handling_by_the_SR (per_second).
	double i_Ca_L; // in component L_type_Ca_channel (picoA).
	double r_inact; // in component Ca_handling_by_the_SR (per_second).
	double E_K; // in component Ca_independent_transient_outward_K_current (millivolt).
	double i_t; // in component Ca_independent_transient_outward_K_current (picoA).
	double i_Kur; // in component ultra_rapid_K_current (picoA).
	double i_K1; // in component inward_rectifier (picoA).
	double i_Ks; // in component delayed_rectifier_K_currents (picoA).
	double pip; // in component delayed_rectifier_K_currents_pi_gate (dimensionless).
	double i_Kr; // in component delayed_rectifier_K_currents (picoA).
	double i_B_Na; // in component background_currents (picoA).
	double E_Ca; // in component background_currents (millivolt).
	double i_B_Ca; // in component background_currents (picoA).
	double pow_Na_i_15; // in component sodium_potassium_pump (millimolar15).
	double i_NaK; // in component sodium_potassium_pump (picoA).
	double i_CaP; // in component sarcolemmal_calcium_pump_current (picoA).
	double i_NaCa; // in component Na_Ca_ion_exchanger_current (picoA).
	double i_KACh; // in component ACh_dependent_K_current (picoA).
	double i_di; // in component intracellular_ion_concentrations (picoA).
	double I; // in component membrane (pA_per_nF).
	double J_O_C; // in component intracellular_Ca_buffering (per_second).
	double J_O_TC; // in component intracellular_Ca_buffering (per_second).
	double J_O_TMgC; // in component intracellular_Ca_buffering (per_second).
	double J_O; // in component intracellular_Ca_buffering (per_second).
	double i_rel_f2; // in component Ca_handling_by_the_SR (dimensionless).
	double i_rel_factor; // in component Ca_handling_by_the_SR (dimensionless).
	double i_rel; // in component Ca_handling_by_the_SR (picoA).
	double i_up; // in component Ca_handling_by_the_SR (picoA).
	double i_tr; // in component Ca_handling_by_the_SR (picoA).
	double J_O_Calse; // in component Ca_handling_by_the_SR (per_second).

    double i_Ca_L_ToR_ORd; // (A_per_F)
};

struct Constants {
	double R; // in component membrane (millijoule_per_mole_kelvin).
	double T; // in component membrane (kelvin).
	double F; // in component membrane (coulomb_per_mole).
	double Cm; // in component membrane (nanoF).
	double stim_offset; // in component membrane (second).
	double stim_period; // in component membrane (second).
	double stim_duration; // in component membrane (second).
	double P_Na; // in component sodium_current (nanolitre_per_second).
	double g_Ca_L; // in component L_type_Ca_channel (nanoS).
	double E_Ca_app; // in component L_type_Ca_channel (millivolt).
	double k_Ca; // in component L_type_Ca_channel (millimolar).
	double g_t; // in component Ca_independent_transient_outward_K_current (nanoS).
	double g_kur; // in component ultra_rapid_K_current (nanoS).
	double g_K1; // in component inward_rectifier (nanoS).
	double g_Ks; // in component delayed_rectifier_K_currents (nanoS).
	double g_Kr; // in component delayed_rectifier_K_currents (nanoS).
	double g_B_Na; // in component background_currents (nanoS).
	double g_B_Ca; // in component background_currents (nanoS).
	double K_NaK_K; // in component sodium_potassium_pump (millimolar).
	double i_NaK_max; // in component sodium_potassium_pump (picoA).
	double pow_K_NaK_Na_15; // in component sodium_potassium_pump (millimolar15).
	double i_CaP_max; // in component sarcolemmal_calcium_pump_current (picoA).
	double k_CaP; // in component sarcolemmal_calcium_pump_current (millimolar).
	double K_NaCa; // in component Na_Ca_ion_exchanger_current (picoA_per_millimolar_4).
	double d_NaCa; // in component Na_Ca_ion_exchanger_current (per_millimolar_4).
	double gamma_Na; // in component Na_Ca_ion_exchanger_current (dimensionless).
	double ACh; // in component ACh_dependent_K_current (millimolar).
	double phi_Na_en; // in component intracellular_ion_concentrations (picoA).
	double Vol_i; // in component intracellular_ion_concentrations (nanolitre).
	double Vol_d; // in component intracellular_ion_concentrations (nanolitre).
	double tau_di; // in component intracellular_ion_concentrations (second).
	double Mg_i; // in component intracellular_Ca_buffering (millimolar).
	double Vol_c; // in component cleft_space_ion_concentrations (nanolitre).
	double tau_Na; // in component cleft_space_ion_concentrations (second).
	double tau_K; // in component cleft_space_ion_concentrations (second).
	double tau_Ca; // in component cleft_space_ion_concentrations (second).
	double Na_b; // in component cleft_space_ion_concentrations (millimolar).
	double Ca_b; // in component cleft_space_ion_concentrations (millimolar).
	double K_b; // in component cleft_space_ion_concentrations (millimolar).
	double I_up_max; // in component Ca_handling_by_the_SR (picoA).
	double k_cyca; // in component Ca_handling_by_the_SR (millimolar).
	double k_srca; // in component Ca_handling_by_the_SR (millimolar).
	double k_xcs; // in component Ca_handling_by_the_SR (dimensionless).
	double alpha_rel; // in component Ca_handling_by_the_SR (picoA_per_millimolar).
	double Vol_up; // in component Ca_handling_by_the_SR (nanolitre).
	double Vol_rel; // in component Ca_handling_by_the_SR (nanolitre).
	double r_recov; // in component Ca_handling_by_the_SR (per_second).
	double tau_tr; // in component Ca_handling_by_the_SR (second).
	double k_rel_i; // in component Ca_handling_by_the_SR (millimolar).
	double k_rel_d; // in component Ca_handling_by_the_SR (millimolar).

    double pca_ord; // for ical_ord()
};

int action_potential(struct State *initial_state, double *scaling_coefficients, double *AP, float CL, float amp,
                     int current_time, int iso, int baseline_index, int amount_of_baselines, int amount_of_genes);

#define GA_MALECKAR_H

#endif //GA_MALECKAR_H
