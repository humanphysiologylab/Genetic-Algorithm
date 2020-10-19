#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "maleckar.h"


void initialize_states_default(State *S) {
	S->V = -74.031982;
	S->Na_c = 130.022096;
	S->Na_i = 8.516766;
	S->m = 0.003289;
	S->h1 = 0.877202;
	S->h2 = 0.873881;
	S->Ca_d = 7.1e-5;
	S->d_L = 0.000014;
	S->f_L1 = 0.998597;
	S->f_L2 = 0.998586;
	S->K_c = 5.560224;
	S->K_i = 129.485991;
	S->r = 0.001089;
	S->s = 0.948597;
	S->a_ur = 0.000367;
	S->i_ur = 0.96729;
	S->n = 0.004374;
	S->pa = 0.000053;
	S->Ca_c = 1.815768;
	S->Ca_i = 6.5e-5;
	S->O_C = 0.026766;
	S->O_TC = 0.012922;
	S->O_TMgC = 0.190369;
	S->O_TMgMg = 0.714463;
	S->O = 1.38222;
	S->Ca_rel = 0.632613;
	S->Ca_up = 0.649195;
	S->O_Calse = 0.431547;
	S->F1 = 0.470055;
	S->F2 = 0.002814;

    S->d_ord = 0.;
    S->ff = 1.;
    S->fs = 1.;
    S->fcaf = 1.;
    S->fcas = 1.;
    S->jca = 1.;
    S->ffp = 1.;
    S->fcafp = 1.;
    S->nca = 0.;
}


void initialize_states(State *S, State *initial_state) {
	S->V = initial_state->V;
	S->Na_c = initial_state->Na_c;
	S->Na_i = initial_state->Na_i;
	S->m = initial_state->m;
	S->h1 = initial_state->h1;
	S->h2 = initial_state->h2;
	S->Ca_d = initial_state->Ca_d;
	S->d_L = initial_state->d_L;
	S->f_L1 = initial_state->f_L1;
	S->f_L2 = initial_state->f_L2;
	S->K_c = initial_state->K_c;
	S->K_i = initial_state->K_i;
	S->r = initial_state->r;
	S->s = initial_state->s;
	S->a_ur = initial_state->a_ur;
	S->i_ur = initial_state->i_ur;
	S->n = initial_state->n;
	S->pa = initial_state->pa;
	S->Ca_c = initial_state->Ca_c;
	S->Ca_i = initial_state->Ca_i;
	S->O_C = initial_state->O_C;
	S->O_TC = initial_state->O_TC;
	S->O_TMgC = initial_state->O_TMgC;
	S->O_TMgMg = initial_state->O_TMgMg;
	S->O = initial_state->O;
	S->Ca_rel = initial_state->Ca_rel;
	S->Ca_up = initial_state->Ca_up;
	S->O_Calse = initial_state->O_Calse;
	S->F1 = initial_state->F1;
	S->F2 = initial_state->F2;

    S->d_ord = initial_state->d_ord;
    S->ff = initial_state->ff;
    S->fs = initial_state->fs;
    S->fcaf = initial_state->fcaf;
    S->fcas = initial_state->fcas;
    S->jca = initial_state->jca;
    S->ffp = initial_state->ffp;
    S->fcafp = initial_state->fcafp;
    S->nca = initial_state->nca;
}

void initialize_constants_default(Constants *C, double CL, double *scaling_coefficients) {

    double *sc = scaling_coefficients;

	C->R = 8314.;
	C->T = 306.15;
	C->F = 96487.;
	C->Cm = 50.;
	C->stim_offset = 0.;
	C->stim_period = CL; // sec
	C->stim_duration = 0.001;
	C->P_Na = 0.0018 * sc[0];
	C->g_Ca_L = 6.75 * sc[1];
	C->E_Ca_app = 60.;
	C->k_Ca = 0.025;
	C->g_t = 8.25 * sc[2];
	C->g_kur = 2.25 * sc[3];
	C->g_K1 = 3.1 * sc[4];
	C->g_Ks = 1. * sc[5];
	C->g_Kr = 0.5 * sc[6];
	C->g_B_Na = 0.060599 * sc[7];
	C->g_B_Ca = 0.078681 * sc[8];
	C->K_NaK_K = 1.;
	C->i_NaK_max = 68.55 * sc[9];
	C->pow_K_NaK_Na_15 = 36.4829;
	C->i_CaP_max = 4. * sc[10];
	C->k_CaP = 0.0002;
	C->K_NaCa = 0.0374842 * sc[11];
	C->d_NaCa = 0.0003;
	C->gamma_Na = 0.45;
	C->ACh = 1e-24;
	C->phi_Na_en = 0.;
	C->Vol_i = 0.005884;
	C->Vol_d = 0.00011768;
	C->tau_di = 0.01;
	C->Mg_i = 2.5;
	C->Vol_c = 0.000800224;
	C->tau_Na = 14.3;
	C->tau_K = 10.;
	C->tau_Ca = 24.7;
	C->Na_b = 140.;
	C->Ca_b = 2.;
	C->K_b = 4.;
	C->I_up_max = 2800. * sc[12];
	C->k_cyca = 0.0003;
	C->k_srca = 0.5;
	C->k_xcs = 0.4;
	C->alpha_rel = 200000. * sc[13];
	C->Vol_up = 0.0003969;
	C->Vol_rel = 0.0000441;
	C->r_recov = 0.815;
	C->tau_tr = 0.01;
	C->k_rel_i = 0.0003;
	C->k_rel_d = 0.003;

	C->pca_ord = 5e-5 * sc[14];
}


double calculate_ICaL(State *S, State *R, Algebraic *A, Constants *C, const bool GHK_FLAG = false) {

    if (GHK_FLAG) {
        A->d_L_infinity = (S->V <= 50.0) ? 1.9 * exp(-2.47567057 * exp(-0.02690762 * S->V)) : 1.;
    } else {
        A->d_L_infinity = 1.0 / (1.0 + exp((S->V + 9.0) / -5.8));
    }
    A->tau_d_L = 0.0027 * exp(-pow((S->V + 35.0) / 30.0, 2)) + 0.002;
    R->d_L = (A->d_L_infinity - S->d_L) / A->tau_d_L;

    A->f_L_infinity = 1.0 / (1.0 + exp((S->V + 27.4) / 7.1));
    A->tau_f_L1 = 0.161  * exp(-pow((S->V + 40.0) / 14.4, 2)) + 0.01;
    R->f_L1 = (A->f_L_infinity - S->f_L1) / A->tau_f_L1;

    A->tau_f_L2 = 1.3323 * exp(-pow((S->V + 40.0) / 14.2, 2)) + 0.0626;
    R->f_L2 = (A->f_L_infinity - S->f_L2) / A->tau_f_L2;

    A->f_Ca = S->Ca_d / (S->Ca_d + C->k_Ca);

    double ICaL_GHK = 0.;
    if (GHK_FLAG) {
        double gamma_ca_i = 0.61, gamma_ca_o = 0.66;
        double phi_CaL = 4.0 * (S->V * C->F * C->F / (C->R * C->T)) *
                         (gamma_ca_i * S->Ca_d * exp(2.0 * (S->V * C->F / (C->R * C->T))) - gamma_ca_o * S->Ca_c) /
                         (exp(2.0 * (S->V * C->F / (C->R * C->T))) - 1.0);
        ICaL_GHK = 0.004 * (C->g_Ca_L / 6.75) * phi_CaL * S->d_L * (A->f_Ca * S->f_L1 + (1 - A->f_Ca) * S->f_L2);
    } else {
        ICaL_GHK = C->g_Ca_L * S->d_L * (A->f_Ca * S->f_L1 + (1 - A->f_Ca) * S->f_L2) * (S->V - C->E_Ca_app);
    }

    return ICaL_GHK;
}


inline double calculate_dss(double v) {
    //return 1.0 / (1.0 + exp((-(v + 3.940)) / 4.230)); // default ORd model
    return (v < 31.4978) ? 1.0763 * exp(-1.0070 * exp(-0.0829 * v)) : 1; // ToR_ORd
}

inline double calculate_fss(double v) {
    return 1.0 / (1.0 + exp((v + 19.58) / 3.696));
}

inline double calculate_tff(double v) {
    return 7.0 + 1.0 / (0.0045 * exp(-(v + 20.0) / 10.0) +
                        0.0045 * exp((v + 20.0) / 10.0));
}

inline double calculate_tfs(double v) {
    return 1000.0 + 1.0 / (0.000035 * exp(-(v + 5.0) / 4.0) +
                           0.000035 * exp((v + 5.0) / 6.0));
}

inline double calculate_tfcaf(double v) {
    return 7.0 + 1.0 / (0.04 * exp(-(v - 4.0) / 7.0) +
                        0.04 * exp((v - 4.0) / 7.0));
}

inline double calculate_tfcas(double v) {
    return 100.0 + 1.0 / (0.00012 * exp(-v / 3.0) + 0.00012 * exp(v / 7.0));
}

inline double calculate_Afcaf(double v) {
    return 0.3 + 0.6 / (1.0 + exp((v - 10.0) / 10.0));
}

inline double calculate_td(double v) {
    return 0.6 + 1.0 / (exp(-0.05 * (v + 6.0)) + exp(0.09 * (v + 14.0)));
}

double calculate_ICaL_ToR_ORd(State *S, State *R, Algebraic *A, Constants *C) {

    double dss = calculate_dss(S->V); // dimensionless
    double fss = calculate_fss(S->V); // dimensionless
    double td = calculate_td(S->V) /*ms*/ / 1000.; // sec
    double tff = calculate_tff(S->V) /*ms*/ / 1000.; // sec
    double tfs = calculate_tfs(S->V) /*ms*/ / 1000.; // sec
    double tfcaf = calculate_tfcaf(S->V) /*ms*/ / 1000.; // sec
    double tfcas = calculate_tfcas(S->V) /*ms*/ / 1000.; // sec
    double Afcaf = calculate_Afcaf(S->V); // dimensionless

    double Aff = 0.6;
    double Afs = 1.0 - Aff;
    double f = Aff * S->ff + Afs * S->fs;
    double fcass = fss;
    double Afcas = 1.0 - Afcaf;
    double fca = Afcaf * S->fcaf + Afcas * S->fcas;
    double tjca = 75.0 /*ms*/ / 1000.; // sec
    double tffp = 2.5 * tff;
    double tfcafp = 2.5 * tfcaf;
    double Kmn = 0.002; // mM
    double k2n = 1000.0  /*1/ms*/ * 1000.; // sec
    double km2n = S->jca * 1.0  /1/*ms*/ * 1000.; // sec
    double anca = 1.0 / (k2n / km2n + pow(1.0 + Kmn / S->Ca_d, 4.0));

    double gamma_ca_i = 0.61, gamma_ca_o = 0.66;
    double phi_CaL = 4.0 * (S->V * C->F * C->F / (C->R * C->T)) *
                     (gamma_ca_i * S->Ca_d * exp(2.0 * (S->V * C->F / (C->R * C->T))) - gamma_ca_o * S->Ca_c) /
                     (exp(2.0 * (S->V * C->F / (C->R * C->T))) - 1.0);

    // Rush-Larsen method was used in ORd
    // Here I use simple Euler (like everywhere else)
    R->d_ord = (dss - S->d_ord) / td; // NextSt->d = Sc->dss - (Sc->dss - S->d_ord) * exp(-Par->dt / td);
    R->ff = (fss - S->ff) / tff; // NextSt->ff = Sc->fss - (Sc->fss - S->ff) * exp(-Par->dt / tff);
    R->fs = (fss - S->fs) / tfs; // NextSt->fs = Sc->fss - (Sc->fss - S->fs) * exp(-Par->dt / tfs);
    R->fcaf = (fcass - S->fcaf) / tfcaf; // NextSt->fcaf = fcass - (fcass - S->fcaf) * exp(-Par->dt / tfcaf);
    R->fcas = (fcass - S->fcas) / tfcas; // NextSt->fcas = fcass - (fcass - S->fcas) * exp(-Par->dt / tfcas);
    R->jca = (fcass - S->jca) / tjca; // NextSt->jca = fcass - (fcass - S->jca) * exp(-Par->dt / tjca);
    R->ffp = (fss - S->ffp) / tffp; // NextSt->ffp = Sc->fss - (Sc->fss - S->ffp) * exp(-Par->dt / tffp);
    R->fcafp = (fcass - S->fcafp) / tfcafp; // NextSt->fcafp = fcass - (fcass - S->fcafp) * exp(-Par->dt / tfcafp);
    R->nca = anca * k2n - S->nca * km2n; // NextSt->nca = anca * k2n / km2n - (anca * k2n / km2n - S->nca) * exp(-km2n * Par->dt);

    double ICaL_ORd_max = C->pca_ord * phi_CaL; // Col / m^3
    double ICaL_ORd = ICaL_ORd_max * S->d_ord * (f  * (1.0 - S->nca) + S->jca * fca  * S->nca);

    ICaL_ORd = ICaL_ORd * C->Cm /*nF*/; // in ORd model [ICaL] = [uA/uF]

    return ICaL_ORd;
}


void compute_rates(double VOI, State *S, State *R, Algebraic *A, Constants *C, double dt, double stim_amplitude) {

    const bool IMPLICIT_FLAG = true;
    // A->J_X = kon * C->Mg_ww * ((1.0 - S->Y) - S->X) - koff * S->X;
    // A->J_X = ((S->X + dt * (kon* C->Mg_ww * (1.0 - S->Y))) / (1 + dt * (kon* C->Mg_ww + koff)) - S->X) / dt;

    if (IMPLICIT_FLAG) {
        A->J_O_TMgMg = ((S->O_TMgMg + dt * (2000.0 * C->Mg_i * (1.0 - S->O_TMgC))) / (1. + dt * (2000.0 * C->Mg_i + 666.0)) - S->O_TMgMg) / dt;
    } else {
        A->J_O_TMgMg = 2000.00 * C->Mg_i * ((1.00000 - S->O_TMgC) - S->O_TMgMg) - 666.000 * S->O_TMgMg;
    }
	R->O_TMgMg = A->J_O_TMgMg;

	A->tau_r =  0.00350000*exp((( - S->V*S->V)/30.0000)/30.0000)+0.00150000;
	A->r_infinity = 1.00000/(1.00000+exp((S->V - 1.00000)/- 11.0000));
	R->r = (A->r_infinity - S->r)/A->tau_r;
	A->a_ur_infinity = 1.00000/(1.00000+exp(- (S->V+6.00000)/8.60000));
	A->tau_a_ur = 0.00900000/(1.00000+exp((S->V+5.00000)/12.0000))+0.000500000;
	R->a_ur = (A->a_ur_infinity - S->a_ur)/A->tau_a_ur;
	A->i_ur_infinity = 1.00000/(1.00000+exp((S->V+7.50000)/10.0000));
	A->tau_i_ur = 0.590000/(1.00000+exp((S->V+60.0000)/10.0000))+3.05000;
	R->i_ur = (A->i_ur_infinity - S->i_ur)/A->tau_i_ur;
	A->m_infinity = 1.00000/(1.00000+exp((S->V+27.1200)/- 8.21000));
	A->m_factor = (S->V+25.5700)/28.8000;
	A->tau_m =  4.20000e-05*exp( - A->m_factor*A->m_factor)+2.40000e-05;
	R->m = (A->m_infinity - S->m)/A->tau_m;
	A->h_infinity = 1.00000/(1.00000+exp((S->V+63.6000)/5.30000));
	A->h_factor = 1.00000/(1.00000+exp((S->V+35.1000)/3.20000));
	A->tau_h1 =  0.0300000*A->h_factor+0.000300000;
	R->h1 = (A->h_infinity - S->h1)/A->tau_h1;
	A->tau_h2 =  0.120000*A->h_factor+0.00300000;
	R->h2 = (A->h_infinity - S->h2)/A->tau_h2;
	A->s_factor = (S->V+52.4500)/15.8827;
	A->tau_s =  0.0256350*exp( - A->s_factor*A->s_factor)+0.0141400;
	A->s_infinity = 1.00000/(1.00000+exp((S->V+40.5000)/11.5000));
	R->s = (A->s_infinity - S->s)/A->tau_s;
	A->n_factor = (S->V - 20.0000)/20.0000;
	A->tau_n = 0.700000+ 0.400000*exp( - A->n_factor*A->n_factor);
	A->n_infinity = 1.00000/(1.00000+exp((S->V - 19.9000)/- 12.7000));
	R->n = (A->n_infinity - S->n)/A->tau_n;
	A->pa_factor = (S->V+20.1376)/22.1996;
	A->tau_pa = 0.0311800+ 0.217180*exp( - A->pa_factor*A->pa_factor);
	A->p_a_infinity = 1.00000/(1.00000+exp((S->V+15.0000)/- 6.00000));
	R->pa = (A->p_a_infinity - S->pa)/A->tau_pa;
	A->r_Ca_d_term = S->Ca_d/(S->Ca_d+C->k_rel_d);
	A->r_Ca_d_factor =  A->r_Ca_d_term*A->r_Ca_d_term*A->r_Ca_d_term*A->r_Ca_d_term;
	A->r_Ca_i_term = S->Ca_i/(S->Ca_i+C->k_rel_i);
	A->r_Ca_i_factor =  A->r_Ca_i_term*A->r_Ca_i_term*A->r_Ca_i_term*A->r_Ca_i_term;
	A->r_act =  203.800*(A->r_Ca_i_factor+A->r_Ca_d_factor);
	A->r_inact = 33.9600+ 339.600*A->r_Ca_i_factor;

    if (IMPLICIT_FLAG) {
        R->F1 = ((S->F1 + dt * (C->r_recov * (1.0 - S->F2))) / (1. + dt * (C->r_recov + A->r_act)) - S->F1) / dt;
        R->F2 = ((S->F2 + dt * A->r_act * S->F1) / (1. + dt * A->r_inact) - S->F2) / dt;
    } else {
        R->F1 = C->r_recov * ((1.00000 - S->F1) - S->F2) - A->r_act * S->F1;
        R->F2 = A->r_act * S->F1 - A->r_inact * S->F2;
    }

	A->E_K =  (( C->R*C->T)/C->F)*log(S->K_c/S->K_i);
	A->i_t =  C->g_t*S->r*S->s*(S->V - A->E_K);
	A->i_Kur =  C->g_kur*S->a_ur*S->i_ur*(S->V - A->E_K);
	A->i_K1 = ( C->g_K1*pow(S->K_c/1.00000, 0.445700)*(S->V - A->E_K))/(1.00000+exp(( 1.50000*((S->V - A->E_K)+3.60000)*C->F)/( C->R*C->T)));
	A->pip = 1.00000/(1.00000+exp((S->V+55.0000)/24.0000));
	A->i_Kr =  C->g_Kr*S->pa*A->pip*(S->V - A->E_K);
	A->i_Ks =  C->g_Ks*S->n*(S->V - A->E_K);
	A->pow_Na_i_15 = pow(S->Na_i, 1.50000);
	A->i_NaK = ( (( (( C->i_NaK_max*S->K_c)/(S->K_c+C->K_NaK_K))*A->pow_Na_i_15)/(A->pow_Na_i_15+C->pow_K_NaK_Na_15))*(S->V+150.000))/(S->V+200.000);
	A->past =  floor(VOI/C->stim_period)*C->stim_period;
	A->i_Stim = (VOI - A->past>=C->stim_offset&&VOI - A->past<=C->stim_offset+C->stim_duration ? stim_amplitude : 0.00000);
	R->K_i = - (((A->i_t+A->i_Kur+A->i_K1+A->i_Ks+A->i_Kr) -  2.00000*A->i_NaK)+ A->i_Stim*C->Cm)/( C->Vol_i*C->F);
	R->K_c = (C->K_b - S->K_c)/C->tau_K+((A->i_t+A->i_Kur+A->i_K1+A->i_Ks+A->i_Kr) -  2.00000*A->i_NaK)/( C->Vol_c*C->F);
	A->E_Na =  (( C->R*C->T)/C->F)*log(S->Na_c/S->Na_i);
	A->i_Na = ( (( C->P_Na*S->m*S->m*S->m*( 0.900000*S->h1+ 0.100000*S->h2)*S->Na_c*S->V*C->F*C->F)/( C->R*C->T))*(exp(( (S->V - A->E_Na)*C->F)/( C->R*C->T)) - 1.00000))/(exp(( S->V*C->F)/( C->R*C->T)) - 1.00000);
	A->i_B_Na =  C->g_B_Na*(S->V - A->E_Na);
	A->i_NaCa = ( C->K_NaCa*( S->Na_i*S->Na_i*S->Na_i*S->Ca_c*exp(( C->F*S->V*C->gamma_Na)/( C->R*C->T)) -  S->Na_c*S->Na_c*S->Na_c*S->Ca_i*exp(( (C->gamma_Na - 1.00000)*S->V*C->F)/( C->R*C->T))))/(1.00000+ C->d_NaCa*( S->Na_c*S->Na_c*S->Na_c*S->Ca_i+ S->Na_i*S->Na_i*S->Na_i*S->Ca_c));
	R->Na_i = - (A->i_Na+A->i_B_Na+ 3.00000*A->i_NaCa+ 3.00000*A->i_NaK+C->phi_Na_en)/( C->Vol_i*C->F);

	A->i_Ca_L = calculate_ICaL(S, R, A, C, true);
    A->i_Ca_L_ToR_ORd = calculate_ICaL_ToR_ORd(S, R, A, C);

	A->E_Ca =  (( C->R*C->T)/( 2.00000*C->F))*log(S->Ca_c/S->Ca_i);
	A->i_B_Ca =  C->g_B_Ca*(S->V - A->E_Ca);
	A->i_CaP = ( C->i_CaP_max*S->Ca_i)/(S->Ca_i+C->k_CaP);
	R->Ca_c = (C->Ca_b - S->Ca_c)/C->tau_Ca+((A->i_Ca_L+A->i_Ca_L_ToR_ORd+A->i_B_Ca+A->i_CaP) -  2.00000*A->i_NaCa)/( 2.00000*C->Vol_c*C->F);
	R->Na_c = (C->Na_b - S->Na_c)/C->tau_Na+(A->i_Na+A->i_B_Na+ 3.00000*A->i_NaCa+ 3.00000*A->i_NaK+C->phi_Na_en)/( C->Vol_c*C->F);
	A->i_di = ( (S->Ca_d - S->Ca_i)*2.00000*C->Vol_d*C->F)/C->tau_di;
	R->Ca_d = - (A->i_Ca_L+A->i_Ca_L_ToR_ORd+A->i_di)/( 2.00000*C->Vol_d*C->F);
	A->i_KACh =  (10.0000/(1.00000+( 9.13652*pow(1.00000, 0.477811))/pow(C->ACh, 0.477811)))*(0.0517000+0.451600/(1.00000+exp((S->V+59.5300)/17.1800)))*(S->V - A->E_K)*C->Cm;
	A->I = (A->i_Na+A->i_Ca_L+A->i_Ca_L_ToR_ORd+A->i_t+A->i_Kur+A->i_K1+A->i_Kr+A->i_Ks+A->i_B_Na+A->i_B_Ca+A->i_NaK+A->i_CaP+A->i_NaCa+A->i_KACh)/C->Cm+A->i_Stim;
	R->V =  - A->I*1000.00;

    if (IMPLICIT_FLAG) {
        A->J_O_C = ((S->O_C + dt * (200000. * S->Ca_i)) / (1. + dt * (200000. * S->Ca_i + 476.000)) - S->O_C) / dt;
        A->J_O_TC = ((S->O_TC + dt * (78400.0 * S->Ca_i)) / (1. + dt * (78400.0 * S->Ca_i + 392.000)) - S->O_TC) / dt;
        A->J_O_TMgC = ((S->O_TMgC + dt * (200000. * S->Ca_i * (1.0 - S->O_TMgMg))) / (1. + dt * (200000. * S->Ca_i + 6.60000)) - S->O_TMgC) / dt;
    } else {
        A->J_O_C = 200000. * S->Ca_i * (1.00000 - S->O_C) - 476.000 * S->O_C;
        A->J_O_TC = 78400.0 * S->Ca_i * (1.00000 - S->O_TC) - 392.000 * S->O_TC;
        A->J_O_TMgC = 200000. * S->Ca_i * ((1.00000 - S->O_TMgC) - S->O_TMgMg) - 6.60000 * S->O_TMgC;
    }
    R->O_C = A->J_O_C;
    R->O_TC = A->J_O_TC;
    R->O_TMgC = A->J_O_TMgC;
    A->J_O = 0.08 * A->J_O_TC + 0.16 * A->J_O_TMgC + 0.045 * A->J_O_C;
	R->O = A->J_O;

	A->i_up = ( C->I_up_max*(S->Ca_i/C->k_cyca - ( C->k_xcs*C->k_xcs*S->Ca_up)/C->k_srca))/((S->Ca_i+C->k_cyca)/C->k_cyca+( C->k_xcs*(S->Ca_up+C->k_srca))/C->k_srca);
	A->i_rel_f2 = S->F2/(S->F2+0.250000);
	A->i_rel_factor =  A->i_rel_f2*A->i_rel_f2;
	A->i_rel =  C->alpha_rel*A->i_rel_factor*(S->Ca_rel - S->Ca_i);
	R->Ca_i = - ((A->i_B_Ca+A->i_CaP+A->i_up) - (A->i_di+A->i_rel+ 2.00000*A->i_NaCa))/( 2.00000*C->Vol_i*C->F) -  1.00000*A->J_O;
	A->i_tr = ( (S->Ca_up - S->Ca_rel)*2.00000*C->Vol_rel*C->F)/C->tau_tr;
	R->Ca_up = (A->i_up - A->i_tr)/( 2.00000*C->Vol_up*C->F);

    if (IMPLICIT_FLAG) {
        A->J_O_Calse = ((S->O_Calse + dt * (480.000 * S->Ca_rel)) / (1. + dt * (480.000 * S->Ca_rel + 400.000)) - S->O_Calse) / dt;
    } else {
        A->J_O_Calse = 480.000 * S->Ca_rel * (1.00000 - S->O_Calse) - 400.000 * S->O_Calse;
    }
	R->O_Calse = A->J_O_Calse;

	R->Ca_rel = (A->i_tr - A->i_rel)/( 2.00000*C->Vol_rel*C->F) -  31.0000*A->J_O_Calse;
}

void euler(double dt, State *S, State *R) {
	S->V += dt * R->V;
	S->Na_c += dt * R->Na_c;
	S->Na_i += dt * R->Na_i;
	S->m += dt * R->m;
	S->h1 += dt * R->h1;
	S->h2 += dt * R->h2;
	S->Ca_d += dt * R->Ca_d;
	S->d_L += dt * R->d_L;
	S->f_L1 += dt * R->f_L1;
	S->f_L2 += dt * R->f_L2;
	S->K_c += dt * R->K_c;
	S->K_i += dt * R->K_i;
	S->r += dt * R->r;
	S->s += dt * R->s;
	S->a_ur += dt * R->a_ur;
	S->i_ur += dt * R->i_ur;
	S->n += dt * R->n;
	S->pa += dt * R->pa;
	S->Ca_c += dt * R->Ca_c;
	S->Ca_i += dt * R->Ca_i;
	S->O_C += dt * R->O_C;
	S->O_TC += dt * R->O_TC;
	S->O_TMgC += dt * R->O_TMgC;
	S->O_TMgMg += dt * R->O_TMgMg;
	S->O += dt * R->O;
	S->Ca_rel += dt * R->Ca_rel;
	S->Ca_up += dt * R->Ca_up;
	S->O_Calse += dt * R->O_Calse;
	S->F1 += dt * R->F1;
	S->F2 += dt * R->F2;

    S->d_ord  += dt * R->d_ord;
    S->ff += dt * R->ff;
    S->fs += dt * R->fs;
    S->fcaf += dt * R->fcaf;
    S->fcas += dt * R->fcas;
    S->jca += dt * R->jca;
    S->ffp += dt * R->ffp;
    S->fcafp += dt * R->fcafp;
    S->nca += dt * R->nca;
}


int action_potential(struct State *initial_state, double *scaling_coefficients, double *AP, float CL, float amp,
                     int AP_length, int iso, int baseline_index, int amount_of_baselines, int amount_of_genes) {

    const int chain_length = 1; //?? maybe it should not be hardcoded 
    const int target_cell_index = chain_length / 2;

    const int number_of_stimuli = 1; //?? maybe it should not be hardcoded 
    const double ft = number_of_stimuli * CL; // ms

    const double dt = 1e-2; // ms //?? maybe it should not be hardcoded
    const double baseline_step = 0.1; //?? maybe it should not be hardcoded
     
    const int skip = baseline_step / dt; // number of timesetps to skip in sampling of data in output file

    const int genes_without_concentrations = amount_of_genes - 3 * amount_of_baselines;

    State *S = new State[chain_length];
    State *R = new State[chain_length];
    Algebraic *A = new Algebraic[chain_length];
    Constants *C = new Constants;

    initialize_constants_default(C, CL / 1000., scaling_coefficients);
    for (int i = 0; i < chain_length; ++i) {
        initialize_states(S + i, initial_state);

        S[i].Na_i = scaling_coefficients[genes_without_concentrations + baseline_index];
        S[i].Ca_rel = scaling_coefficients[genes_without_concentrations + amount_of_baselines + baseline_index];
        S[i].Ca_up = scaling_coefficients[genes_without_concentrations + amount_of_baselines + baseline_index];
        S[i].K_i = scaling_coefficients[genes_without_concentrations + 2 * amount_of_baselines + baseline_index];
    }

    double t = 0; // ms
    int AP_length_current = 0, dt_counter = 0;

    while (t <= ft) {

        if (t >= ft - CL && t < ft - CL + AP_length) {
            if (dt_counter % skip == 0) {
                if (AP_length_current < AP_length) {
                    AP[AP_length_current] = S[target_cell_index].V; //?? maybe rewrite because of lots of LL cache write misses
                }
                AP_length_current += 1;
            }
            dt_counter++;
        }

        double VOI = t / 1000.;

        for (int i = 0; i < chain_length; ++i) {
            compute_rates(VOI, S + i, R + i, A + i, C, dt / 1000., (i) ? 0 : amp);

            double g_gap_junc = 1.0;
            double I_gap_junc = 0;
            if (i < chain_length - 1) {
                I_gap_junc += -g_gap_junc * (S[i + 1].V - S[i].V);
            }
            if (i > 0) {
                I_gap_junc += -g_gap_junc * (S[i - 1].V - S[i].V);
            }
            R[i].V -= I_gap_junc * 1000.00;
        }

        for (int i = 0; i < chain_length; ++i) {
            euler(dt / 1000., S + i, R + i);
        }
        t += dt;
    }

    //Rewriting states!
    scaling_coefficients[genes_without_concentrations + baseline_index] = S[target_cell_index].Na_i;
    scaling_coefficients[genes_without_concentrations + amount_of_baselines + baseline_index] = S[target_cell_index].Ca_rel;
    scaling_coefficients[genes_without_concentrations + 2 * amount_of_baselines + baseline_index] = S[target_cell_index].K_i;
    initialize_states(initial_state, S + target_cell_index);

    delete [] S;
    delete [] R;
    delete [] A;
    delete C;

    return 0;
}


