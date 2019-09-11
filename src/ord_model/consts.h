#ifndef H_CONSTS
#define H_CONSTS
#define STRUCT_SIZE 40

// Physical  Constants
static const double pi         = 3.14159265358979;
static const double R=8314.;      		// Universal Gas Constant (J/kmol*K)
static const double FRD=96485.; 		// Faraday's Constant (C/mol)
static const double TEMP=310.;   		// Temperature (K)
static const double RTF=26.71;			//R*TEMP/FRD;	// const* 26.7123387

// Capacitances are used by rhs as well as diffusion functions. Global for now.
static const double capacit_at=1.1e-4; 	//uF
static const double capacit_san=5.7e-5;	//uF
static const double capacit_fb=6.3e-6;	//uF


struct State
{
    double E;                                     // Potential
    double m, h, j;                                // Activation and Inactivation Components for Fast Na Current
//    double mp, hp, jp;							//phosph.
    double xkr, xks;                              // Components of Rapidly Activating K Current
//    double xkrp, xks1p, xks2p;							//for phosphorylated channels
    double xtos, ytos, xtof, ytof;                // Slow and Fast Components for Transient Outward K Current
    double d, f; 
//dp, fp;                                  //
    double fCaBJunc, fCaBSl; 
//fCaBJuncp, fCaBSlp;                      //
    double Csqnb;                                 // SR Calcium Buffer
    double RyRr, RyRo, RyRi;                      //
    double NaJ, NaSl, Nai, CaJ, CaSl, Cai, CaSR;  // Na and Ca Concentrations
    double NaBj, NaBSl;                           // Na Buffers Concentrations
    double a_kur,i_kur;


//mssp, hssp, jssp;//ina
   
    double TnCl, TnChc, TnChm, CaM, Myoc, Myom, SRB;                              // Cytosolic Ca Buffers
    double SLLj, SLLsl, SLHj, SLHsl;                                              // Junctional and SL Ca Buffers
    

    
};


struct Intracellat{ double ENaJunc, ENaSl, EK, ECl; };				//Nernst potentials. Don't need to keep ECa in memory since it is used only in IbgCa






struct Isat
{
    double INaKJunc, INaKSl;                                                      // Na/K Pump Currents in the Junctional Cleft and the Subsarcolemmal Space
    double INaJunc, INaSl;                                                        // Fast Na Currents
    double INaBkJunc, INaBkSl;                                                    // Background Na Currents
    double ICaBkJunc, ICaBkSl;                                                    // Background Ca Currents
    double ICaJunc, ICaSl, ICaJuncp, ICaSlp;                                                        // L-type Calcium Currents
    double ICaK, ICaNaJunc, ICaNaSl, ICaP, ICaKp,ICaNaJuncp,ICaNaSlp;                                              // Ca-K, Ca-Na Currents
    double IpCaJunc, IpCaSl;                                                      // Sarcolemmal Ca Pump Currents
    double IncxJunc, IncxSl;                                                      // Na/Ca Exchanger Currents
    double INaTotJunc, INaTotSl;                                                  // Total Sodium Currents
    double ICaTotJunc, ICaTotSl;                                                  // Total Calcium Currents
    double ItoSlow, ItoFast;
    double Ikr, Iks;
    double Iki, Ikp;
    double IClCa, IBgCl;
    double Iks_ph, Iks_tot, INaJuncp, INaSlp, INaKJuncp, INaKSlp;					
    double I_kur;

};

struct Jsat
{
    double JSRCaRel, JserCa, JserCaP, JSRLeak;                                            
    double JCaBCytosol, JCaBJunc, JCaBSl;                                        
};


struct Const_at
{ 
//Fractional Currents
 double Fjunc,Fsl,FjuncCaL,FslCaL, Kmr;
 double Cmem;
 double cellLength;
 double cellRadius;
 double juncLength;
 double juncRadius;
 double distSLcyto;
 double distJuncSL;
 double Vcell;
 double Vmyo;
 double Vsr;
 double Vsl;
 double Vjunc;
 double JCaJuncSl;
  double JCaSlMyo;
  double JNaJuncSl;
  double JNaSlMyo;
  double Cli;
  double Mgi;
  double Ki;
  double Clo;
  double Ko;
  double Nao;
  double Cao;
  double gNa;
  double gNaB;
  double iNaK ;
  double KmNaip;
  double KmKo;
  double pNa;
  double pCa;
  double pK;
  double iNCX ;
  double iPMCA ;
  double KmCai;
  double KmCao;
  double KmNai;
  double KmNao;
  double KmPCa;
  double ksat;
  double nu ;
  double KdAct;
  double gCaB ;
  double gClCa;
  double gClB ;
  double KdClCa ;
  double pNaK  ;
  double gkp   ;
  double gki	;
  double gksJunc  ;
  double gksSl    ;
  double gto;
  double gkur;
  double gkr;
  double VmaxSRCaP;
  double Kmf ;
  double mr   ;
  double hillSRCaP;
  double ks   ;
  double koCa_ctrl ;
  double koCa_ISO     ;
  double koM_ctrl     ;
  double koM_ISO;
  double kiCa  ;
  double kiM      ;
  double ec50SR   ;
  double BmaxNaJ  ;
  double BmaxNaSl ;
  double KoffNa    ;
  double KonNa    ;
  double BmaxTnClow ;
  double KoffTnCl;
  double KonTnCl ;
  double BmaxTnChigh;
  double KoffTnChCa;
  double KonTnChCa;
  double KoffTnChMg ;
  double KonTnChMg  ;
  double BmaxCaM   ;
  double KoffCaM   ;
  double KonCaM   ;
  double BmaxMyosin ;
  double KoffMyoCa ;
  double KonMyoCa ;
  double KoffMyoMg ;
  double KonMyoMg ;
  double BmaxSR   ;
  double KoffSR   ;
  double KonSR    ;
  double BmaxSLlowS ;
  double BmaxSLlowJ;
  double KoffSll  ;
  double KonSll   ;
  double BmaxSLhighS;
  double BmaxSLhighJ;
  double KoffSlh   ;
  double KonSlh   ;
  double BmaxCsqn  ;
  double KoffCsqn  ;
  double KonCsqn ;
};

struct san_consts 
{

  double Cmem;
  double L_cell;
  double L_sub;
  double R_cell;
  double V_jsr_part;
  double V_i_part;
  double V_nsr_part;
  double V_cell;
  double V_sub;
  double V_i;
  double V_jsr;
  double V_nsr;
  double Cao;
  double Ki;
  double Ko;
  double Nao;
  double Mgi;
  double g_f_Na;
  double g_f_K;
  double P_CaL;
  double P_CaT;
  double g_Kr;
  double g_ks;
  double g_KACh;
  double g_to;
  double g_Na;
  double I_NaKmax;
  double K_NaCa;
  double g_Kur;
  double Km_fCa;
  double Km_Kp;
  double Km_Nap;
  double alfa_fCa;
  double K1_ni;
  double K1_no;
  double K2_ni;
  double K2_no;
  double K3_ni;
  double K3_no;
  double K_ci;
  double K_cni;
  double K_co;
  double Q_ci;
  double Q_co;
  double Q_n;
  double  tau_dif_Ca;
  double tau_tr;
  double K_up;
  double P_up;
  double slope_up;
  double kiCa;
  double kim;
  double koCa;
  double kom;
  double ks;
  double EC50_SR;
  double HSR;
  double MaxSR;
  double MinSR;
  double CM_tot;
  double CQ_tot;
  double TC_tot;
  double TMC_tot;
  double kb_CM;
  double kb_CQ;
  double kb_TC;
  double kb_TMC;
  double kkb_TMM;
  double kf_CM;
  double kf_CQ;
  double kf_TC;
  double kf_TMC;
  double kf_TMM;
} ;


struct Const_fb {double Vi, cmfb, gbnafb, I_nakfb, kmko, kmnai, gk1, gkv, gto, naofb, kofb;};// *Cfb;


#endif					//ifndef H_CONSTS
