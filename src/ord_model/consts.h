#pragma once
#define STRUCT_SIZE 41	//We need this to update initial states. 41 for O'Hara-Rudy
typedef double REAL;

struct State {
    double v;
    double nai;
    double nass;
    double ki;
    double kss;
    double cai;
    double cass;
    double cansr;
    double cajsr;
    double m;
    double hf;
    double hs;
    double j;
    double hsp;
    double jp;
    double mL;
    double hL;
    double hLp;
    double a;
    double iF;
    double iS;
    double ap;
    double iFp;
    double iSp;
    double d;
    double ff;
    double fs;
    double fcaf;
    double fcas;
    double jca;
    double nca;
    double ffp;
    double fcafp;
    double xrf;
    double xrs;
    double xs1;
    double xs2;
    double xk1;
    double Jrelnp;
    double Jrelp;
    double CaMKt;
};

//constants
extern double const nao;//extracellular sodium in mM
extern double const cao;//extracellular calcium in mM
extern double const ko;//extracellular potassium in mM

//buffer paramaters
extern double const BSRmax;
extern double const KmBSR;
extern double const BSLmax;
extern double const KmBSL;
extern double const cmdnmax;
extern double const kmcmdn;
extern double const trpnmax;
extern double const kmtrpn;
extern double const csqnmax;
extern double const kmcsqn;

//CaMK paramaters
extern double const aCaMK;
extern double const bCaMK;
extern double const CaMKo;
extern double const KmCaM;
extern double const KmCaMK;

//physical constants
extern double const R;
extern double const T;
extern double const F;

//cell geometry
extern double const L;
extern double const rad;
extern double const vcell;
extern double const Ageo;
extern double const Acap;
extern double const vmyo;
extern double const vmito;
extern double const vsr;
extern double const vnsr;
extern double const vjsr;
extern double const vss;

//introduce varaibles for reversal potentials, currents, fluxes, and CaMK
extern double ENa,EK,EKs;
extern double INa,INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1,INaCa_i,INaCa_ss,INaCa,INaK,IKb,INab,IpCa,ICab,Ist;
extern double Jrel,Jup,Jtr,Jdiff,JdiffNa,JdiffK,Jleak;
extern double CaMKa,CaMKb;

//introduce APD, timing, and counting parameters
extern int APD_flag;
extern double APD;
extern double t_vdot_max;
extern double vrest;
extern double vo;
extern double dt;
extern double t0;
extern double t;
extern double dto;
extern double vdot_old;
extern double vdot;
extern double vdot_max;
extern int p;
extern int n_stim;
extern int Count;


extern const double amp;
extern const double start;
extern const double duration;

extern double g_gap_junc;

