#define READIS 0
#define WRITEIS 0
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "consts.h"

#define DURSTIM 1




//------------------Output Files--------------------------------
/*FILE *E, *txr;                                     // Potential
FILE *m, *h, *j;                               // Activation and Inactivation Components for Fast Na Current
FILE *xkr, *xks;                              // Components of Rapidly Activating K Current
FILE *xtos, *ytos, *xtof, *ytof;                // Slow and Fast Components for Transient Outward K Current
FILE *d, *f, *dp, *fp;                                  //
FILE *fCaBJunc, *fCaBSl;                      //
FILE *fCaBJuncp, *fCaBSlp;                      //
FILE *Csqnb;                                 // SR Calcium Buffer
FILE *RyRr, *RyRo, *RyRi;                      //
FILE *NaJ, *NaSl, *Nai, *CaJ, *CaSl, *Cai, *CaSR;  // Na and Ca Concentrations
FILE *NaBj, *NaBSl;
FILE *I_kur;

FILE *ItoSlow, *ItoFast, *Iks, *Iki, *INaKJunc, *INaKJuncPlusINaKSl, *INaKSl, *IncxJunc, *IncxSl,*IncxJuncPlusIncxSl;
FILE *INaBkJunc, *INaBkSl, *ICaBkJunc, *ICaBkSl, *ICaJunc, *ICaSl, *ICaK, *ICaNaJunc, *ICaNaSl,
        *INaTotJunc, *INaTotSl,*ICaTotJunc, *ICaTotSl, *Ikr, *Ikp, *IClCa, *IBgCl;

//FILE *inafile, *ICaNa, *INaK, *CaJt, *jsrCaRel;

FILE *Iepi, *Iendo, *Ikl, *icalfile, *Iks_tot, *xkrp, *xksp, *Iksp;

FILE *gkrf, *rkrf, *tauxkr, *xkrf, *xkr1f, *xrssf;


//for iki
FILE *aki, *bki,*kiss;

//for ical
FILE *Dss, *Fss,*Td,*Tf, *icaj, *icas, *inaj, *inas, *ik, *Ca_i, *Caitmax;*/
//  bool el_only,ficalp;
double Cai_t_max, casr;






// Model Parameters

void set_at_consts(struct Const_at *Cat, double* scaling_coefficients)
{
//Fractional Currents
  Cat->Fjunc      = 0.11;
  Cat->Fsl =0.89;// 1.0 - Fjunc;
  Cat->FjuncCaL   = 0.9;
  Cat->FslCaL     =0.1;// 1.0-FjuncCaL;
  Cat->Cmem      = 1.1e-10;
// Enviromental Parameters
  Cat->cellLength = 100.0;                 // (um)
  Cat->cellRadius = 10.25;               // (um)
  Cat->juncLength = 160e-3;              // (um)
  Cat->juncRadius = 15e-3;               // (um)
  Cat->distSLcyto = 0.45;                // (um)
  Cat->distJuncSL = 0.5;                 // (um)
  //Cat->DcaJuncSL  = 1.64e-6;             // (cm^2/sec)
  //Cat->DcaSLcyto  = 1.22e-6;             // (cm^2/sec)
  //Cat->DnaJuncSL  = 1.09e-5;             // (cm^2/sec)
  //Cat->DnaSLcyto  = 1.79e-5;             // (cm^2/sec)
  Cat->Vcell      = 33e-12;//(pi*cellRadius*cellRadius*cellLength*1e-15);//33e-12;//      // (L)
  Cat->Vmyo       = 21.45e-12;//(0.65*Vcell);          // (um^2)
  Cat->Vsr        = 1.155e-12;//(0.035*Vcell);         //(0.023*Vcell); // (um^2)
  Cat->Vsl        = 0.66e-12;//(0.02*Vcell);          //(0.013*Vcell);     // (um^2)
  Cat->Vjunc      = 0.178e-12;//0.0539*0.1*Vcell; //(0.0004*Vcell);  //(5.39e-4*Vcell);       // (um^2)
//  Cat->SAjunc     = 303.7;//(40.3e3*pi*juncLength*juncRadius);                // (um^2)
//  Cat->SAsl       = 6437;//(pi*2.0*cellLength*cellRadius);                   // (um^2)
  Cat->JCaJuncSl  = 8.2413e-13;          // (L/ms)
  Cat->JCaSlMyo   = 3.2743e-12;          // (L/ms)
  Cat->JNaJuncSl  =6.1043e-13;// 1.83128e-14;         // (L/ms)
  Cat->JNaSlMyo   =5.4621e-11;// 1.63863e-12;         // (L/ms)

// Ion Concentrations
  Cat->Cli        = 15.0;                  // (mM)
  Cat->Mgi        = 1.0;                   // (mM)
  Cat->Ki         = 120.0;                 // (mM)
  Cat->Clo        = 150.0;                 // (mM)
  Cat->Ko         = 5.4;                 // (mM)
  Cat->Nao        = 140.0;                 // (mM)
  Cat->Cao        = 1.8;                 // (mM)

// Na Transport
  Cat->gNa        = 23.0*scaling_coefficients[3];                  // Conductivity (mS/uF)
  Cat->gNaB       = 0.597e-3;            // Conductivity (mS/uF)
  Cat->iNaK       = 1.26*scaling_coefficients[7];//1.8;                 // Current (A/F)
  Cat->KmNaip     = 11.0;                  // Michaelis  ant (mM)
  Cat->KmKo       = 1.5;                 // Michaelis  ant (mM)

// Ca transport
  Cat->pNa        = (0.5*1.5e-8)*scaling_coefficients[5];          // (cm/sec)
  Cat->pCa        = (0.5*5.4e-4)*scaling_coefficients[5];          // (cm/sec)
  Cat->pK         = (0.5*2.7e-7)*scaling_coefficients[5];          // (cm/sec)
  Cat->iNCX       =3.15*scaling_coefficients[6];// 4.5;                 // Current (A/F)
  Cat->iPMCA      = 0.0673*scaling_coefficients[8];              // Current (A/F)
  Cat->KmCai      = 3.59e-3;             // Michaelis  ant (mM)
  Cat->KmCao      = 1.3;                 // Michaelis  ant (mM)
  Cat->KmNai      = 12.29;               // Michaelis  ant (mM)
  Cat->KmNao      = 87.5;                // Michaelis  ant (mM)
  Cat->KmPCa      = 0.5e-3;              // Michaelis  ant (mM)
  Cat->ksat       = 0.27;//0.32;
  Cat->nu         =0.32;// 0.27;
  Cat->KdAct      = 0.384e-3;//1.5*0.15e-3;             // Dissociation  ant (mM)
  Cat->gCaB       = 6.0643e-4;//5.513e-4;            // (A/F)



// Cl Currents
  Cat->gClCa      = (0.5*0.109625);   // Conductivity (mS/uF)
  Cat->gClB       = 9e-3*scaling_coefficients[12];                 // Conductivity (mS/uF)
  Cat->KdClCa     = 0.1;                 // Dissociation  ant(mM)

// K Currents
  Cat->pNaK       = 0.01833;
  Cat->gkp        = 0.002;               // Conductivity (mS/uF)
  Cat->gksJunc    =0.0035*scaling_coefficients[2];              // Conductivity (mS/uF)
  Cat->gksSl      =0.0035*scaling_coefficients[2];              // Conductivity (mS/uF)

  Cat->gto= 0.165*scaling_coefficients[4];              // Conductivity (mS/uF)
  Cat->gkur=0.045;			// nS/pF
  Cat->gkr=0.035*scaling_coefficients[1];			// nS/pF
  Cat->gki=0.0525*scaling_coefficients[0];			//nS/pF

// SR Ca Fluxes
  Cat->VmaxSRCaP  = 5.3114e-3*scaling_coefficients[9];           // (mM/ms)
  Cat->Kmf        = 2.5*0.246e-3;            // Michaelis  ant (mM)
  Cat->Kmr        = 1.7;                 // Michaelis  ant (mM)
  Cat->hillSRCaP  = 1.787;
  Cat->ks         = 25.0*scaling_coefficients[11];                  // (1/ms)
  Cat->koCa_ctrl       = 10.0;             // (mM^2/ms)
  Cat->koCa_ISO        = 15.0;              // (mM^2/ms)
  Cat->koM_ctrl        = 0.06;             // (1/ms)
  Cat->koM_ISO         = 0.09;              // (1/ms)
  Cat->kiCa       = 0.5;                 // (1/(ms*mM))
  Cat->kiM        = 0.005;               // (1/ms)
  Cat->ec50SR     = 0.45;                // (mM)

// Buffering
  Cat->BmaxNaJ    = 7.561;               // (mM)
  Cat->BmaxNaSl   = 1.65;                // (mM)
  Cat->KoffNa     = 1e-3;                // (1/ms)
  Cat->KonNa      = 1e-4;                // (1/(ms*mM))
  Cat->BmaxTnClow = 7e-2;                // (mM)
  Cat->KoffTnCl   = 0.0196;              // (1/ms)
  Cat->KonTnCl    = 32.7;                // (1/(ms*mM))
  Cat->BmaxTnChigh= 0.14;                // (mM)
  Cat->KoffTnChCa = 0.32e-4;             // (1/ms)
  Cat->KonTnChCa  = 2.37;                // (1/(ms*mM))
  Cat->KoffTnChMg = 3.33e-3;             // (1/ms)
  Cat->KonTnChMg  = 3e-3;                // (1/(ms*mM))
  Cat->BmaxCaM    = 0.024*scaling_coefficients[10];               // (mM)
  Cat->KoffCaM    = 0.238;               // (1/ms)
  Cat->KonCaM     = 34.0;                  // (1/(ms*mM))
  Cat->BmaxMyosin = 0.14;                // (mM)
  Cat->KoffMyoCa  = 0.46e-3;             // (1/ms)
  Cat->KonMyoCa   = 13.8;                // (1/(ms*mM))
  Cat->KoffMyoMg  = 0.57e-4;             // (1/ms)
  Cat->KonMyoMg   = 0.0157;              // (1/(ms*mM))
  Cat->BmaxSR     = (19.0*0.9e-3);         // (mM)
  Cat->KoffSR     = 6e-2;                // (1/ms)
  Cat->KonSR      = 100.0;                 // (1/(ms*mM))
  Cat->BmaxSLlowS = 1.2155;//((Vmyo*1e13)/(Vsl*1e13)*0.0374);   // (mM)
  Cat->BmaxSLlowJ = 0.0554;//(4.6e-4*(Vmyo*1e13)/(Vjunc*1e13));   // (mM)
  Cat->KoffSll    = 1.3;                 // (1/ms)
  Cat->KonSll     = 100.0;                 // (1/(ms*mM))
  Cat->BmaxSLhighS= 0.4355;//(13.4e-3*(Vmyo*1e13)/(Vsl*1e13));    // (mM)
  Cat->BmaxSLhighJ= 0.0199;//(1.65e-4*(Vmyo*1e13)/(Vjunc*1e13));  // (mM)
  Cat->KoffSlh    = 0.03;                // (1/ms)
  Cat->KonSlh     = 100.0;                 // (1/(ms*mM))
  Cat->BmaxCsqn   = 2.6;//(0.14*(Vmyo*1e13)/(Vsr*1e13));       // (mM)
  Cat->KoffCsqn   = 65.0;                  // (1/ms)
  Cat->KonCsqn    = 100.0;
}

int cabuf(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Jsat *J, struct Const_at *C)
{
   // Cytosolic Ca Buffers Time Derivative
   D->TnCl  = C->KonTnCl   * S->Cai * (C->BmaxTnClow  - S->TnCl) - C->KoffTnCl * S->TnCl;
   D->TnChc = C->KonTnChCa * S->Cai * (C->BmaxTnChigh - S->TnChc - S->TnChm) - C->KoffTnChCa * S->TnChc;
   D->TnChm = C->KonTnChMg * C->Mgi    * (C->BmaxTnChigh - S->TnChc - S->TnChm) - C->KoffTnChMg * S->TnChm;
   D->CaM   = C->KonCaM    * S->Cai * (C->BmaxCaM     - S->CaM)  - C->KoffCaM * S->CaM;
   D->Myoc  = C->KonMyoCa  * S->Cai * (C->BmaxMyosin  - S->Myoc  - S->Myom)  - C->KoffMyoCa * S->Myoc;
   D->Myom  = C->KonMyoMg  * C->Mgi    * (C->BmaxMyosin  - S->Myoc  - S->Myom)  - C->KoffMyoMg * S->Myom;
   D->SRB   = C->KonSR     * S->Cai * (C->BmaxSR      - S->SRB)  - C->KoffSR * S->SRB;
    
   // Cytosolic Ca Buffers
   J->JCaBCytosol = D->TnCl + D->TnChc + D->TnChm + D->CaM + D->Myoc + D->Myom + D->SRB;
    
   // Junctional and SL Ca Buffers Time Derivative
   D->SLLj  = C->KonSll * S->CaJ  * (C->BmaxSLlowJ  - S->SLLj)  - C->KoffSll * S->SLLj;
   D->SLLsl = C->KonSll * S->CaSl * (C->BmaxSLlowS  - S->SLLsl) - C->KoffSll * S->SLLsl;
   D->SLHj  = C->KonSlh * S->CaJ  * (C->BmaxSLhighJ - S->SLHj)  - C->KoffSlh * S->SLHj;
   D->SLHsl = C->KonSlh * S->CaSl * (C->BmaxSLhighS - S->SLHsl) - C->KoffSlh * S->SLHsl;
   // Junctional and SL Ca Buffers
    
   J->JCaBJunc = D->SLLj  + D->SLHj;
   J->JCaBSl   = D->SLLsl + D->SLHsl;
    
    return 0;
} 
double caconcat(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Jsat *J, struct Const_at *C)
{
    // SR Calcium Buffer Calculation
    D->Csqnb = C->KonCsqn * S->CaSR * (C->BmaxCsqn - S->Csqnb) - C->KoffCsqn * S->Csqnb;
   
    // Total Calcium Current Calculation in the Junctional Cleft and the Subsarcolemmal Space
    I->ICaTotJunc = I->ICaJunc + I->ICaBkJunc + I->IpCaJunc - 2.0 * I->IncxJunc;
    I->ICaTotSl    = I->ICaSl + I->ICaBkSl + I->IpCaSl - 2.0 * I->IncxSl;
    
    // Calcium Concentrations Time Derivatives Calculation
    D->CaJ  = (-I->ICaTotJunc * C->Cmem/(2.0 * C->Vjunc * FRD) + C->JCaJuncSl/C->Vjunc * (S->CaSl - S->CaJ) - J->JCaBJunc + J->JSRCaRel * C->Vsr/C->Vjunc + J->JSRLeak * C->Vmyo/C->Vjunc);
    D->CaSl =(-I->ICaTotSl   * C->Cmem/(2.0 * C->Vsl * FRD)   + C->JCaJuncSl/C->Vsl * (S->CaJ - S->CaSl) + C->JCaSlMyo/C->Vsl * (S->Cai - S->CaSl) - J->JCaBSl);
    D->Cai  = (-J->JserCa * C->Vsr/C->Vmyo - J->JCaBCytosol + C->JCaSlMyo/C->Vmyo * (S->CaSl - S->Cai));
    D->CaSR =( J->JserCa - J->JSRLeak * C->Vmyo/C->Vsr -J->JSRCaRel - D->Csqnb);

     
    // Return of Total Calcium Current
    return I->ICaTotJunc + I->ICaTotSl;
}


double naconcat(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Jsat *J, struct Const_at *C)
{
    // Total Na Current Calculation in the junctional cleft and the subsarcolemmal space

    I->INaTotJunc = I->INaJunc + I->INaBkJunc + 3.0 * I->IncxJunc + 3.0 * I->INaKJunc + I->ICaNaJunc;

    I->INaTotSl    = I->INaSl + I->INaBkSl + 3.0 * I->IncxSl + 3.0 * I->INaKSl + I->ICaNaSl;
    
    // Sodium Buffer Calculation
    D->NaBj  = C->KonNa * S->NaJ  * (C->BmaxNaJ  - S->NaBj)  - C->KoffNa * S->NaBj;
    D->NaBSl = C->KonNa * S->NaSl * (C->BmaxNaSl - S->NaBSl) - C->KoffNa * S->NaBSl;
    
   // Sodium Concentrations Calculation
   D->NaJ  = -I->INaTotJunc * C->Cmem/(C->Vjunc * FRD) + C->JNaJuncSl/C->Vjunc * (S->NaSl - S->NaJ) - D->NaBj;
   D->NaSl = -I->INaTotSl   * C->Cmem/(C->Vsl * FRD)   + C->JNaJuncSl/C->Vsl * (S->NaJ - S->NaSl) + C->JNaSlMyo/C->Vsl * (S->Nai - S->NaSl) - D->NaBSl;
   D-> Nai  = C->JNaSlMyo/C->Vmyo * (S->NaSl - S->Nai);
   return I->INaTotJunc + I->INaTotSl;
} /** naconcat **/


int srflux(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Jsat *J, struct Const_at *C)
{
      double MaxSR = 15.0;
      double MinSR = 1.0;
    double  kCaSR;
    
    
    kCaSR  = MaxSR - (MaxSR - MinSR)/(1.0 + pow(C->ec50SR/S->CaSR, 2.5));
    double koSRCa = C->koCa_ctrl/kCaSR;
    double kiSRCa = C->kiCa * kCaSR;
    double RI   = 1.0 - S->RyRr - S->RyRo - S->RyRi;
    D->RyRr =((C->kiM * RI - kiSRCa * S->CaJ * S->RyRr) - (koSRCa * S->CaJ * S->CaJ * S->RyRr - C->koM_ctrl * S->RyRo));
    D->RyRo =((koSRCa * S->CaJ * S->CaJ * S->RyRr - C->koM_ctrl * S->RyRo) - (kiSRCa * S->CaJ * S->RyRo - C->kiM * S->RyRi));
    D->RyRi =((kiSRCa * S->CaJ * S->RyRo - C->kiM * S->RyRi) - (C->koM_ctrl * S->RyRi - koSRCa * S->CaJ * S->CaJ * RI));



    
    J->JSRCaRel = C->ks * S->RyRo * (S->CaSR - S->CaJ);
    J->JserCa   = C->VmaxSRCaP * (pow(S->Cai/C->Kmf, C->hillSRCaP) -  pow(S->CaSR/C->Kmr, C->hillSRCaP))/(1.0 + pow(S->Cai/C->Kmf, C->hillSRCaP) +  pow(S->CaSR/C->Kmr, C->hillSRCaP));
    J->JSRLeak  = 5.348e-6 * (S->CaSR - S->CaJ);
    
    return 0;
} /** srflux **/

//  Background Na Current
double ibgna(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
    I->INaBkJunc = C->Fjunc * C->gNaB * (S->E - Sc->ENaJunc);
    I->INaBkSl   = C->Fsl   * C->gNaB * (S->E - Sc->ENaSl);
    
    return I->INaBkJunc + I->INaBkSl;
} /** ibgna **/


// Backgound Cl Current
double ibgcl(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{

    I->IBgCl = C->gClB * (S->E - Sc->ECl);
    return I->IBgCl;
} /** ibgcl **/


// Backgound Ca Current
double ibgca(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{

    double ECaJunc = RTF/ 2.0 * log(C->Cao/S->CaJ);
    double ECaSl   = RTF/ 2.0 * log(C->Cao/S->CaSl);
    I->ICaBkJunc = C->Fjunc * C->gCaB * (S->E - ECaJunc);
    I->ICaBkSl   = C->Fsl * C->gCaB * (S->E - ECaSl);
    
    return I->ICaBkJunc + I->ICaBkSl;
} /** ibgca **/



double ical(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
    double ICaJ, ICaS, INaJ, INaS, IK;
    double ICa, ICaNa;
    double dss, fss;
    double td, tf; 
    // d calculation
    dss   = 1.0/(1.0 + exp(-(S->E + 9)/6));
    if (S->E == -5.){ td = dss*1/(6*0.035);}
    else {td = dss * (1.0 - exp(-(S->E + 9)/6.0))/(0.035 * (S->E + 9));}
    D->d=(dss-S->d)/td;

    // f calculation
    fss   = 1.0/(1.0 + exp((S->E + 30.0)/7.0)) + 0.2/(1.0 + exp(-(S->E - 50.0)/20.0));
    tf    = 1.0/(0.0197 * exp(-0.0337 * (S->E + 25) * (S->E + 25)) + 0.02);
    D->f=(fss-S->f)/tf;

    D->fCaBJunc = (1.7 * S->CaJ * (1.0 - S->fCaBJunc) - 11.9E-3 * S->fCaBJunc);
    D->fCaBSl   = (1.7 * S->CaSl * (1.0 - S->fCaBSl)   - 11.9E-3 * S->fCaBSl);


    
    // Currents' Amplitude Calculation
    ICaJ = 1.364 * C->pCa * FRD/RTF * S->E * (S->CaJ  * exp(2.0 * (S->E)/RTF) - C->Cao)/(exp(2.0 * (S->E)/RTF) - 1.0);
    ICaS = 1.364 * C->pCa * FRD/RTF * S->E * (S->CaSl * exp(2.0 * (S->E)/RTF) - C->Cao)/(exp(2.0 * (S->E)/RTF) - 1.0);
    INaJ = 0.75 * C->pNa * FRD/RTF * S->E * (S->NaJ  * exp((S->E)/RTF) - C->Nao)/(exp((S->E)/RTF) - 1.0);
    INaS = 0.75 * C->pNa * FRD/RTF * S->E * (S->NaSl * exp((S->E)/RTF) - C->Nao)/(exp((S->E)/RTF) - 1.0);
    IK   =  0.75 * C->pK  * FRD/RTF * S->E * (C->Ki * exp((S->E)/RTF) - C->Ko)/(exp((S->E)/RTF) - 1.0);		//Ki is constant here!
                       
    // Ca, Ca-K and Ca-Na Currents Calculation
    I->ICaJunc   = 0.45 * C->FjuncCaL * ICaJ * S->d * S->f * (1.0 - S->fCaBJunc);
    I->ICaSl     = 0.45 * C->FslCaL   * ICaS * S->d * S->f * (1.0 - S->fCaBSl);
    ICa  =( I->ICaJunc + I->ICaSl);//2
    I->ICaK   = IK * S->d * S->f * (C->FjuncCaL * (1.0 - S->fCaBJunc) + C->FslCaL * (1 - S->fCaBSl)) * 0.45;

    I->ICaNaJunc = C->FjuncCaL * INaJ * S->d * S->f * (1 - S->fCaBJunc) * 0.45;
    I->ICaNaSl = C->FslCaL *INaS * S->d * S->f * (1.0 - S->fCaBSl) * 0.45;

    ICa +=  (I->ICaK + I->ICaNaJunc + I->ICaNaSl);
    
    return ICa;
} /** ical **/


double iclca(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
    double        IClCaJunc, IClCaSl;
    
   
    IClCaJunc = C->Fjunc * C->gClCa * (S->E - Sc->ECl)/(1.0 + C->KdClCa/S->CaJ);
    IClCaSl   = C->Fsl   * C->gClCa * (S->E - Sc->ECl)/(1.0 + C->KdClCa/S->CaSl);
    I->IClCa = IClCaJunc + IClCaSl;
    
    return I->IClCa;
} /** iclca **/





double iki(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{    
    double        aKi, bKi;
    double        Kiss;

    aKi = 1.02 / (1.0 + exp(0.2385 * (S->E - Sc->EK - 59.215)));
    bKi = (0.49214 * exp(0.08032 * (S->E - Sc->EK + 5.476)) + exp(0.06175 * (S->E - Sc->EK - 594.31))) /
              (1.0 + exp(-0.5143 * (S->E - Sc->EK + 4.753)));
    Kiss = aKi / (aKi + bKi);
    I->Iki =C->gki*sqrt(C->Ko / 5.4) * Kiss * (S->E - Sc->EK);

    return I->Iki;
} /** iki **/

double ikp(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
    double        IkpJunc, IkpSl;
    double        kp;
    
    kp      = 1.0/(1.0 + exp(7.488 - S->E/5.98));
    IkpJunc = C->Fjunc * C->gkp * kp * (S->E - Sc->EK);
    IkpSl   = C->Fsl   * C->gkp * kp * (S->E - Sc->EK);

    I->Ikp = IkpJunc + IkpSl;
    
    return I->Ikp;
} /** ikp **/

double ikr_at(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
    double        rkr;


    double xrss =   1.0/(1.0 + exp(-(S->E + 10.0)/5.0));
    double txr  = 550.0/(1.0 + exp(-(S->E + 22.0)/9.0)) * 6.0/(1.0 + exp((S->E + 11.0)/9.0)) + 230.0/(1.0 + exp( (S->E + 40.0)/20.0));

  

    rkr  =   1.0/(1.0 + exp((S->E + 74.0)/24.0));
    D->xkr=(xrss-S->xkr)/txr;
    I->Ikr = C->gkr *sqrt(C->Ko/5.4)* S->xkr * rkr * (S->E - Sc->EK);
    return I->Ikr;

} /** ikr_at **/


double iks_at(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
	double EKsJunc = RTF * log((C->Ko + C->pNaK * C->Nao)/(C->Ki + C->pNaK * S->NaJ));
	double EKsSl   = RTF * log((C->Ko + C->pNaK * C->Nao)/(C->Ki + C->pNaK * S->NaSl));
	double xsss = 1 / (1+exp(-S->E+3.8)/14.25);
	double txs     = 990.1/(1.0 + exp(-(S->E + 2.436)/14.12));
	D->xks= (xsss-S->xks)/txs;
	double I_ks_junc = C->Fjunc*C->gksJunc*S->xks*S->xks*(S->E - EKsJunc);
	double I_ks_sl = C->Fsl*C->gksSl*S->xks*S->xks*(S->E - EKsSl);                                                                                                                                  
	I->Iks_tot = I_ks_junc+I_ks_sl;
return I->Iks_tot;
} /** iks_at **/


 double ikur(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{  
	double a_inf=1/(1+exp((S->E+6)/-8.6));
	double i_inf=1/(1+exp((S->E+7.5)/10));
	double tau_a=(9/(1+exp((S->E+5)/12))+0.5);
	double tau_i=(590/(1+exp((S->E+60)/10))+3050);
	
	D->a_kur =(a_inf-S->a_kur)/tau_a;
	D->i_kur =(i_inf-S->i_kur)/tau_i;

	I->I_kur=C->gkur*S->a_kur*S->i_kur*(S->E-Sc->EK);
	return I->I_kur;

} 


double ina(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
 
 double        ah, bh, aj, bj;

    // m calculation
    double mss   = pow(1.0/(1.0 + exp(-(56.86 + S->E)/9.03)), 2);

    double tm    = 0.1292 * exp(-pow((S->E + 45.79)/15.54, 2)) + 0.06487 * exp(-pow((S->E + 4.823)/51.12, 2));


    // h calculation
    double hss   = pow(1.0/(1.0 + exp((71.55 + S->E)/7.43)), 2);
    if (S->E < -40.0)
    {
        ah = 0.057 * exp(-(S->E + 80)/6.8);
        bh = 2.7   * exp(0.079 * S->E) + 3.1 * 1E5 * exp(0.3485*S->E);
    }
    else
    {
        ah = 0;
        bh = 0.77/(0.13 * (1 + exp(-(S->E + 10.66)/11.1)));
    }
    double th    = 1.0/(ah + bh);

    
    // j calculation
    double jss = pow(1.0/(1.0 + exp((S->E + 71.55)/7.43)), 2);
    if (S->E < -40.0)
    {
        aj = (-2.5428 * 1E4 * exp(0.2444 * S->E) - 6.9481 * 1E-6 * exp(-0.04391 * S->E)*(S->E + 37.78))/(1.0 + exp(0.311 * (S->E + 79.23)));
        bj =  0.02424 * exp(-0.01052 * S->E)/(1.0 + exp(-0.1378 * (S->E + 40.14)));
    }
    else
    {
        aj = 0;
        bj = 0.6 * exp(0.057 * S->E)/(1.0 + exp(-0.1 * (S->E + 32)));
    }
    double tj    = 1.0/(aj + bj);

    //Derivatives
    D->m = (mss - S->m)/tm;//rl
    D->h = (hss - S->h)/th;//rl
    D->j = (jss - S->j)/tj;//rl
    
    // Fast Na Current Calculation in the Junctional Cleft and the Subsarcolemmal Space
    I->INaJunc = C->Fjunc * C->gNa * pow(S->m, 3) * S->h * S->j * (S->E - Sc->ENaJunc);
    I->INaSl   = C->Fsl   * C->gNa * pow(S->m, 3) * S->h * S->j * (S->E - Sc->ENaSl);

    
    return I->INaJunc + I->INaSl;
} /** ina **/


double inak(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
    double sigma, fNaK;
    

    sigma    = (exp(C->Nao/67.3) - 1.0)/7.0;
    fNaK     = 1.0/(1.0 + 0.1245*exp((-0.1*S->E)/RTF) + 0.0365*sigma*exp((-S->E)/RTF));

    I->INaKJunc =(C->Fjunc * C->iNaK * fNaK/(1.0 + pow(C->KmNaip/S->NaJ, 4))  * C->Ko/(C->Ko + C->KmKo));
    I->INaKSl   = (C->Fsl   * C->iNaK * fNaK/(1.0 + pow(C->KmNaip/S->NaSl, 4)) * C->Ko/(C->Ko + C->KmKo));

    return (I->INaKJunc + I->INaKSl);

} /** inak **/


double incx(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
    double KaJunc, KaSl;
    double s1junc, s2junc, s3junc;
    double s1sl,   s2sl,   s3sl;
    
    KaJunc = 1.0/(1.0 + (C->KdAct/S->CaJ)  * (C->KdAct/S->CaJ));
    KaSl   = 1.0/(1.0 + (C->KdAct/S->CaSl) * (C->KdAct/S->CaSl));
    
    s1junc = exp(C->nu * S->E/RTF) * (S->NaJ * S->NaJ * S->NaJ) * C->Cao;
    s1sl   = exp(C->nu * S->E/RTF) * (S->NaSl * S->NaSl * S->NaSl) * C->Cao;
    
    s2junc = exp((C->nu - 1) * S->E/RTF) * (C->Nao*C->Nao*C->Nao) * S->CaJ;
    s2sl   = exp((C->nu - 1) * S->E/RTF) * (C->Nao*C->Nao*C->Nao) * S->CaSl;

   
    s3junc = (C->Nao*C->Nao*C->Nao) * C->KmCai * (1.0 + (S->NaJ/C->KmNai) * (S->NaJ/C->KmNai) * (S->NaJ/C->KmNai)) + (C->KmNao*C->KmNao*C->KmNao) * S->CaJ * (1.0 + S->CaJ/C->KmCai) +
	    (C->KmCao + C->Cao) * (S->NaJ * S->NaJ * S->NaJ) +  S->CaJ * (C->Nao*C->Nao*C->Nao);
    s3sl   = (C->Nao*C->Nao*C->Nao) * C->KmCai * (1.0 + (S->NaSl/C->KmNai) * (S->NaSl/C->KmNai) * (S->NaSl/C->KmNai)) + (C->KmNao*C->KmNao*C->KmNao) * S->CaSl * (1.0 + S->CaSl/C->KmCai) +
	    (C->KmCao + C->Cao) * (S->NaSl * S->NaSl * S->NaSl) +  S->CaSl * (C->Nao*C->Nao*C->Nao);
    

    I->IncxJunc =((C->Fjunc * C->iNCX * KaJunc * (s1junc - s2junc))/(s3junc * (1.0 + C->ksat * exp((C->nu - 1) * (S->E)/RTF)))); 
    I->IncxSl   =((C->Fsl * C->iNCX * KaSl * (s1sl - s2sl))/(s3sl * (1.0 + C->ksat * exp((C->nu - 1) * (S->E)/RTF)))); 
    
    return I->IncxJunc + I->IncxSl;
} /** incx **/


double ipca(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
    I->IpCaJunc = C->Fjunc * C->iPMCA * pow(S->CaJ,  1.6)/(pow(C->KmPCa, 1.6) + pow(S->CaJ,  1.6));
    I->IpCaSl   = C->Fsl   * C->iPMCA * pow(S->CaSl, 1.6)/(pow(C->KmPCa, 1.6) + pow(S->CaSl, 1.6));
    
    return I->IpCaJunc + I->IpCaSl;
} /** ipca **/


double ito(struct State *S, struct State *D,struct Intracellat *Sc, struct Isat *I, struct Const_at *C)
{
   
   
  
double x_inf=1/(1+exp(-(S->E+1)/11.));
double y_inf=1/(1+exp((S->E+40.5)/11.5));
double tau_y=(25.635*exp(-pow((S->E+52.45)/15.8827,2))+24.14);
double tau_x=(3.5*exp(-pow(((S->E)/30),2))+1.5);



D->ytof =(y_inf-S->ytof)/tau_y;
D->xtof =(x_inf-S->xtof)/tau_x;

I->ItoFast=C->gto*S->xtof*S->ytof*(S->E-Sc->EK);





    return I->ItoFast; 

} /** ito **/

// Nernst Potentials Calculation Using Intra- and Extra- Concentrations





void potentialsat(struct State *Stat, struct Intracellat *Scat, struct Const_at *C)
{
    static int first =1;
    Scat->ENaJunc = RTF * log(C->Nao/Stat->NaJ);
    Scat->ENaSl 	= RTF * log(C->Nao/Stat->NaSl);
    for( ; first; first = 0 )
    {
        Scat->EK = RTF * log(C->Ko/C->Ki);	//Ki is constant for now 
        Scat->ECl = RTF * log(C->Cli/C->Clo);
    }

} /** nernstPotentials **/



void Initialize_state_at(double* y,struct State *State)
{
    y[0]          =State->E;
    y[1]          = State->m;
    y[2]          = State->h;
    y[3]        = State->j;
    y[4]        = State->xkr;
    y[5]	= State->xks;
    y[6]	= State->a_kur;
    y[7]	= State->i_kur;    
    y[8]       = State->xtos;
    y[9]       = State->ytos;
    y[10]       = State->xtof;
    y[11]       = State->ytof;
    y[12]        = State->d;
    y[13]        = State->f;
    y[14]   = State->fCaBJunc;
    y[15]     = State->fCaBSl;
    y[16]      = State->Csqnb;
    y[17]       = State->RyRr;
    y[18]       = State->RyRo;
    y[19]       = State->RyRi;
    
    y[20]    = State->NaJ;
    y[21]   = State->NaSl;
    y[22]    =State->Nai;
    y[23]    =State->CaJ;
    y[24]   = State->CaSl;
    y[25]   = State->CaSR;
    y[26]   = State->Cai;
    y[27]   = State->NaBj;
    y[28]  = State->NaBSl;

    y[29]    = State->TnCl;
    y[30]   = State->TnChc;
    y[31]   = State->TnChm;
    y[32]    = State->CaM;
    y[33]    = State->Myoc;
    y[34]    = State->Myom;
    y[35]    = State->SRB;
    y[36]    = State->SLLj;
    y[37]   = State->SLLsl;
    y[38]    = State->SLHj;
    y[39]   = State->SLHsl;
} /** stateInitialisation **/






void Y_to_state(struct State *State, double *y)
{
	State->E=    y[0];		//E
 	State->m=   y[1];          //m
 	State->h=   y[2];          //h
 	State->j=   y[3];        //j
 	State->xkr=   y[4];      //xkr
 	State->xks=   y[5];	//xks	
 	State->a_kur= y[6];			//akur
	State->i_kur=y[7];			//ikur
	State->xtos=   y[8];    //xtos
 	State->ytos=   y[9];    //ytos
 	State->xtof=   y[10];    //xtof
 	State->ytof=   y[11];    //ytof
 	State->d=   y[12];      //d
 	State->f=   y[13];      //f
 	State->fCaBJunc=y[14];	//fCaBJunc
 	State->fCaBSl=   y[15];  //fCaBSL
 	State->Csqnb=   y[16];  //Csqnb
 	State->RyRr=  y[17];    //RyRr
 	State->RyRo=   y[18];   //RyRo
 	State->RyRi=   y[19];   //RyRi
    
 	State->NaJ=   y[20];    //NaJ
 	State->NaSl=   y[21];   //NaSl
 	State->Nai=   y[22];    //Nai
 	State->CaJ=   y[23];    //CaJ
 	State->CaSl=   y[24];   //CaSl
 	State->CaSR=   y[25];   //CaSR
 	State->Cai=   y[26];   //Cai
 	State->NaBj=   y[27];  //NaBj
 	State->NaBSl=   y[28]; //NaBSl

 	State->TnCl=   y[29];  //TnCl
 	State->TnChc=   y[30]; //TnChc
 	State->TnChm=   y[31]; //TnChm
 	State->CaM=   y[32];   //CaM
 	State->Myoc=   y[33];  //Myoc
 	State->Myom=   y[34];  //Myom
 	State->SRB=   y[35];   //SRB
 	State->SLLj=   y[36];  //SLLj
 	State->SLLsl=   y[37]; //SLLsl
 	State->SLHj=   y[38];  //SLHj
 	State->SLHsl=   y[39]; //SLHsl

}



double fitotalat(double* y, double *dy, struct Const_at *C){
  
  
  static struct Jsat Js;
  static struct Isat Is;
  
  static struct State Stat, dState;
  
  static struct Intracellat Scat;
  
  Y_to_state(&Stat, y);


  potentialsat(&Stat, &Scat, C);
 
  double Itot = 0;
  double Ikur=ikur(&Stat, &dState, &Scat, &Is, C);
  double Ical=ical(&Stat, &dState, &Scat, &Is, C);
  double Incx=incx(&Stat, &dState, &Scat, &Is, C);
  double Inak=inak(&Stat, &dState, &Scat, &Is, C);
  double Ito=ito(&Stat, &dState, &Scat, &Is, C);
  double Iki=iki(&Stat, &dState, &Scat, &Is, C);
  double Ikp=ikp(&Stat, &dState, &Scat, &Is, C);
  double Ipca=ipca(&Stat, &dState, &Scat, &Is, C);
  double Ikr=ikr_at(&Stat, &dState, &Scat, &Is, C);
  double Iks=iks_at(&Stat, &dState, &Scat, &Is, C);
  double Iclca=iclca(&Stat, &dState, &Scat, &Is, C);
  double Ibgcl=ibgcl(&Stat, &dState, &Scat, &Is, C);
  double Ibgca=ibgca(&Stat, &dState, &Scat, &Is, C);
  double Ibgna=ibgna(&Stat, &dState, &Scat, &Is, C);
  double Ina=ina(&Stat, &dState, &Scat, &Is, C);
    
  srflux(&Stat, &dState, &Scat, &Is, &Js, C);
  cabuf(&Stat, &dState, &Scat, &Is, &Js,  C);
  double ICatot=caconcat(&Stat, &dState, &Scat, &Is, &Js, C);
  double INatot=naconcat(&Stat, &dState, &Scat, &Is, &Js, C);

  Itot = INatot+ICatot+(Is.ItoSlow + Is.ItoFast) + Is.Ikr + Is.Iki - 2.0 * (Is.INaKJunc + Is.INaKSl) + Is.Ikp +
            Is.ICaK+Is.Iks_tot+Is.IClCa+Is.IBgCl+Is.I_kur;

  /////////////////State-to-array//////////////////////////////////
	

        dy[0] = -Itot; // E	
 
 	dy[1]=dState.m;          //m
 	dy[2]=dState.h;          //h
 	dy[3]=dState.j;        //j
 	dy[4]=dState.xkr;      //xkr
 	dy[5]=dState.xks;	//xks	
 	dy[6]=dState.a_kur;			//akur
	dy[7]=dState.i_kur;			//ikur
	dy[8]=0;    //xtos No slow component in Grandi atrial model.
 	dy[9]=0;    //ytos No slow component in Grandi atrial model.

 	dy[10]=dState.xtof;    //xtof
 	dy[11]=dState.ytof;    //ytof
 	dy[12]=dState.d;      //d
 	dy[13]=dState.f;      //f
 	dy[14]=dState.fCaBJunc;	//fCaBJunc
 	dy[15]=dState.fCaBSl;  //fCaBSL
 	dy[16]=dState.Csqnb;  //Csqnb
 	dy[17]=dState.RyRr;    //RyRr
 	dy[18]=dState.RyRo;   //RyRo
 	dy[19]=dState.RyRi;   //RyRi
    
 	dy[20]=dState.NaJ;    //NaJ
 	dy[21]=dState.NaSl;   //NaSl
 	dy[22]=dState.Nai;    //Nai
 	dy[23]=dState.CaJ;    //CaJ
 	dy[24]=dState.CaSl;   //CaSl
 	dy[25]=dState.CaSR;   //CaSR
 	dy[26]=dState.Cai;   //Cai
 	dy[27]=dState.NaBj;  //NaBj
 	dy[28]=dState.NaBSl; //NaBSl

 	dy[29]=dState.TnCl;  //TnCl
 	dy[30]=dState.TnChc; //TnChc
 	dy[31]=dState.TnChm; //TnChm
 	dy[32]=dState.CaM;   //CaM
 	dy[33]=dState.Myoc;  //Myoc
 	dy[34]=dState.Myom;  //Myom
 	dy[35]=dState.SRB;   //SRB
 	dy[36]=dState.SLLj;  //SLLj
 	dy[37]=dState.SLLsl; //SLLsl
 	dy[38]=dState.SLHj;  //SLHj
 	dy[39]=dState.SLHsl; //SLHsl

   int i;
//   for(i=0;i<=STRUCT_SIZE;i++) dy[i]=dy[i]*1000; // !!!!!!!!!!! ms->s

  return dy[0];


}


void Euler(double* Var, double* dVar, double timeStep, int structsize)
{
	for(int i=0;i<structsize;i++)
	{
		Var[i]+=dVar[i]*timeStep;
	}
}







int atrium(double* y,double* dy,struct Const_at *Cat){
 
    static int first=1;	//for the first call of function
    fitotalat(y, dy, Cat);

    return 0;
}

int action_potential(struct State *initial_state, double *scaling_coefficients, double *AP, float CL, float amp, int current_time, int iso, int baseline_index, int amount_of_baselines, int amount_of_genes)
{
    const double ft = 9 * CL;
	double t=0;
    double timeStep=1e-3;
	int cnt=0, Count=0;

	const int skip = 0.1/timeStep;//200;//number of timesetps to skip in sampling of data in output file

	int chain_length = 1;
    int target_cell = 0;
	int z;

    int genes_without_concentrations = amount_of_genes - 2 * amount_of_baselines;


	struct Const_at Cat;
 	set_at_consts(&Cat, scaling_coefficients);

	double *Var, *dVar;

	if ((Var = (double*)calloc(chain_length*STRUCT_SIZE, sizeof(double))) == NULL) {printf ("not enough memory\n"); exit(-1);}
	if ((dVar = (double*)calloc(chain_length*STRUCT_SIZE, sizeof(double))) == NULL) {printf ("not enough memory\n"); exit(-1);}
//	if ((DiffusionCurrent = (double*)calloc(chain_length, sizeof(double))) == NULL) {printf ("not enough memory\n"); exit(-1);}

    for(z=0;z<chain_length;z++)
    {
		
    	
	Initialize_state_at(&Var[STRUCT_SIZE*z], initial_state);
	
        initial_state->Nai = scaling_coefficients[genes_without_concentrations + baseline_index];
        initial_state->NaSl = scaling_coefficients[genes_without_concentrations + baseline_index];
	initial_state->NaJ = scaling_coefficients[genes_without_concentrations + baseline_index];
        initial_state->CaSR = scaling_coefficients[genes_without_concentrations + amount_of_baselines + baseline_index];
    }


    while (t<=ft)
    {
	if (t>=ft-CL && t<ft-CL+current_time)
            {
                if (Count%skip==0)
                {
                    if (cnt<current_time) AP[cnt] = Var[STRUCT_SIZE*target_cell];
                    cnt+=1;
                }
                Count++;
            }

    	for(z=0;z<chain_length;z++)
    	{
		atrium(&Var[STRUCT_SIZE*z], &dVar[STRUCT_SIZE*z], &Cat);
		Euler(&Var[STRUCT_SIZE*z],&dVar[STRUCT_SIZE*z],timeStep,STRUCT_SIZE);
	}
	//STIMULUS
        if(fmod(t,CL)<DURSTIM) Var[0] -=amp*timeStep;
	t=t+timeStep;
    }

    
    //Rewriting states!
    Y_to_state(initial_state, &Var[STRUCT_SIZE*target_cell]);
    scaling_coefficients[genes_without_concentrations + baseline_index] = initial_state->Nai;
    scaling_coefficients[genes_without_concentrations + amount_of_baselines + baseline_index] = initial_state->CaSR;
    return 0;
}
