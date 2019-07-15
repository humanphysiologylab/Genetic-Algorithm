// Copyright (c) 2011-2015 by Thomas O'Hara, Yoram Rudy,
//                            Washington University in St. Louis.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the names of the copyright holders nor the names of its
// contributors may be used to endorse or promote products derived from 
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//

// C++ Implementation of the O'Hara-Rudy dynamic (ORd) model for the
// undiseased human ventricular action potential and calcium transient
//
// The ORd model is described in the article "Simulation of the Undiseased
// Human Cardiac Ventricular Action Potential: Model Formulation and
// Experimental Validation"
// by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
//
// The article and supplemental materails are freely available in the
// Open Access jounal PLoS Computational Biology
// Link to Article:
// http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061
// 
// Email: tom.ohara@gmail.com / rudy@wustl.edu
// Web: http://rudylab.wustl.edu
// 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "consts.h"

#define CHAIN 0 // 28.06.2017
#define ADAPTIVE_STEP 1 // 06.07.2017

const int skip = 1./dt;//200;//number of timesetps to skip in sampling of data in output file
const double safetime = 25.0;//time from the beginning of each beat during which dt is fixed to small values
const double beatssave = 1;//number of beats to save in the output

#define IKATP 0

const int celltype = 0;  //endo = 0, epi = 1, M = 2

void RGC(struct State *CurrentState, struct State *NextState,  int z, int chain_length, double *scaling_coefficients, float amp, float CL, int iso)
{
    for (z=0; z<chain_length; z++)
    {
        CurrentState[z]=NextState[z];
        
        //---------OLD: REVPORTS()---------//
        ENa=(R*T/F)*log(nao/CurrentState[z].nai);
        EK=(R*T/F)*log(ko/CurrentState[z].ki);
        EKs=(R*T/F)*log((ko + 0.01833 * nao)/(CurrentState[z].ki + 0.01833 * CurrentState[z].nai));
        
        //-------------CURRENTS CALCULATIONS--------------//
        
        /* 1. Calculate INa_fast*/
            CaMKb = CaMKo * scaling_coefficients[12] * (1.0 - CurrentState[z].CaMKt)/(1.0+KmCaM/CurrentState[z].cass);
            CaMKa = CaMKb + CurrentState[z].CaMKt;
            double vffrt = CurrentState[z].v * F * F/(R*T);
            double vfrt = CurrentState[z].v * F/(R*T);
            double kmtrpn=0.0005;
            double hss, hssp;
            double mss;
        
            mss=1.0/(1.0+exp((-(CurrentState[z].v + 39.57))/9.871));
            double tm=1.0/(6.765 * exp((CurrentState[z].v + 11.64)/34.77) + 8.552 * exp(-(CurrentState[z].v + 77.42)/5.955));
            hss=1.0/(1+exp((CurrentState[z].v+82.90)/6.086));
        
            double thf=1.0/(1.432e-5 * exp(-(CurrentState[z].v + 1.196)/6.285)+6.149 * exp((CurrentState[z].v + 0.5096)/20.27));
            double ths=1.0/(0.009794*exp(-(CurrentState[z].v + 17.95)/28.05)+0.3343*exp((CurrentState[z].v + 5.730)/56.66));
            double tj=2.038+1.0/(0.02136*exp(-(CurrentState[z].v +100.6)/8.281)+0.3052*exp((CurrentState[z].v + 0.9941)/38.45));
        
            hssp=1.0/(1+exp((CurrentState[z].v + 89.1)/6.086));
        
            double Ahf=0.99;
            double Ahs=1.0-Ahf;
            double jss=hss;
            double thsp=3.0*ths;
            double tjp=1.46*tj;
        
            NextState[z].m = mss - (mss - CurrentState[z].m)*exp(-dt/tm);
            NextState[z].hf = hss - (hss - CurrentState[z].hf) * exp(-dt/thf);
            NextState[z].hs = hss - (hss - CurrentState[z].hs) * exp(-dt/ths);
            NextState[z].j = jss - ( jss - CurrentState[z].j) * exp(-dt/tj);
            NextState[z].hsp = hssp -(hssp - CurrentState[z].hsp) * exp(-dt/thsp);
            NextState[z].jp = jss - (jss - CurrentState[z].jp) * exp(-dt/tjp);
        
            double h = Ahf * CurrentState[z].hf + Ahs * CurrentState[z].hs;
            double hp = Ahf * CurrentState[z].hf + Ahs * CurrentState[z].hsp;
        
            double GNa=75.;
            GNa*=scaling_coefficients[3];
        
            double fINap=(1.0/(1.0+KmCaMK/CaMKa));
        
            INa = GNa * (CurrentState[z].v-ENa) * CurrentState[z].m * CurrentState[z].m * CurrentState[z].m * ((1.0-fINap) * h * CurrentState[z].j + fINap * hp * CurrentState[z].jp);
        
        /* 2. Calculate INa_late*/
            double mLss=1.0/(1.0+exp((-(CurrentState[z].v + 42.85))/5.264));
            double hLss=1.0/(1.0+exp((CurrentState[z].v + 87.61)/7.488));
            double hLssp=1.0/(1.0+exp((CurrentState[z].v +93.81)/7.488));
            double tmL=tm;
            double thL=200.0;
            double thLp=3.0*thL;
        
            NextState[z].mL = mLss - (mLss - CurrentState[z].mL) * exp(-dt/tmL);
            NextState[z].hL = hLss - (hLss - CurrentState[z].hL) * exp(-dt/thL);
            NextState[z].hLp = hLssp - (hLssp - CurrentState[z].hLp) * exp(-dt/thLp);
        
            double GNaL=0.0075*scaling_coefficients[3];
            if (celltype==1)
            {
                GNaL*=0.6;
            }
            double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
            INaL = GNaL * (CurrentState[z].v - ENa) * CurrentState[z].mL * ((1.0-fINaLp) * CurrentState[z].hL + fINaLp * CurrentState[z].hLp);
        
        /* 3. Calculate Ito*/
            double ass=1.0/(1.0+exp((-(CurrentState[z].v - 14.34))/14.82));
            double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(CurrentState[z].v -18.4099)/29.3814)))+3.5/(1.0+exp((CurrentState[z].v + 100.0)/29.3814)));
            double iss=1.0/(1.0+exp((CurrentState[z].v + 43.94)/5.711));
            double delta_epi;
            if (celltype==1)
            {
                delta_epi=1.0-(0.95/(1.0+exp((CurrentState[z].v + 70.0)/5.0)));
            }
            else
            {
                delta_epi=1.0;
            }
            double tiF=4.562+1/(0.3933*exp((-(CurrentState[z].v +100.0))/100.0)+0.08004*exp((CurrentState[z].v + 50.0)/16.59));
            double tiS=23.62+1/(0.001416*exp((-(CurrentState[z].v + 96.52))/59.05)+1.780e-8*exp((CurrentState[z].v + 114.1)/8.079));
            tiF *= delta_epi;
            tiS *= delta_epi;
        
            double assp=1.0/(1.0+exp((-(CurrentState[z].v -24.34))/14.82));
            double dti_develop=1.354+1.0e-4/(exp((CurrentState[z].v - 167.4)/15.89)+exp(-(CurrentState[z].v - 12.23)/0.2154));
            double dti_recover=1.0-0.5/(1.0+exp((CurrentState[z].v + 70.0)/20.0));
            double tiFp=dti_develop*dti_recover*tiF;
            double tiSp=dti_develop*dti_recover*tiS;
  
        
            double AiF=1.0/(1.0+exp((CurrentState[z].v - 213.6)/151.2));
            double AiS=1.0-AiF;
            double i = AiF * CurrentState[z].iF + AiS * CurrentState[z].iS;
            double ip = AiF * CurrentState[z].iFp + AiS * CurrentState[z].iSp;
        
            NextState[z].a = ass - (ass - CurrentState[z].a) * exp(-dt/ta);
            NextState[z].iF = iss - (iss - CurrentState[z].iF) * exp(-dt/tiF);
            NextState[z].iS = iss - (iss - CurrentState[z].iS) * exp(-dt/tiS);
            NextState[z].ap = assp -(assp - CurrentState[z].ap) * exp(-dt/ta);
            NextState[z].iFp = iss - (iss - CurrentState[z].iFp) * exp(-dt/tiFp);
            NextState[z].iSp = iss - (iss - CurrentState[z].iSp) * exp(-dt/tiSp);
        
            double Gto=0.02 * scaling_coefficients[4];
            if (celltype==1)
            {
                Gto*=4.0;
            }
            if (celltype==2)
            {
                Gto*=4.0;
            }
            double fItop=(1.0/(1.0+KmCaMK/CaMKa));
            Ito = Gto * (CurrentState[z].v - EK) * ((1.0 - fItop) * CurrentState[z].a * i + fItop * CurrentState[z].ap * ip);
        
        /* 4-6. Calculate ICaL, ICaNa, ICaK */
            double dss, fss;
            dss=1.0/(1.0+exp((-(CurrentState[z].v+3.940))/4.230));
            double td=0.6+1.0/(exp(-0.05*(CurrentState[z].v +6.0))+exp(0.09*(CurrentState[z].v +14.0)));
            fss=1.0/(1.0+exp((CurrentState[z].v + 19.58)/3.696));
        
            double tff=7.0+1.0/(0.0045*exp(-(CurrentState[z].v + 20.0)/10.0)+0.0045*exp((CurrentState[z].v + 20.0)/10.0));
            double tfs=1000.0+1.0/(0.000035*exp(-(CurrentState[z].v+5.0)/4.0)+0.000035*exp((CurrentState[z].v + 5.0)/6.0));
            double tfcaf=7.0+1.0/(0.04*exp(-(CurrentState[z].v - 4.0)/7.0)+0.04 * exp((CurrentState[z].v - 4.0)/7.0));
            double tfcas=100.0+1.0/(0.00012 * exp(-CurrentState[z].v /3.0) + 0.00012 * exp(CurrentState[z].v /7.0));
            double Afcaf=0.3+0.6/(1.0+exp((CurrentState[z].v - 10.0)/10.0));
        
            double Aff=0.6;
            double Afs=1.0-Aff;
            double f = Aff * CurrentState[z].ff + Afs * CurrentState[z].fs;
            double fcass=fss;
            double Afcas=1.0-Afcaf;
            double fca = Afcaf * CurrentState[z].fcaf + Afcas * CurrentState[z].fcas;
            double tjca=75.0;
            double tffp=2.5*tff;
            double fp = Aff * CurrentState[z].ffp + Afs * CurrentState[z].fs;
            double tfcafp=2.5*tfcaf;
            double fcap = Afcaf * CurrentState[z].fcafp + Afcas * CurrentState[z].fcas;
            double Kmn=0.002;
            double k2n=1000.0;
            double km2n = CurrentState[z].jca * 1.0;
            double anca=1.0/(k2n/km2n+pow(1.0+Kmn/CurrentState[z].cass,4.0));
            double PhiCaL;
        
            if (CurrentState[z].cass<0.03)  PhiCaL = 4.0 * vffrt * (CurrentState[z].cass * exp(2.0*vfrt) - 0.341 * cao)/(exp(2.0 * vfrt)-1.0);
            else  PhiCaL = 4.0 * vffrt * (0.03 * exp(2.0*vfrt) - 0.341 * cao)/(exp(2.0 * vfrt)-1.0);
        
            double PhiCaNa=1.0*vffrt*(0.75 * CurrentState[z].nass * exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
            double PhiCaK = 1.0 * vffrt * (0.75 * CurrentState[z].kss * exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
            double zca=2.0;
        
            NextState[z].d = dss - (dss - CurrentState[z].d) * exp(-dt/td);
            NextState[z].ff = fss - (fss - CurrentState[z].ff) * exp(-dt/tff);
            NextState[z].fs = fss - (fss - CurrentState[z].fs) * exp(-dt/tfs);
            NextState[z].fcaf = fcass - (fcass - CurrentState[z].fcaf) * exp(-dt/tfcaf);
            NextState[z].fcas = fcass - (fcass - CurrentState[z].fcas) * exp(-dt/tfcas);
            NextState[z].jca = fcass - (fcass - CurrentState[z].jca) * exp(-dt/tjca);
            NextState[z].ffp = fss - (fss - CurrentState[z].ffp) * exp(-dt/tffp);
            NextState[z].fcafp = fcass - (fcass - CurrentState[z].fcafp) * exp(-dt/tfcafp);
            NextState[z].nca = anca * k2n/km2n - (anca * k2n/km2n - CurrentState[z].nca) * exp(-km2n * dt);
        
            double PCa=0.0001 * scaling_coefficients[5];
            if (celltype==1)
            {
                PCa*=1.2;
            }
            if (celltype==2)
            {
                PCa*=2.5;
            }
            double PCap=1.1*PCa;
            double PCaNa=0.00125*PCa;
            double PCaK=3.574e-4*PCa;
            double PCaNap=0.00125*PCap;
            double PCaKp=3.574e-4*PCap;
            double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
        
            ICaL=(1.0-fICaLp) * PCa * PhiCaL * CurrentState[z].d * (f * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fca * CurrentState[z].nca) + fICaLp * PCap * PhiCaL * CurrentState[z].d * (fp * (1.0-CurrentState[z].nca) + CurrentState[z].jca * fcap * CurrentState[z].nca);
            ICaNa=(1.0-fICaLp) * PCaNa * PhiCaNa * CurrentState[z].d * (f * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fca * CurrentState[z].nca)+fICaLp * PCaNap * PhiCaNa * CurrentState[z].d * (fp * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fcap * CurrentState[z].nca);
            ICaK=(1.0-fICaLp) * PCaK * PhiCaK * CurrentState[z].d * (f * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fca * CurrentState[z].nca)+fICaLp * PCaKp * PhiCaK * CurrentState[z].d * (fp * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fcap * CurrentState[z].nca);
        
        
        /* 7. Calculate IKr */
            double GKr = 0.046 * scaling_coefficients[1];
            double xrss=1.0/(1.0+exp((-(CurrentState->v + 8.337))/6.789));
            double txrf=12.98+1.0/(0.3652*exp((CurrentState->v - 31.66)/3.869)+4.123e-5 * exp((-(CurrentState->v -47.78))/20.38));
            double txrs=1.865+1.0/(0.06629*exp((CurrentState->v -34.70)/7.355)+1.128e-5 * exp((-(CurrentState->v - 29.74))/25.94));
            double Axrf=1.0/(1.0+exp((CurrentState->v +54.81)/38.21));
            double Axrs=1.0-Axrf;
            NextState[z].xrf = xrss - (xrss - CurrentState[z].xrf) * exp(-dt/txrf);
            NextState[z].xrs = xrss - (xrss - CurrentState[z].xrs) * exp(-dt/txrs);
            double xr = Axrf * CurrentState[z].xrf + Axrs * CurrentState[z].xrs;
            double rkr=1.0/(1.0+exp((CurrentState->v + 55.0)/75.0))*1.0/(1.0+exp((CurrentState->v - 10.0)/30.0));
        
            if (celltype==1)
            {
                GKr*=1.3;
            }
            if (celltype==2)
            {
                GKr*=0.8;
            }
            IKr=GKr*sqrt(ko/5.4)*xr*rkr*(CurrentState[z].v - EK);
        
        /* 8. Calculate IKs */
            double xs1ss=1.0/(1.0+exp((-(CurrentState[z].v + 11.60))/8.932));
            double txs1=817.3+1.0/(2.326e-4*exp((CurrentState[z].v + 48.28)/17.80)+0.001292*exp((-(CurrentState[z].v + 210.0))/230.0));
            double txs2=1.0/(0.01*exp((CurrentState[z].v - 50.0)/20.0)+0.0193 * exp((-(CurrentState[z].v + 66.54))/31.0));
            double xs2ss=xs1ss;
        
            NextState[z].xs1 = xs1ss - (xs1ss - CurrentState[z].xs1) * exp(-dt/txs1);
            NextState[z].xs2 = xs2ss - (xs2ss - CurrentState[z].xs2) * exp(-dt/txs2);
        
            double KsCa=1.0+0.6/(1.0+pow(3.8e-5/CurrentState[z].cai,1.4));
            double GKs=0.0034 * scaling_coefficients[2];
        
            if (celltype==1)
            {
                GKs*=1.4;
            }
            IKs = GKs * KsCa * CurrentState[z].xs1 * CurrentState[z].xs2 * (CurrentState[z].v-EKs);
        
        /* 9. Calculate IK1 */
            double xk1ss=1.0/(1.0+exp(-(CurrentState[z].v + 2.5538 * ko + 144.59)/(1.5692 * ko + 3.8115)));
            double txk1=122.2/(exp((-(CurrentState[z].v + 127.2))/20.36)+exp((CurrentState[z].v + 236.8)/69.33));
            double rk1=1.0/(1.0+exp((CurrentState[z].v +105.8-2.6*ko)/9.493));
        
            NextState[z].xk1 = xk1ss - (xk1ss - CurrentState[z].xk1) * exp(-dt/txk1);
            double GK1 = 0.1908 * scaling_coefficients[0];
            if (celltype==1)
            {
                GK1*=1.2;
            }
            if (celltype==2)
            {
                GK1*=1.3;
            }
            IK1 = GK1 * sqrt(ko) * rk1 * CurrentState[z].xk1 * (CurrentState[z].v - EK);
        
        /* 10. Calculate INaCa */
            double kna1=15.0;
            double kna2=5.0;
            double kna3=88.12;
            double kasymm=12.5;
            double wna=6.0e4;
            double wca=6.0e4;
            double wnaca=5.0e3;
            double kcaon=1.5e6;
            double kcaoff=5.0e3;
            double qna=0.5224;
            double qca=0.1670;
            double hca = exp((qca * CurrentState[z].v * F)/(R * T));
            double hna = exp((qna * CurrentState[z].v * F)/(R * T));
        
            double h1 = 1 + CurrentState[z].nai / kna3 * (1 + hna);
            double h2=(CurrentState[z].nai * hna)/(kna3*h1);
            double h3=1.0/h1;
            double h4 = 1.0 + CurrentState[z].nai/kna1 * (1+ CurrentState[z].nai/kna2);
            double h5=CurrentState[z].nai * CurrentState[z].nai/(h4*kna1*kna2);
            double h6=1.0/h4;
            double h7=1.0+nao/kna3*(1.0+1.0/hna);
            double h8=nao/(kna3*hna*h7);
            double h9=1.0/h7;
            double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
            double h11=nao*nao/(h10*kna1*kna2);
            double h12=1.0/h10;
        
            double k1=h12*cao*kcaon;
            double k2=kcaoff;
            double k3p=h9*wca;
            double k3pp=h8*wnaca;
            double k3=k3p+k3pp;
            double k4p=h3*wca/hca;
            double k4pp=h2*wnaca;
            double k4=k4p+k4pp;
            double k5=kcaoff;
            double k6 = h6 * CurrentState[z].cai * kcaon;
            double k7=h5*h2*wna;
            double k8=h8*h11*wna;
            double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
            double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
            double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
            double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
            double E1=x1/(x1+x2+x3+x4);
            double E2=x2/(x1+x2+x3+x4);
            double E3=x3/(x1+x2+x3+x4);
            double E4=x4/(x1+x2+x3+x4);
            double KmCaAct=150.0e-6;
            double allo=1.0/(1.0+pow(KmCaAct/CurrentState[z].cai,2.0));
            double zna=1.0;
            double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
            double JncxCa=E2*k2-E1*k1;
            double Gncx = 0.0008 * scaling_coefficients[6];
            if (celltype==1)
            {
                Gncx*=1.1;
            }
            if (celltype==2)
            {
                Gncx*=1.4;
            }
            INaCa_i = 0.8 * Gncx * allo * (zna*JncxNa+zca*JncxCa);
        
            h1 = 1 + CurrentState[z].nass/kna3*(1 + hna);
            h2 = (CurrentState[z].nass * hna)/(kna3 * h1);
            h3 = 1.0/h1;
            h4 = 1.0 + CurrentState[z].nass/kna1*(1 + CurrentState[z].nass/kna2);
            h5 = CurrentState[z].nass * CurrentState[z].nass/(h4*kna1*kna2);
            h6 = 1.0/h4;
            h7 = 1.0+nao/kna3*(1.0+1.0/hna);
            h8 = nao/(kna3*hna*h7);
            h9 = 1.0/h7;
            h10= kasymm+1.0+nao/kna1*(1+nao/kna2);
            h11= nao*nao/(h10*kna1*kna2);
            h12= 1.0/h10;
        
            k1 = h12*cao*kcaon;
            k2 = kcaoff;
            k3p= h9*wca;
            k3pp= h8*wnaca;
            k3 =k3p+k3pp;
            k4p=h3*wca/hca;
            k4pp=h2*wnaca;
            k4=k4p+k4pp;
            k5=kcaoff;
            k6 = h6 * CurrentState[z].cass * kcaon;
            k7=h5*h2*wna;
            k8=h8*h11*wna;
        
            x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
            x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
            x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
            x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
        
            E1=x1/(x1+x2+x3+x4);
            E2=x2/(x1+x2+x3+x4);
            E3=x3/(x1+x2+x3+x4);
            E4=x4/(x1+x2+x3+x4);
            KmCaAct=150.0e-6;
            allo = 1.0/(1.0 + pow(KmCaAct/CurrentState[z].cass,2.0));
            JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
            JncxCa=E2*k2-E1*k1;
            INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);
        
            INaCa=INaCa_i+INaCa_ss;
        
        /* 11. Calculate INaK*/
            double k1p=949.5;
            double k1m=182.4;
            double k2p=687.2;
            double k2m=39.4;
            k3p=1899.0;
            double k3m=79300.0;
            k4p=639.0;
            double k4m=40.0;
            double Knai0=9.073;
            double Knao0=27.78;
            double delta=-0.1550;
            double Knai = Knai0 * exp((delta * CurrentState[z].v * F)/(3.0 * R * T));
            if (iso==1)
            {
                Knai*=0.72;//0.7;
            }
            double Knao = Knao0 * exp(((1.0-delta) * CurrentState[z].v * F)/(3.0 * R * T));
            double Kki = 0.5;
            double Kko = 0.3582;
            double MgADP = 0.05;
            double MgATP = 9.8;
            double Kmgatp = 1.698e-7;
            double H = 1.0e-7;
            double eP = 4.2;
            double Khp = 1.698e-7;
            double Knap = 224.0;
            double Kxkur = 292.0;
            double P = eP/(1.0+H/Khp + CurrentState[z].nai/Knap + CurrentState[z].ki/Kxkur);
        
            double a1 = (k1p*pow(CurrentState[z].nai/Knai,3.0))/(pow(1.0+CurrentState[z].nai/Knai,3.0)+pow(1.0 + CurrentState[z].ki/Kki,2.0)-1.0);
            double b1 = k1m*MgADP;
            double a2 = k2p;
            double b2 = (k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
            double a3 = (k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
            double b3 = (k3m*P*H)/(1.0+MgATP/Kmgatp);
            double a4 = (k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
            double b4 = (k4m * pow(CurrentState[z].ki/Kki,2.0))/(pow(1.0 + CurrentState[z].nai/Knai,3.0)+pow(1.0 + CurrentState[z].ki/Kki,2.0)-1.0);
        
            x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
            x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
            x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
            x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
            E1=x1/(x1+x2+x3+x4);
            E2=x2/(x1+x2+x3+x4);
            E3=x3/(x1+x2+x3+x4);
            E4=x4/(x1+x2+x3+x4);
        
            double zk=1.0;
            double JnakNa=3.0*(E1*a3-E2*b3);
            double JnakK=2.0*(E4*b1-E3*a1);
            double Pnak=30 * scaling_coefficients[7];
            if (celltype==1)
            {
                Pnak*=0.9;
            }
            if (celltype==2)
            {
                Pnak*=0.7;
            }
            INaK = Pnak * (zna * JnakNa + zk * JnakK);
        
        /* 12-14. Calculate INab, ICab, IpCa*/
            double xkb=1.0/(1.0+exp(-(CurrentState[z].v - 14.48)/18.34));
            double GKb = 0.003 * 1;
            if (iso==1)
            {
                GKb*=3;//2.5;
            }
            else GKb*=1.;
            if (celltype==1)
            {
                GKb*=0.6;
            }
            IKb = GKb * xkb * (CurrentState[z].v - EK);
        
            double PNab=3.75e-10;
            INab = 1*PNab * vffrt * (CurrentState[z].nai * exp(vfrt)-nao)/(exp(vfrt)-1.0);
        
            double PCab=2.5e-8;
            ICab= 1*PCab*4.0*vffrt*(CurrentState[z].cai * exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
        
            double GpCa=0.0005;
            GpCa *= scaling_coefficients[8];
            IpCa = GpCa * CurrentState[z].cai/(0.0005 + CurrentState[z].cai);
        
        /* 15. Calculate IKatp*/
        #if IKATP
            double GKatp = 155; // nS/pF
            IKatp = GKatp * f_atp * pow(ko / 5.4, 0.3) * (1. / (40. + 3.5 * exp(0.025 * CurrentState[z].v)) * (CurrentState[z].v - EK));
        #endif
        
        
        //----------OLD: STIMULUS()---------//
        
            if ((t>(start+n_stim*CL) && t<(start+duration+n_stim*CL-dt)))
            {
                if (Ist==0) vrest=CurrentState[z].v;
                Ist=amp;
            }
            else if (t>(start+duration+n_stim*CL-dt))
            {
                Ist = 0.0;
                n_stim = n_stim+1;
            }
            vo = CurrentState[z].v;
        
        //----------OLD: VOLTAGE()---------//
            double I_gap_junc = 0;

            NextState[z].v = CurrentState[z].v-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab);
            #if CHAIN
                if(z<(chain_length-1))	I_gap_junc += -g_gap_junc*(CurrentState[z+1].v-CurrentState[z].v);
                if(z>0)	I_gap_junc += -g_gap_junc*(-CurrentState[z].v + CurrentState[z-1].v);
                NextState[z].v -= dt*I_gap_junc;
                if (z==0) NextState[z].v -= dt*Ist;
            #else
                NextState[z].v -= dt*Ist;
            #endif
        
            #if IKATP
                NextState[z].v -= dt*IKatp;
            #endif
        
        //----------OLD: dVdt_APD()---------//
        
            vdot_old = vdot;
            vdot = (CurrentState[z].v - vo)/dt;
            if (APD_flag == 0 && CurrentState[z].v > -40 && vdot < vdot_old)
            {
                vdot_max = vdot_old;
                t_vdot_max = t - dt;
                APD_flag = 1;
            }
            if	(APD_flag==1 && CurrentState[z].v < 0.9 * vrest)
            {
                APD = t - t_vdot_max;
                APD_flag = 0;
            }
        
        //------------OLD: FBS------------//
        
            double CaMKb = CaMKo * scaling_coefficients[12]*(1.0 - CurrentState[z].CaMKt)/(1.0 + KmCaM/CurrentState[z].cass);
            CaMKa = CaMKb + CurrentState[z].CaMKt;
            NextState[z].CaMKt = CurrentState[z].CaMKt+dt * (aCaMK * CaMKb * (CaMKb + CurrentState[z].CaMKt) - bCaMK * CurrentState[z].CaMKt);
        
            JdiffNa = (CurrentState[z].nass - CurrentState[z].nai)/2.0;
            JdiffK=(CurrentState[z].kss - CurrentState[z].ki)/2.0;
            Jdiff = (CurrentState[z].cass - CurrentState[z].cai)/0.2;
        
            double bt=4.75;
            double a_rel=0.5*bt;
            a_rel *= scaling_coefficients[11];
            double Jrel_inf=a_rel*(-ICaL)/(1.0+pow(1.5/CurrentState[z].cajsr,8.0));
            if (celltype==2)
            {
                Jrel_inf*=1.7;
            }
            double tau_rel=bt/(1.0+0.0123/CurrentState[z].cajsr);
        
            if (tau_rel<0.005)
            {
                tau_rel=0.005;
            }
        
            NextState[z].Jrelnp = Jrel_inf - (Jrel_inf - CurrentState[z].Jrelnp) * exp(-dt/tau_rel);
        
            double btp=1.25*bt;
            double a_relp=0.5*btp;
            a_relp *= scaling_coefficients[11];
        
            double Jrel_infp=a_relp*(-ICaL)/(1.0+pow(1.5/CurrentState[z].cajsr,8.0));
            if (celltype==2)
            {
                Jrel_infp*=1.7;
            }
            double tau_relp=btp/(1.0+0.0123/CurrentState[z].cajsr);
            if (iso==1)
            {
                tau_relp*=0.5;//0.5;
            }
        
            if (tau_relp<0.005)
            {
                tau_relp=0.005;
            }
        
        
            NextState[z].Jrelp = Jrel_infp - (Jrel_infp - CurrentState[z].Jrelp) * exp(-dt/tau_relp);
            double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
            Jrel = (1.0-fJrelp) * CurrentState[z].Jrelnp + fJrelp * CurrentState[z].Jrelp;
        
            double Jupnp;
            double Jupp;
            Jupnp = scaling_coefficients[9]* 0.004375 * CurrentState[z].cai/(CurrentState[z].cai + 0.00092);
            Jupp = scaling_coefficients[9] *2.75 * 0.004375 * CurrentState[z].cai/(CurrentState[z].cai + 0.00092 - 0.00017);
            if (celltype==1)
            {
                Jupnp*=1.3;
                Jupp*=1.3;
            }
        
            double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
            Jleak = 0.0039375 * CurrentState[z].cansr/15.0;
            Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;
            Jtr = (CurrentState[z].cansr - CurrentState[z].cajsr)/100.0;
        
            NextState[z].nai = CurrentState[z].nai +dt * (-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo);
            NextState[z].nass = CurrentState[z].nass +dt * (-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);
            NextState[z].ki = CurrentState[z].ki+ dt * (-(Ito+IKr+IKs+IK1+IKb+Ist-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo);
            NextState[z].kss = CurrentState[z].kss+ dt * (-(ICaK) * Acap/(F * vss) - JdiffK);
        
            double Bcai;
            Bcai=1.0/(1.0+scaling_coefficients[10]*cmdnmax*kmcmdn/pow(kmcmdn + CurrentState[z].cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn + CurrentState[z].cai,2.0));
        
            NextState[z].cai = CurrentState[z].cai+ dt * (Bcai * (-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));
        
            double Bcass=1.0/(1.0+BSRmax*KmBSR/pow(KmBSR + CurrentState[z].cass,2.0)+BSLmax*KmBSL/pow(KmBSL + CurrentState[z].cass,2.0));
            NextState[z].cass = CurrentState[z].cass+dt * (Bcass * (-(ICaL - 2.0 * INaCa_ss) * Acap/(2.0 * F * vss) + Jrel * vjsr/vss - Jdiff));
            NextState[z].cansr = CurrentState[z].cansr+dt * (Jup - Jtr * vjsr/vnsr);
        
            double Bcajsr = 1.0/(1.0 + csqnmax * kmcsqn/pow(kmcsqn + CurrentState[z].cajsr,2.0));
            NextState[z].cajsr = CurrentState[z].cajsr+ dt * (Bcajsr * (Jtr - Jrel));
    }

}

//initial values for state variables, there are 41 of them
void stateInitialization(struct State *CurrentState){
    CurrentState->v         = -87.5;
    CurrentState->nai       = 7.;
    CurrentState->nass      = 7.;
    CurrentState->ki        = 145.;
    CurrentState->kss       = 145.;
    CurrentState->cai       = 1.0e-4;
    CurrentState->cass      = 1.0e-4;
    CurrentState->cansr     = 1.2;
    CurrentState->cajsr     = 1.2;
    CurrentState->m         = 0;
    CurrentState->hf        = 1;
    CurrentState->hs        = 1;
    CurrentState->j         = 1;
    CurrentState->hsp       = 1;
    CurrentState->jp        = 1;
    CurrentState->mL        = 0;
    CurrentState->hL        = 1;
    CurrentState->hLp       = 1;
    CurrentState->a         = 0;
    CurrentState->iF        = 1;
    CurrentState->iS        = 1;
    CurrentState->ap        = 0;
    CurrentState->iFp       = 1;
    CurrentState->iSp       = 1;
    CurrentState->d         = 0;
    CurrentState->ff        = 1;
    CurrentState->fs        = 1;
    CurrentState->fcaf      = 1;
    CurrentState->fcas      = 1;
    CurrentState->jca       = 1;
    CurrentState->nca       = 0;
    CurrentState->ffp       = 1;
    CurrentState->fcafp     = 1;
    CurrentState->xrf       = 0;
    CurrentState->xrs       = 0;
    CurrentState->xs1       = 0;
    CurrentState->xs2       = 0;
    CurrentState->xk1       = 1;
    CurrentState->Jrelnp    = 0;
    CurrentState->Jrelp     = 0;
    CurrentState->CaMKt     = 0;
}


float dv_max(struct State *CurrentState, struct State *NextState, int chain_length)
{
    int z;
    float max=0;
    for (z=0; z<chain_length; z++)
    {
        if(fabs(NextState[z].v-CurrentState[z].v)>max) max=fabs(NextState[z].v-CurrentState[z].v);
    }
    return max;
}

//value holders for state varaibles in the case that the increase in dt was too aggressive, so a smaller one can be taken
double nai0,nass0,ki0,kss0,cai0,cass0,cansr0,cajsr0,m0,hf0,hs0,jO,hsp0,jp0,mL0,hL0,hLp0,a0,iF0,iS0,ap0,iFp0,iSp0,d0,ff0,fs0,fcaf0,fcas0,jca0,nca0,ffp0,fcafp0,xrf0,xrs0,xs10,xs20,xk10,Jrelnp0,Jrelp0,CaMKt0;

void save_state(struct State *initial_state, struct State *LastState)	//O'Hara-Rudy version
{
	initial_state->v=LastState->v;
	initial_state->nai=LastState->nai;
	initial_state->nass=LastState->nass;
	initial_state->ki=LastState->ki;
	initial_state->kss=LastState->kss;
	initial_state->cai=LastState->cai;
	initial_state->cass=LastState->cass;
	initial_state->cansr=LastState->cansr;
	initial_state->cajsr=LastState->cajsr;
	initial_state->m=LastState->m;
    initial_state->hf=LastState->hf;
	initial_state->hs=LastState->hs;
	initial_state->j=LastState->j;
	initial_state->hsp=LastState->hsp;
	initial_state->jp=LastState->jp;
	initial_state->mL=LastState->mL;
	initial_state->hL=LastState->hL;
	initial_state->hLp=LastState->hLp;
	initial_state->a=LastState->a;
	initial_state->iF=LastState->iF;
	initial_state->iS=LastState->iS;
	initial_state->ap=LastState->ap;
	initial_state->iFp=LastState->iFp;
	initial_state->iSp=LastState->iSp;
	initial_state->d=LastState->d;
	initial_state->ff=LastState->ff;
	initial_state->fs=LastState->fs;
	initial_state->fcaf=LastState->fcaf;
	initial_state->fcas=LastState->fcas;
	initial_state->jca=LastState->jca;
	initial_state->nca=LastState->nca;
	initial_state->ffp=LastState->ffp;
	initial_state->fcafp=LastState->fcafp;
	initial_state->xrf=LastState->xrf;
	initial_state->xrs=LastState->xrs;
	initial_state->xs1=LastState->xs1;
	initial_state->xs2=LastState->xs2;
	initial_state->xk1=LastState->xk1;
	initial_state->Jrelnp=LastState->Jrelnp;
	initial_state->Jrelp=LastState->Jrelp;
	initial_state->CaMKt=LastState->CaMKt;

}

int action_potential(struct State *initial_state, double *scaling_coefficients, double *AP, float CL, float amp,const char *statedat_name, int current_time, int iso, int baseline_index, int amount_of_baselines, int amount_of_genes)
{
    const double ft = 1 * CL;
    int chain_length = 1;
    int target_cell = 0;
    
    int z;
    float dv=0;
    float t_output=0;
    struct State *CurrentState, *NextState;
    if ((CurrentState = (struct State*)calloc( chain_length, sizeof(struct State))) == NULL) {printf ("not enough memory for CurrentState\n"); exit(-1);}
    if ((NextState = (struct State*)calloc( chain_length, sizeof(struct State))) == NULL) {printf ("not enough memory for CurrentState\n"); exit(-1);}
    int genes_without_concentrations = amount_of_genes - 2 * amount_of_baselines;
    
    Count = 0;
    t = 0;
    n_stim = 0;
    
    for(z=0;z<chain_length;z++)
    {
        CurrentState[z] = *initial_state;
        NextState[z] = *initial_state;
        
        CurrentState[z].nai = scaling_coefficients[genes_without_concentrations + baseline_index];
        CurrentState[z].nass = scaling_coefficients[genes_without_concentrations + baseline_index];
        CurrentState[z].cansr = scaling_coefficients[genes_without_concentrations + amount_of_baselines + baseline_index];
        
        NextState[z].nai = CurrentState[z].nai;
        NextState[z].nass = CurrentState[z].nass;
        NextState[z].cansr = CurrentState[z].cansr;
    }
    
    vo = CurrentState[0].v;
    
    int cnt = 0;
    while (t<=ft)
    {
        #if ADAPTIVE_STEP
            if (t>=ft-CL && t<ft-CL+current_time)
            {
                if (Count%skip==0)
                {
                    AP[cnt] = CurrentState[target_cell].v;
                    cnt+=1;
                }
                dt=0.005;
                Count++;
            }
            else if (((t>=(start+n_stim*CL-2)) && (t<(start+duration+n_stim*CL)))||(t<(start+duration+(n_stim-1)*CL+safetime)))
            {
                dt=0.005;
            }
            else if((dv>0.08))
            {
                t=t-dt;
                dt=fabs(0.02*dt/dv);
                RGC(CurrentState,NextState, z, chain_length, scaling_coefficients, amp, CL, iso);
                t=t+dt;
                while((dv_max(CurrentState,NextState, chain_length)>0.08)&&(dt>0.005))
                {
                    t=t-dt;
                    dt=dt/2.0;
                    if(dt<0.005) dt=0.005;
                    RGC(CurrentState,NextState, z, chain_length, scaling_coefficients, amp, CL, iso);
                    t=t+dt;
                }
            }
            else if (dv<0.02)
            {
                dt=fabs(0.08*dt/dv);
                if (dt>0.09)
                {
                    dt=0.09;
                }
            }
        #endif
        
        for (z=0; z<chain_length; z++) CurrentState[z]=NextState[z];
        RGC(CurrentState,NextState, z, chain_length, scaling_coefficients, amp, CL, iso);
        t=t+dt;
        
        dv=dv_max(CurrentState,NextState, chain_length);
    }
    
    //Rewriting states!
    save_state(initial_state, &NextState[target_cell]);
    scaling_coefficients[genes_without_concentrations + baseline_index] = initial_state->nai;
    scaling_coefficients[genes_without_concentrations + amount_of_baselines + baseline_index] = initial_state->cansr;
    
    return 0;
}

