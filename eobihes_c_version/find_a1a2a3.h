//
//  find_a1a2a3.h
//  EOB_ihes
//
//  Created by Philipp Fleig on 23/01/2017.
//  Copyright © 2017 Philipp Fleig. All rights reserved.
//

#ifndef find_a1a2a3_h
#define find_a1a2a3_h

#include <stdio.h>
#include <vector>

#include <ios>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "input_struc.h"
#include "multipole_index.h"
#include  "s_D1.h"
//#include "nos_ddotr.h"
//#include "s_ddotr.h"
#include "dtnqc_fit.h"
//#include "get_Omg_orb.h"

#include <limits>

typedef std::numeric_limits< double > dbl;

using namespace::std;

vector<gsl_complex> find_a1a2a3_hlm(vector<double> T,vector<double> r,vector<double> w,vector<double> pph,vector<double> pr_star,vector<double> hlm_phase,std::ofstream& n1245,vector<double> Omg_orb,vector<double> A,vector<double> ddotr,void *params) {

/*
 % determine NQC parameters
 % it uses aprocedure based on
 % i) the matching of the 22 maximum in case of EOB dynamics
 % ii) the matching ad different macthing times in case of EMR dynamics
 % hard-coded parameters (for procedure (ii)
 %    based   : base dir for EMR data = './EMRdata/'
 %    use_tab : table of EMR data
 %    use_tm  : matching time
 % they can be set via global variables
 %    global_EMRvEOB_Set_EMRDataDir
 %    global_EMRvEOB_Set_EMRDataTab
 %    global_EMRvEOB_Set_EMRMatchT
 */

/*
if fitted_nqc_point

%========================================
% all quantities are inizialized to zero.
% this is only limited to the 22 mode.
% TO BE CHECKED/IMPROVED
%========================================
*/

double nu = (*(input *)params).nu;
double chi1 = (*(input *)params).chi1;
double aKerr = (*(input *)params).aK;

    
long int t_length = T.size();

for (int j=t_length;j--;){
    A[j] = A[j]/sqrt(24.); // in general divide by sqrt( (l+2) (l+1) l (l-1) ) Use mutlipole structure to get the correct L
}

vector<double> omg = s_D1(T, hlm_phase, t_length-1);
vector<double> domg = s_D1(T, omg, t_length-1);

/*
 global representation of the NQC-extraction point on
 the NR waveform through fit. (2,2) multipole only

 IMPROVED GLOBAL FIT USING NEW SXS DATA AND CUBIC
 BEHAVIOR
*/
    
double A_tmp, dA_tmp, omg_tmp, domg_tmp;

//        case 'NQC_fit_hybrid' :
            if (nu==0.25) {
                
                //Amplitude
                vector<double> pA(5);
                
                pA[0] = 0.00178195;
                pA[1] = 0.00435589;
                pA[2] = 0.00344489;
                pA[3] = -0.00076165;
                pA[4] = 0.31973334;
                A_tmp = pA[0]*pow(aKerr,4) + pA[1]*pow(aKerr,3) + pA[2]*pow(aKerr,2) + pA[3]*aKerr + pA[4];
                
                //Derivative of amplitude
                vector<double> pdA(5);
                
                pdA[0] = 0.00000927;
                pdA[1] = -0.00024550;
                pdA[2] = 0.00012469;
                pdA[3] = 0.00123845;
                pdA[4] = -0.00195014;
                dA_tmp = pdA[0]*pow(aKerr,4) + pdA[1]*pow(aKerr,3) + pdA[2]*pow(aKerr,2) + pdA[3]*aKerr + pdA[4];
                
                //Frequency
                vector<double> pomg(5);
                
                pomg[0] = 0.00603482;
                pomg[1] = 0.01604555;
                pomg[2] = 0.02290799;
                pomg[3] = 0.07084587;
                pomg[4] = 0.38321834;
                omg_tmp = pomg[0]*pow(aKerr,4) + pomg[1]*pow(aKerr,3) + pomg[2]*pow(aKerr,2) + pomg[3]*aKerr + pomg[4];
                
                //Derivative of the frequency
                vector<double> pdomg(5);
                
                pdomg[0] = 0.00024066;
                pdomg[1] = 0.00038123;
                pdomg[2] = -0.00049714;
                pdomg[3] = 0.00041219;
                pdomg[4] = 0.01190548;
                domg_tmp = pdomg[0]*pow(aKerr,4) + pdomg[1]*pow(aKerr,3) + pdomg[2]*pow(aKerr,2) + pdomg[3]*aKerr + pdomg[4];
                
            } else {
                
                //Amplitude
                vector<double> p1v(2),p2v(2),p3v(2);
                
                p1v[0] = 0.05385059;
                p1v[1] = -0.00890942;
                p2v[0] = -0.07942102;
                p2v[1] = 0.02152423;
                p3v[0] = 0.14805262;
                p3v[1] = 0.28210487;
                double p1    = p1v[0]*nu + p1v[1];
                double p2    = p2v[0]*nu + p2v[1];
                double p3    = p3v[0]*nu + p3v[1];
                A_tmp = p1*pow(aKerr,2) + p2*aKerr + p3;
                
                //Derivative of amplitude
                vector<double> pdA1v(2),pdA2v(2),pdA3v(2);
                
                pdA1v[0] = 0.00248472;
                pdA1v[1] = -0.00033422;
                pdA2v[0] = 0.00105298;
                pdA2v[1] = 0.00085160;
                pdA3v[0] = -0.00339257;
                pdA3v[1] = -0.00110932;
                double pdA1    = pdA1v[0]*nu + pdA1v[1];
                double pdA2    = pdA2v[0]*nu + pdA2v[1];
                double pdA3    = pdA3v[0]*nu + pdA3v[1];
                dA_tmp  = pdA1*pow(aKerr,2) + pdA2*aKerr + pdA3;
                
                //Frequency
                vector<double> pn0(2),pd1(2);
                pn0[0]  = 0.45584139;
                pn0[1]  = 0.27315247;
                pd1[0]  = 0.75276414;
                pd1[1]  = -0.40081625;
                double n0      = pn0[0]*nu + pn0[1];
                double d1      = pd1[0]*nu + pd1[1];
                omg_tmp = n0/(1. + d1*aKerr);

                //Derivative of the frequency
                vector<double> ppdomg1(2),ppdomg2(2);
                ppdomg1[0] = -0.00177362;
                ppdomg1[1] = 0.00123900;
                ppdomg2[0] = 0.02424739;
                ppdomg2[1] = 0.00566504;
                double pdomg1     = ppdomg1[0]*nu + ppdomg1[1];
                double pdomg2     = ppdomg2[0]*nu + ppdomg2[1];
                domg_tmp   = pdomg1*aKerr + pdomg2;
                
            }
                                  
vector<double> max_A(35),max_dA(35),d2max(35),d3max(35),max_omg(35),max_domg(35),maxd2omg(35),DeltaT(35);

for (int k=35;k--;){
    max_A[k]    = 0.;
    max_dA[k]   = 0.;
    max_omg[k]  = 0.;
    max_domg[k] = 0.;
}
// switch on the 22  values (only)
max_A[1]    = A_tmp;
max_dA[1]   = dA_tmp;
max_omg[1]  = omg_tmp;
max_domg[1] = domg_tmp;
    
// NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)

// NQC basis for (2,2) waveform : AMPLITUDE
// note: n3 and n6 are not used

vector<double> n1(t_length), n2(t_length), n3(t_length), n4(t_length), n5(t_length), n6(t_length);
                                  
for (int j=t_length;j--;) {
    n1[j]  = pow(pr_star[j]/(r[j]*w[j]),2);  // [pr*/(r Omg)]^2
    n2[j]  = ddotr[j]/(r[j]*w[j]*w[j]);      // [ddot{r}/(r Omg^2)]
    n3[j]  = n1[j]*pr_star[j]*pr_star[j];            // [pr*/(r Omg)]^2 *(pr*)^2

    //====================================
    // NQC basis for (2,2) waveform: PHASE
    //------------------------------------

    n4[j]  = pr_star[j]/(r[j]*w[j]);           //  pr*/(r Omg)
    n5[j]  = n4[j]*pow(r[j]*w[j],2);           // (pr*)*(r Omg)
    n6[j]  = n5[j]*pow(pr_star[j],2);          // (pr*^3)*(r Omg)
}
                                  
// take derivatives: needed for the phase
vector<double> d_n4  = s_D1(n4,T,t_length-1);
vector<double> d_n5  = s_D1(n5,T,t_length-1);
vector<double> d2_n4 = s_D1(d_n4,T,t_length-1);
vector<double> d2_n5 = s_D1(d_n5,T,t_length-1);

int Omgmax_index = 0;
double Omg_max = Omg_orb[0];
int i = 1;
while (Omg_orb[i] > Omg_max){
    Omg_max = Omg_orb[i];
    Omgmax_index = i;
    //printf(" %s %d %.20e \n","Omg_orb",i,Omg_orb[i]);

    i++;
}
//printf(" %s %d %.20e \n","Omg_orb",Omgmax_index,Omg_orb[Omgmax_index+1]);

double tOmgOrb_pk = T[Omgmax_index];

double DeltaT_nqc = 0.;
if (chi1 >= 0.8498){
    /* Interpolating fit for Deltat_NQC. See Eq.(21) of arXiv:1506.08457 */
    DeltaT_nqc = dtnqc_fit(chi1,0.8498);
}else {
    DeltaT_nqc = 1.;
}

double tNQC    = tOmgOrb_pk - DeltaT_nqc;

//int jmax    = find(T>=tNQC,1,'first');  // IN GENERAL: jmax(l,m+1)!!!
i = 0;
int jmax = 0;
while (T[i] < tNQC){
    jmax = i;
    //printf(" %s %d %.20e %.20e \n","tNQC",i,tNQC,T[i]);
    
    i++;
}

//printf("\n %s %.20e \n","tOmgOrb_pk = ",tOmgOrb_pk);
//printf(" %s %.20e \n ","DeltaT_nqc = ",DeltaT_nqc);
//printf(" %s %.20e \n","tNQC = ",tNQC);
//printf(" %s %.20e \n","T[jmax] = ",T[jmax]);

//==========================================================
// Determination of NQC correction: solving a linear systems
//==========================================================

vector<double> m11(t_length),m12(t_length),m13(t_length);

for (int j=t_length; j--;) {

    // Determination of NQC correction: solving a linear systems
    
    // matrix elements
    m11[j] = n1[j]*A[j]; // waveform amplitude at all points
    m12[j] = n2[j]*A[j];
    
}
vector<double> m21 = s_D1(m11,T,t_length-1);
vector<double> m22 = s_D1(m12,T,t_length-1);

vector<double> p1tmp = A;
vector<double> p2tmp = s_D1(p1tmp,T,t_length-1);

vector<double> ai(3);
vector<double> bi(3);

vector<double> P(2);
vector<double> M(4);
    
double detM = 1.;
//for (int k=35;k--;){
    int k=1;
    
    /* computation of ai coefficients */
    P[0] = max_A[k]   - p1tmp[jmax];
    P[1] = max_dA[k]  - p2tmp[jmax];
    
    M[0] = m11[jmax];
    M[1] = m12[jmax];
    M[2] = m21[jmax];
    M[3] = m22[jmax];
    
    detM = M[0]*M[3]-M[1]*M[2];
    ai[0] = (M[3]*P[0] - M[1]*P[1])/detM;
    ai[1] = (M[0]*P[1] - M[2]*P[0])/detM;
    ai[2] = 0.;
    
    /* computation of bi coefficients */
    P[0] = omg[jmax] - max_omg[k];
    P[1] = domg[jmax]  - max_domg[k];
    
    printf(" %s %.20e \n","omg jmax ",omg[jmax]);
    printf(" %s %.20e \n","domg jmax ",domg[jmax]);
    printf(" %s %.20e \n","max omg ",max_omg[k]);
    printf(" %s %.20e \n","max domg ",max_domg[k]);
    
    M[0] = d_n4[jmax];
    M[1] = d_n5[jmax];
    M[2] = d2_n4[jmax];
    M[3] = d2_n5[jmax];

    detM = M[0]*M[3]-M[1]*M[2];
    bi[0] = (M[3]*P[0] - M[1]*P[1])/detM;
    bi[1] = (M[0]*P[1] - M[2]*P[0])/detM;
    bi[2] = 0.;
//}

vector<gsl_complex> o(t_length);
//for (int k=35;k--;){
    for (int j=t_length; j--;) {
        o[j].dat[0] = 1. + ai[0]*n1[j] + ai[1]*n2[j] + ai[2]*n3[j];
        o[j].dat[1] = bi[0]*n4[j] + bi[1]*n5[j] + bi[2]*n6[j];
        n1245 << T[j] << "\t" << d_n4[j] << "\t" << d_n5[j] << "\t" << d2_n4[j] << "\t" << d2_n5[j] << endl;
    }

    printf("\n %s %.20e \n","ai[0] = ",ai[0]);
    printf(" %s %.20e \n","ai[1] = ",ai[1]);
    printf(" %s %.20e \n","bi[0] = ",bi[0]);
    printf(" %s %.20e \n","bi[1] = ",bi[1]);

//}

return o;

}

#endif /* find_a1a2a3_h */

/* DEFINITION OF MERGER TIME: THE PEAK OF (2,2) MULTIPOLE */
//double tmrg    = tNQC-2.;
//i = 0;
////int jmrg    = find(T>=tmrg,1,'first');
//int jmrg = 0;
//while (T[i] < tmrg){
//    jmrg = i;
//    i++;
//}
////checks of relocating the points //NOT used anywhere! Just for reference here
//double chk_nqc = T[jmax]-tNQC;  // jmax should be called jNQC
//double chk_mrg = T[jmrg]-tmrg; // only for checks

/* s.ell(l).emm(m).A     =  abs(s.ell(l).emm(m).psi);
psi = h/sqrt( (l+2) (l+1) l (l-1) ) */

// 4117.7448

//-0.01139205	   1.50758908	   0.15879853 	   1.81653810

//vector<double> A(t_length);
//vector<double> B(t_length);
//vector<double> metric(4);
//for (int j=t_length;j--;){
//    metric = s_Metric(r[j], params); //{A,B,dA,d2A} data[0]=A; data[1]=A_dr; data[2]=A_du; data[3]=B; data[4]=B_dr;
//    A[j] = metric[0];
//    B[j] = metric[1];
//}
//vector<double> Omg_orb = get_Omg_orb(r, pph, pr_star, A, B, params);


//int wmax_index = 0;
//double wmax = 0.;
//int i = w.size()-1;
//while (w[i] > wmax){
//    wmax = w[i];
//    i--;
//    wmax_index = i;
//}
//double tOmg_pk        = T[wmax_index];
// time difference between peaks
//double DtOmg_pk     = tOmgOrb_pk-tOmg_pk;


//vector<double> ddotr(t_length);
//if (S1!=0. && S2!=0.){
//    //ddotr = EOB_ddotrSpin(T,r,pph,pr_star);
//    for(int j=t_length;j--;){
//        ddotr[j] = s_ddotr(T[j],r[j],pph[j],pr_star[j],params);
//    }
//} else {
//    //ddotr = EOB_ddotr(T,r,pph,pr_star);
//    for(int j=t_length;j--;){
//        ddotr[j] = nos_ddotr(T[j],r[j],pph[j],pr_star[j],params);
//    }
//}


//    switch (<#expression#>) {
//        case ... : //this case is not relevant for the moment
//            /*  */
//
//            vector<double> p1(2),p2(2),p3(2),p4(2);
//
//            p1[0] = 0.04680896;
//            p1[1] = -0.00632114;
//            p2[0] = 0.06586192;
//            p2[1] = -0.01180039;
//            p3[0] = -0.11617413;
//            p3[1] = 0.02704959;
//            p4[0] = 0.15597465;
//            p4[1] = 0.28034978;
//            double p1    = p1[0]*nu + p1[1];
//            double p2    = p2[0]*nu + p2[1];
//            double p3    = p3[0]*nu + p3[1];
//            double p4    = p4[0]*nu + p4[1];
//            A_tmp     = p1*pow(aKerr,3) + p2*pow(aKerr,2) + p3*aKerr + p4;
//
//            //Derivative of amplitude
//
//            vector<double> pdA1(2),pdA2(2),pdA3(2),pdA4(2);
//
//            pdA1[0] = -0.00130824;
//            pdA1[1] = 0.00006202;
//            pdA2[0] = 0.00199855;
//            pdA2[1] = -0.00027474;
//            pdA3[0] = 0.00218838;
//            pdA3[1] = 0.00071540;
//            pdA4[0] = -0.00362779;
//            pdA4[1] = -0.00105397;
//            pdA1   = pdA1[0]*nu + pdA1[1];
//            pdA2   = pdA2[0]*nu + pdA2[1];
//            pdA3   = pdA3[0]*nu + pdA3[1];
//            pdA4   = pdA4[0]*nu + pdA4[1];
//            dA_tmp = pdA1*pow(aKerr,3) + pdA2*pow(aKerr,2) + pdA3*aKerr + pdA4;
//
//            //Frequency
//            vector<double> pn0(2), pd1(2);
//
//            pn0[0]  = 0.46908067;
//            pn0[1]  = 0.27022141;
//            pd1[0]  = 0.64131115;
//            pd1[1]  = -0.37878384;
//            double n0      = pn0[0]*nu + pn0[1];
//            double d1      = pd1[0]*nu + pd1[1];
//            omg_tmp = n0/(1. + d1*aKerr);
//
//            //Derivative of the frequency
//            vector<double> ppdomg1(2), ppdomg2(2);
//
//            ppdomg1[0] = 0.00061175;
//            ppdomg1[1] = 0.00074001;
//            ppdomg2[0] = 0.02504442;
//            ppdomg2[1] = 0.00548217;
//            double pdomg1   = ppdomg1[0]*nu + ppdomg1[1];
//            double pdomg2   = ppdomg2[0]*nu + ppdomg2[1];
//            domg_tmp = pdomg1*aKerr + pdomg2;
//            break;

//....

//            break;
//
//        default:
//            break;
//    }
//



