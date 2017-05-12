//
//  Omg_orb.h
//  EOB_IHES_cpp_bitbucket
//
//  Created by Philipp Fleig on 14/04/2017.
//  Copyright © 2017 Philipp Fleig. All rights reserved.
//

#ifndef Omg_orb_h
#define Omg_orb_h

vector<double> get_Omg_orb(vector<double> r,vector<double> pph,vector<double> pr_star,vector<double> A,vector<double> B,void *params){

    double nu = (*(input *)params).nu;
    double aK2 = (*(input *)params).aK2;
    double S1 = (*(input *)params).S1;
    double S2 = (*(input *)params).S2;
    double a1 = (*(input *)params).a1;
    double a2 = (*(input *)params).a2;
    double X1 = (*(input *)params).X1;
    double X2 = (*(input *)params).X2;
    double chi1 = (*(input *)params).chi1;
    double chi2 = (*(input *)params).chi2;
    double c3 = (*(input *)params).cN3LO;
    
    double S     = S1 + S2;
    double Sstar = X2*a1 + X1*a2;
    double z3 = 2.*nu*(4.-3.*nu);

    long int r_length = r.size();
    vector<double> Omg_orb(r_length);
    for (int i=r_length; i--; ) {

        vector<double> rc_vec = s_get_rc(r[i],aK2,params);
        double rc = rc_vec[0];
        double drc_dr = rc_vec[1];
        double uc              = 1./rc;
        double uc2             = uc*uc;

        vector<double> ggm = s_GS(r[i],rc,drc_dr,aK2,pr_star[i],pph[i],nu,chi1,chi2,X1,X2,c3);
        double GS = ggm[2];
        double GSs = ggm[3];
        
        double prstar2 = pr_star[i]*pr_star[i];
        double prstar4 = prstar2*prstar2;
        double pphi2 =pph[i]*pph[i];
        double Horbeff=sqrt(prstar2+A[i]*(1. + pphi2*uc2 +  z3*prstar4*uc2));

        double Heff = Horbeff + (GS*S + GSs*Sstar)*pph[i];
        double H = sqrt( 1. + 2.*nu*(Heff - 1.) );
        double one_H = 1./H;

        Omg_orb[i] = one_H*pph[i]*A[i]*uc2/Horbeff;
    }

    return Omg_orb;
    
}

#endif /* Omg_orb_h */
