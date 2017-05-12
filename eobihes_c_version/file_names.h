//
//  files_names.h
//  EOB_ihes
//
//  Created by Philipp Fleig on 23/01/2017.
//  Copyright © 2017 Philipp Fleig. All rights reserved.
//

#ifndef files_names_h
#define files_names_h

#include <iostream>
#include <list>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <string>

#include <vector>

#include "input_struc.h"


using namespace::std;

vector<string> file_names (input *params) {

    int lm = (*params).lm;
    double q = (*params).q;
    double r0 = (*params).r0;
    double dt = (*params).dt;
    double chi1 = (*params).chi1;
    double chi2 = (*params).chi2;
    int solver_scheme = (*params).solver_scheme;

    //Output file definitions
    //ofstream Omeg("data/Omega.dat");
    
     std::string conf;
     if ((*params).tidal==true) {
     conf = "BNS";
     } else {
     conf = "BBH";
     }
     std::string q_str = std::to_string (q);
     q_str.erase ( q_str.find_last_not_of('0') + 1, std::string::npos );
     if (q_str.back()=='.') { q_str.pop_back();}
     std::string r0_str = std::to_string (r0);
     r0_str.erase ( r0_str.find_last_not_of('0') + 1, std::string::npos );
     if (r0_str.back()=='.') {r0_str.pop_back();}
     std::string dt_str = std::to_string (dt);
     dt_str.erase ( dt_str.find_last_not_of('0') + 1, std::string::npos );
     if (dt_str.back()=='.') {dt_str.pop_back();}
     std::string chi1_str = std::to_string (chi1);
     chi1_str.erase ( chi1_str.find_last_not_of('0') + 1, std::string::npos );
     if (chi1_str.back()=='.') {chi1_str.pop_back();}
     std::string chi2_str = std::to_string (chi2);
     chi2_str.erase ( chi2_str.find_last_not_of('0') + 1, std::string::npos );
     if (chi2_str.back()=='.') {chi2_str.pop_back();}
     std::string scheme;
     switch (solver_scheme) {
     case 0:
     scheme="adap";
     break;
     case 1:
     scheme="hybrid";
     break;
     case 2:
     scheme="fixed";
     break;
     default:
     scheme="adap";
     break;
     }

    string data = "data/dynamics_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat";
    string wave("data/h_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");
    string QOmg("data/Qomg_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");
    string Afunc("data/Afunc_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");
    string init("data/initialData_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+".dat");
    string time_step("data/Tstep_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");
    string n1245("data/n1245_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");

    return {data,wave,QOmg,Afunc,init,time_step,n1245};
}


#endif /* files_names_h */
