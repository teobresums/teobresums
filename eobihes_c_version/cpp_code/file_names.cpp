/**
 * Copyright (C) 2017 Alessandro Nagar, Gregorio Carullo, Ka Wa Tsang, Philipp Fleig, Sebastiano Bernuzzi, Walter Del Pozzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <iostream>
#include <list>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <string>
#include <vector>
#include <tuple>
#include "input_struc.h"

#include "file_names.h"
#include "multipole_index.h"

using namespace::std;


tuple<vector<string>, vector<string> > file_names (input *params)
{
    int solver_scheme = (*params).solver_scheme;
    int lm            = (*params).lm;
    double q          = (*params).q;
    double r0         = (*params).r0;
    double dt         = (*params).dt;
    double chi1       = (*params).chi1;
    double chi2       = (*params).chi2;

    vector<string> wave(35);
    std::string conf;
    std::string scheme;
    std::string q_str    = std::to_string (q);
    std::string r0_str   = std::to_string (r0);
    std::string dt_str   = std::to_string (dt);
    std::string chi1_str = std::to_string (chi1);
    std::string chi2_str = std::to_string (chi2);

    /**Output file definitions*/
    //ofstream Omeg("data/Omega.dat");
    
     if ((*params).tidal==true){conf = "BNS";}
     else {conf = "BBH";}
    
     q_str.erase ( q_str.find_last_not_of('0') + 1, std::string::npos );
     if (q_str.back()=='.') { q_str.pop_back();}
     r0_str.erase ( r0_str.find_last_not_of('0') + 1, std::string::npos );
     if (r0_str.back()=='.') {r0_str.pop_back();}
     dt_str.erase ( dt_str.find_last_not_of('0') + 1, std::string::npos );
     if (dt_str.back()=='.') {dt_str.pop_back();}
     chi1_str.erase ( chi1_str.find_last_not_of('0') + 1, std::string::npos );
     if (chi1_str.back()=='.') {chi1_str.pop_back();}
     chi2_str.erase ( chi2_str.find_last_not_of('0') + 1, std::string::npos );
     if (chi2_str.back()=='.') {chi2_str.pop_back();}
    
     switch (solver_scheme)
     {
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
    for (int k=35; k--; )
    {
        wave[k] = "data/h_"+conf+"_"+std::to_string(L[k])+std::to_string(M[k])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat";
    }
    //string wave("data/h_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");
    string QOmg("data/Qomg_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");
    string Afunc("data/Afunc_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");
    string init("data/initialData_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+".dat");
    string time_step("data/Tstep_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");
    string runReport("data/runReport_"+conf+"_"+std::to_string(L[lm])+std::to_string(M[lm])+"_q"+q_str+"_r0_"+r0_str+"_dt"+dt_str+"_chi1_"+chi1_str+"_chi2_"+chi2_str+"_"+scheme+".dat");

    //return {data,wave,QOmg,Afunc,init,time_step,runReport};
    vector<string> files = {data,init,QOmg,Afunc,time_step,runReport,"0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"};
    return std::make_tuple(files, wave);
}


