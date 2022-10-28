//
//  params.h
//  latentor
//
//  Created by Eugenio Valdano on 2/14/19.
//  Copyright © 2019 Eugenio Valdano. All rights reserved.
//

#ifndef params_h
#define params_h


// function for computing values out of start,stop,n
std::vector<MYDOUBLE> build_value_array(MYDOUBLE start, MYDOUBLE stop, int n, bool logscale) {
    
    MYDOUBLE ek;
    std::vector<MYDOUBLE> y;
    
    if (logscale) {
        start = log10(start);
        stop = log10(stop);
    }
    
    if (n==1) {
        if (logscale) {
            start = pow(10., start);
        }
        y.push_back(start);
    } else {
        for (int i=0; i<n; ++i) {
            ek = start + i*(stop-start)/MYDOUBLE(n-1);
            if (logscale) {
                ek = pow(10.0, ek);
            }
            y.push_back(ek);
        }
    }
    
    return y;
}


class Params {

public:
    MYDOUBLE mu, epsilon, tau_mu, tau_epsilon;
    unsigned int i_mu, i_epsilon;
    bool latency;
    bool is_sis_estimated;
    Params(double tm, double te, unsigned int im, unsigned int ie, bool lat, bool is_sis_est, MYDOUBLE a, MYDOUBLE b, MYDOUBLE c);
    MYDOUBLE SIS_estimate, SIS_lower, SIS_upper, SEIS_lower, SEIS_upper;
};

Params::Params(double tm, double te, unsigned int im, unsigned int ie, bool lat, bool is_sis_est, MYDOUBLE a, MYDOUBLE b, MYDOUBLE c) {
    tau_mu = tm;
    tau_epsilon = te;
    mu =  1./tau_mu;
    epsilon = 1./tau_epsilon;
    latency = lat;
    i_mu = im;
    i_epsilon = ie;
    is_sis_estimated = is_sis_est;
    
    // estimate SIS & SEIS bounds
    if (is_sis_estimated) {
        double z = 1.0 / log(tau_mu);
        
        SIS_estimate = a + b * z + c * z * z;
        SIS_lower = 0.85 * SIS_estimate;
        SIS_upper = 1.15 * SIS_estimate;
        
        if (tau_mu < 30.) {
            SEIS_lower = (tau_epsilon < 14.) ? 0.92 * SIS_estimate : 0.99 * SIS_estimate;
            SEIS_upper = (tau_epsilon < 14.) ? 1.15 * SIS_estimate : 1.9 * SIS_estimate;
        } else if (tau_mu < 30.*6) {
            SEIS_lower = 0.99 * SIS_estimate;
            SEIS_upper = (tau_epsilon < 14.) ? 1.2 * SIS_estimate : 1.9 * SIS_estimate;
        } else {
            SEIS_lower = 0.99 * SIS_estimate;
            SEIS_upper = (tau_epsilon < 14.) ? 1.1 * SIS_estimate : 1.4 * SIS_estimate;
        }
        
        // clip to 1
        SIS_upper = (SIS_upper) > 1. ? 1. : SIS_upper;
        SEIS_upper = (SEIS_upper) > 1. ? 1. : SEIS_upper;
        
    } else {
        SIS_estimate = -1.0;
        SEIS_lower = -1.0;
        SEIS_upper = -1.0;
    }
}






void create_parametri(rapidjson::Document & d_json, std::vector<Params> &lparams) {
    
    MYDOUBLE taumu_start = d_json["infectious period"]["start"].GetDouble();
    MYDOUBLE taumu_stop = d_json["infectious period"]["stop"].GetDouble();
    MYDOUBLE taumu_n = d_json["infectious period"]["steps"].GetInt();
    bool taumu_logscale = d_json["infectious period"]["log scale"].GetBool();
    MYDOUBLE tauep_start = d_json["latency period"]["start"].GetDouble();
    MYDOUBLE tauep_stop = d_json["latency period"]["stop"].GetDouble();
    MYDOUBLE tauep_n = d_json["latency period"]["steps"].GetInt();
    bool tauep_logscale = d_json["latency period"]["log scale"].GetBool();
    
    std::string path_to_output = d_json["output"]["path"].GetString();
    
    bool with_SIS = d_json["computation"]["SIS"].GetBool();
    bool with_SEIS = d_json["computation"]["SEIS"].GetBool();
    
    // SIS estimate
    bool is_sis_estimated = d_json["computation"]["SIS estimate"]["estimate from fit"].GetBool();
    MYDOUBLE sis_a = d_json["computation"]["SIS estimate"]["fit parameters"]["deg 0"].GetDouble();
    MYDOUBLE sis_b = d_json["computation"]["SIS estimate"]["fit parameters"]["deg 1"].GetDouble();
    MYDOUBLE sis_c = d_json["computation"]["SIS estimate"]["fit parameters"]["deg 2"].GetDouble();
    
    std::vector<MYDOUBLE> l_taumu = build_value_array(taumu_start, taumu_stop, taumu_n, taumu_logscale);
    std::vector<MYDOUBLE> l_tauep = build_value_array(tauep_start, tauep_stop, tauep_n, tauep_logscale);
    
    if (with_SIS) {
        for (unsigned int h=0; h < l_taumu.size(); ++h) {
            lparams.emplace_back(l_taumu[h], 0.5, h, 0, false, is_sis_estimated, sis_a, sis_b, sis_c);
        }
    }
    
    if (with_SEIS) {
        for (unsigned int h=0; h < l_taumu.size(); ++h) {
            for (unsigned int k=0; k < l_tauep.size(); ++k) {
                lparams.emplace_back(l_taumu[h], l_tauep[k], h, k, true, is_sis_estimated, sis_a, sis_b, sis_c);
            }
        }
    }
    
    // write param vectors to files
    std::ostringstream fn1;
    fn1 << path_to_output << "tau_mu.txt";
    std::ofstream sout1(fn1.str().c_str(), std::ios::trunc);
    assert(sout1.is_open());
    sout1 << "mu index,tau\n";
    for (unsigned int h=0; h < l_taumu.size(); ++h) {
        sout1 << h << "," << l_taumu[h] << "\n";
    }
    sout1.close();
    
    std::ostringstream fn2;
    fn2 << path_to_output << "tau_ep.txt";
    std::ofstream sout2(fn2.str().c_str(), std::ios::trunc);
    assert(sout2.is_open());
    sout2 << "epsilon index,tau\n";
    for (unsigned int h=0; h < l_tauep.size(); ++h) {
        sout2 << h << "," << l_tauep[h] << "\n";
    }
    sout2.close();
    
}


#endif /* params_h */
