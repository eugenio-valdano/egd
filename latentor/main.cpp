//
//  main.cpp
//  latentor
//
//  Created by Eugenio Valdano on 2/13/19.
//  Copyright © 2019 Eugenio Valdano. All rights reserved.
//

#include "typedef.h"
#include "my_random.h"
#include "sr.h"
#include "params.h"
//#include "omp.h"

int main(int argc, const char * argv[]) {
    
    // ============================================
    // ============================================
    // ======= COMMAND LINE PARAMETERS: PARAM FILE.
    // ======= full name to JSON file containing parameters (json file MUST BE UTF-8 encoded)
    
    
    argv++;
    std::string name_param = *argv; // json file
    
    // =======
    // ======= READ PARAMETER FILE
    // =======
    rapidjson::Document doc_json;
    FILE* fp_json = fopen(name_param.c_str(), "r"); // non-Windows use "r"
    char readBuffer[65536];
    rapidjson::FileReadStream is_json(fp_json, readBuffer, sizeof(readBuffer));
    rapidjson::Document d_json;
    assert( ~d_json.ParseStream(is_json).HasParseError() );
    fclose(fp_json);
    assert(d_json.IsObject());
    
    // create vector of Parameters objects
    std::vector<Params> lparams;
    create_parametri(d_json, lparams);
    std::string path_to_network = d_json["network"]["path"].GetString();
    unsigned int T = d_json["network"]["T"].GetInt();
    std::string path_to_output = d_json["output"]["path"].GetString();
    MYDOUBLE TOLERANCE = d_json["computation"]["tolerance"].GetDouble();
    MYDOUBLE vmin = d_json["computation"]["vmin"].GetDouble();
    MYDOUBLE vmax = d_json["computation"]["vmax"].GetDouble();
    
    // load network
    std::vector<Snapshot> network;
    for (unsigned int t=0; t<T; ++t) {
        std::ostringstream fn;
        fn << path_to_network << t << "_";
        network.emplace_back(fn.str());
    }
    
    // output vector
    unsigned int NP = (unsigned int) lparams.size();
    std::vector<MYDOUBLE> output (NP);
    
    // BIG LOOP
    //int n_threads = omp_get_max_threads();
    //std::cout << "used threads= " << n_threads << "\n\n";
    //#pragma omp parallel for default(none) shared(lparams,network,NP,output,vmin,vmax,TOLERANCE) num_threads(n_threads)
    for (unsigned int h = 0; h < NP; ++h) {
        MYDOUBLE lc;
        
        // check if it is always below threshold
        lc = compute_sr(1.0, lparams[h].mu, lparams[h].epsilon, network, lparams[h].latency);
        if (lc < 0.) {
            output[h] = 1.0;
            continue;
        }
        
        lc = -1.;
        
        // boundary estimate
        if (lparams[h].is_sis_estimated) {
            std::cout << "enter in sis estimated\n";
            MYDOUBLE vmin_t = (lparams[h].latency) ? lparams[h].SEIS_lower : lparams[h].SIS_lower;
            MYDOUBLE vmax_t = (lparams[h].latency) ? lparams[h].SEIS_upper : lparams[h].SIS_upper;
            lc = findroot_brent(vmin_t, vmax_t, lparams[h].mu, lparams[h].epsilon, network, lparams[h].latency, 1E-7, TOLERANCE);
        }
        
        // check if no threshold computed in boundary estimate
        if (lc < 0.) {
            std::cout << "need to compute in full range " << lparams[h].latency << " " << lparams[h].mu << " " << lparams[h].epsilon << "\n";
            lc = findroot_brent(vmin, vmax, lparams[h].mu, lparams[h].epsilon, network, lparams[h].latency, 1E-7, TOLERANCE);
        }
        
        // assign result
        output[h] = lc;
    }
    
    // write output
    std::ostringstream fn1;
    fn1 << path_to_output << "output.txt";
    std::ofstream sout1(fn1.str().c_str(), std::ios::trunc);
    assert(sout1.is_open());
    sout1 << "mu index,epsilon index,latency,lc\n";
    for (unsigned int h=0; h < NP; ++h) {
        sout1 << lparams[h].i_mu << "," << lparams[h].i_epsilon << "," << lparams[h].latency << "," << output[h] << "\n";
    }
    sout1.close();
    
    return 0;
}
