//
//  sr.h
//  latentor
//
//  Created by Eugenio Valdano on 2/13/19.
//  Copyright © 2019 Eugenio Valdano. All rights reserved.
//

#ifndef sr_h
#define sr_h

typedef long double MYDOUBLE;

// ++++++++++++++++++++++++++
// Class for network handling
// ++++++++++++++++++++++++++

class Snapshot {
    //private:
    
public:
    std::vector<MYDOUBLE> data;
    std::vector<NODE> indices;
    std::vector<NODE> indptr;
    NODE N;
    
    Snapshot (std::string fname);
    void dot (std::vector<MYDOUBLE> & v_input, std::vector<MYDOUBLE> & v_output);  // vnew_i = sum_j A_ij vold_j
    void dot_seis (std::vector<MYDOUBLE> & v_up_input, std::vector<MYDOUBLE> & v_up_output,
                   std::vector<MYDOUBLE> & v_do_input, std::vector<MYDOUBLE> & v_do_output,
                   MYDOUBLE ladda, MYDOUBLE mu, MYDOUBLE epsilon);
    void dot_sis (std::vector<MYDOUBLE> & v_input, std::vector<MYDOUBLE> & v_output, MYDOUBLE ladda, MYDOUBLE mu);
};

Snapshot::Snapshot(std::string finput) { // finput is path to the folder where indices.txt, indptr.txt, data.txt are (you must include last slash).
    
    // +++++
    // read CSR format
    
    // read indptr
    std::ostringstream fn;
    fn << finput << "indptr.txt";
    std::ifstream rin (fn.str().c_str(), std::ios::in);
    assert(rin.is_open());
    NODE x;
    while (true) {
        rin >> x;
        if (rin.eof()) {
            break;
        }
        indptr.push_back(x);
    }
    rin.close();
    
    // number of nodes
    N = NODE(indptr.size())-1;
    
    // read indices
    std::ostringstream fn2;
    fn2 << finput << "indices.txt";
    std::ifstream rin2 (fn2.str().c_str(), std::ios::in);
    assert(rin2.is_open());
    while (true) {
        rin2 >> x;
        if (rin2.eof()) {
            break;
        }
        indices.push_back(x);
    }
    rin2.close();
    
    // read data
    std::ostringstream fn3;
    fn3 << finput << "data.txt";
    std::ifstream rin3 (fn3.str().c_str(), std::ios::in);
    assert(rin3.is_open());
    MYDOUBLE y;
    while (true) {
        rin3 >> y;
        if (rin3.eof()) {
            break;
        }
        data.push_back(y);
    }
    rin3.close();
    
}

void Snapshot::dot(std::vector<MYDOUBLE> & v_input, std::vector<MYDOUBLE> & v_output) {
    
    for (NODE i=0; i<N; ++i) {
        v_output[i] = 0.;
        for (NODE h=indptr[i]; h<indptr[i+1]; ++h) {
            //j = indices[h]
            // A_ij = data[h]
            v_output[i] += data[h] * v_input[indices[h]];
        } // end for h
    }// end for i
}

void Snapshot::dot_seis(std::vector<MYDOUBLE> & v_up_input, std::vector<MYDOUBLE> & v_up_output, std::vector<MYDOUBLE> & v_do_input, std::vector<MYDOUBLE> & v_do_output, MYDOUBLE ladda, MYDOUBLE mu, MYDOUBLE epsilon) {
    
    // matrix M:
    // 1-epsilon | ladda A
    // epsilon   | 1-mu
    
    // UP, MATRIX PART
    dot(v_do_input, v_up_output);
    
    // DOWN whole, and UP multiply by lambda and add other part
    for (NODE i=0; i<N; ++i) {
        v_up_output[i] = ladda * v_up_output[i] + (1.-epsilon) * v_up_input[i]; // multiply by lambda, then add 1-epsilon part
        v_do_output[i] = epsilon * v_up_input[i] + (1.-mu) * v_do_input[i];
    }
}

void Snapshot::dot_sis(std::vector<MYDOUBLE> &v_input, std::vector<MYDOUBLE> &v_output, MYDOUBLE ladda, MYDOUBLE mu) {
    
    //  (1-mu + lambda*A ).v
    
    dot(v_input, v_output); // A.v
    
    for (NODE i=0; i<N; ++i) {
        v_output[i] = ladda * v_output[i] + (1.-mu) * v_input[i]; // multiply by lambda, then add 1-mu part
    }
}

// ****
// **** COMPUTE SPECTRAL RADIUS

// dot product
template <class T>
T dot(std::vector<T> &a, std::vector<T> &b) {
    T ecco = 0;
    for (unsigned int i=0; i < a.size(); ++i) {
        ecco += a[i]*b[i];
    }
    return ecco;
}

// normalize vector to =1
void normalize(std::vector<MYDOUBLE> &a) {
    MYDOUBLE sum = 0.;
    for (unsigned int i = 0; i < a.size(); ++i) {
        sum += a[i]*a[i];
    }
    
    sum = 1.0 / sqrt(sum);
    for (unsigned int i = 0; i < a.size(); ++i) {
        a[i] = a[i] * sum;
    }
}

// normalize vector composed of two vectors to be stacked
void normalize_combo(std::vector<MYDOUBLE> &a, std::vector<MYDOUBLE> &b) {
    MYDOUBLE sum = 0.;
    for (unsigned int i = 0; i < a.size(); ++i) {
        sum += a[i]*a[i];
    }
    for (unsigned int i = 0; i < b.size(); ++i) {
        sum += b[i]*b[i];
    }
    
    sum = 1.0 / sqrt(sum);
    for (unsigned int i = 0; i < a.size(); ++i) {
        a[i] = a[i] * sum;
    }
    for (unsigned int i = 0; i < b.size(); ++i) {
        b[i] = b[i] * sum;
    }
    
}

MYDOUBLE distanceL1(std::vector<MYDOUBLE> &a, std::vector<MYDOUBLE> &b) {
    MYDOUBLE ecco = 0.;
    for (unsigned int i=0; i < a.size(); ++i) {
        ecco += fabs(a[i]-b[i]);
    }
    return ecco;
}



MYDOUBLE compute_sr_sis(MYDOUBLE ladda, MYDOUBLE mu, std::vector<Snapshot> &network,
                        unsigned int valumax=20000, double tolerance=1E-6, unsigned int store=10, MYDOUBLE sr_target=1.0) {
    
    /*
     network is a vector of snapshots
     */
    
    NODE T = NODE( network.size() );
    MYDOUBLE rootT = 1.0 / MYDOUBLE(T);
    NODE N = network[0].N;
    MYDOUBLE fN = 1.0 / MYDOUBLE(N);
    double delta;
    
    std::vector<std::vector<MYDOUBLE>> MV (store, std::vector<MYDOUBLE> (N)); // MV[h, i] with h=store, i=node
    std::vector<MYDOUBLE> v (N);
    std::vector<MYDOUBLE> vnew (N);
    unsigned int c = 0; // register counter
    MYDOUBLE sr; // result
    
    // initialize first vector v, and first vector in register, as normalized foffo v_i = 1/sqrt(N)
    for (NODE i=0; i<N; ++i) {
        MV[0][i] = fN;
        v[i] = MV[0][i];
    }
    
    // big loop
    for (unsigned int k=0; k<T*valumax; ++k) {
        
        network[k%T].dot_sis(v, vnew, ladda, mu); // vnew = (1-mu)v + ladda A.v
        
        if (k % T == T-1) {
            // enter if period has ended
            
            sr = dot(vnew, MV[c%store]);  // spectral radius estimate
            normalize(vnew); // normalize eigenvector
            
            delta = distanceL1(vnew, MV[c%store]);  // compute L1 distance between eigenvectors
            if (delta < tolerance) {
                return pow(sr, rootT) - sr_target;
            }
            
            c++;  // increment register counter
            std::copy(vnew.begin(), vnew.end(), MV[c%store].begin()); // MV[c%store] = vnew
        } // end if k%T
        
        std::copy(vnew.begin(), vnew.end(), v.begin()); // v = vnew
        
    } // end for k
    
    assert(false); // if here, it means that the algorithm has not converged
    return 1E12;
}

MYDOUBLE compute_sr_seis(MYDOUBLE ladda, MYDOUBLE mu, MYDOUBLE epsilon, std::vector<Snapshot> &network,
                         unsigned int valumax=20000, double tolerance=1E-6, unsigned int store=10) {
    
    /*
     network is a vector of snapshots
     */
    
    NODE T = NODE( network.size() );
    MYDOUBLE rootT = 1.0 / MYDOUBLE(T);
    NODE N = network[0].N;
    MYDOUBLE fN = 1.0 / MYDOUBLE(2 * N);
    double delta;
    
    std::vector<std::vector<MYDOUBLE>> MV_up (store, std::vector<MYDOUBLE> (N)); // MV[h, i] with h=store, i=node
    std::vector<std::vector<MYDOUBLE>> MV_do (store, std::vector<MYDOUBLE> (N)); // MV[h, i] with h=store, i=node
    std::vector<MYDOUBLE> v_up (N);
    std::vector<MYDOUBLE> vnew_up (N);
    std::vector<MYDOUBLE> v_do (N);
    std::vector<MYDOUBLE> vnew_do (N);
    unsigned int c = 0; // register counter
    MYDOUBLE sr; // result
    
    // initialize first vector v, and first vector in register, as normalized foffo v_i = 1/sqrt(N)
    for (NODE i=0; i<N; ++i) {
        MV_up[0][i] = fN;
        MV_do[0][i] = fN;
        v_up[i] = MV_up[0][i];
        v_do[i] = MV_do[0][i];
    }
    
    // big loop
    for (unsigned int k=0; k<T*valumax; ++k) {
        
        //network[k%T].dot_sis(v, vnew, ladda, mu); // vnew = (1-mu)v + ladda A.v
        network[k%T].dot_seis(v_up, vnew_up, v_do, vnew_do, ladda, mu, epsilon);
        
        if (k % T == T-1) {
            // enter if period has ended
            
            sr = dot(vnew_up, MV_up[c%store]) + dot(vnew_do, MV_do[c%store]);  // spectral radius estimate
            normalize_combo(vnew_up, vnew_do); // normalize eigenvector
            
            delta = distanceL1(vnew_up, MV_up[c%store]) + distanceL1(vnew_do, MV_do[c%store]);  // compute L1 distance between eigenvectors
            if (delta < tolerance) {
                return pow(sr, rootT) - 1.0;
            }
            
            c++;  // increment register counter
            std::copy(vnew_up.begin(), vnew_up.end(), MV_up[c%store].begin()); // MV_up[c%store] = vnew_up
            std::copy(vnew_do.begin(), vnew_do.end(), MV_do[c%store].begin()); // MV_do[c%store] = vnew_do
        } // end if k%T
        
        std::copy(vnew_up.begin(), vnew_up.end(), v_up.begin()); //cpp  v_up = vnew_up
        std::copy(vnew_do.begin(), vnew_do.end(), v_do.begin()); // v_do = vnew_do
        
    } // end for k
    
    assert(false); // if here, it means that the algorithm has not converged
    return 1E12;
}

// *****
// ***** sr wrapper

MYDOUBLE compute_sr(MYDOUBLE ladda, MYDOUBLE mu, MYDOUBLE epsilon, std::vector<Snapshot> &network, bool latency,
                    unsigned int valumax=20000, double tolerance=1E-6, unsigned int store=10) {
    if (latency) {
        return compute_sr_seis(ladda, mu, epsilon, network, valumax, tolerance, store);
    } else {
        return compute_sr_sis(ladda, mu, network, valumax, tolerance, store);
    }
}

MYDOUBLE findroot_bisect(MYDOUBLE vmin, MYDOUBLE vmax, MYDOUBLE mu, MYDOUBLE epsilon, std::vector<Snapshot> &network, bool latency, MYDOUBLE tol_y=1e-7, MYDOUBLE tol_x=1e-7, unsigned int MAX_ITER = 1000) {
    
    /*
     tol_x: tolerance on the (absolute) width of the interval [a, b].
     tol_y: tolerance on the (absolute) value of the function.
     */
    
    MYDOUBLE a = vmin;
    MYDOUBLE b = vmax;
    MYDOUBLE fa = compute_sr(a, mu, epsilon, network, latency);
    MYDOUBLE fb = compute_sr(b, mu, epsilon, network, latency);
    MYDOUBLE c, fc;
    
    // if wrong order, swap
    if (a > b) {
        std::swap(a,b);
        std::swap(fa,fb);
    }
    
    unsigned int calls = 0;
    
    // if already on a zero
    if (fabs(fa) < tol_y) {
        return a;
    } else if (fabs(fb) < tol_y) {
        return b;
    }
    
    assert (fa*fb < 0.);
    
    while (calls < MAX_ITER) {
        
        c = 0.5 * (a+b);  // new point
        calls++;
        
        if (b-a < tol_x) {  // check if [a,b] small enough // (b-a)/fabs(c) < tol_x
            //std::cout << "calls: " << calls << "\n";
            return c;
        }
        
        fc = compute_sr(c, mu, epsilon, network, latency);
        
        if (fabs(fc) < tol_y) { // check if f is small enough
            //std::cout << "calls: " << calls << "\n";
            return c;
        } else if (fa * fc < 0.) { // choose side
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    } // end while true
    
    assert(false);
    return -1.;
    
} // end function bisect


MYDOUBLE findroot_brent(MYDOUBLE vmin, MYDOUBLE vmax, MYDOUBLE mu, MYDOUBLE epsilon, std::vector<Snapshot> &network, bool latency, MYDOUBLE tol_y=1e-7, MYDOUBLE tol_x=1e-7, unsigned int MAX_ITER = 1000) {
    
    MYDOUBLE a = vmin;
    MYDOUBLE b = vmax;
    MYDOUBLE fa = compute_sr(a, mu, epsilon, network, latency);   // calculated now to save function calls
    MYDOUBLE fb = compute_sr(b, mu, epsilon, network, latency);   // calculated now to save function calls
    MYDOUBLE fs = 0;      // initialize
    
    unsigned int calls = 0;
    
    // if already on a zero
    if (fabs(fa) < tol_y) {
        return a;
    } else if (fabs(fb) < tol_y) {
        return b;
    }
    
    // check that the zero is inside, otherwise return -1
    if (fa*fb > 0.) {
        return -1.;
    }
    //assert (fa*fb < 0.);
    
    if (std::abs(fa) < std::abs(b)) // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
    {
        std::swap(a,b);
        std::swap(fa,fb);
    }
    
    MYDOUBLE c = a;           // c now equals the largest magnitude of the lower and upper bounds
    MYDOUBLE fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
    bool mflag = true;      // boolean flag used to evaluate if statement later on
    MYDOUBLE s = 0;           // Our Root that will be returned
    MYDOUBLE d = 0;           // Only used if mflag is unset (mflag == false)
    
    for (unsigned int iter = 0; iter < MAX_ITER; ++iter)
    {
        // stop if converged on root or error is less than tolerance
        if (std::abs(b-a) < tol_x)
        {
            std::cout << "calls (brent): " << calls << "\n";
            return s;
        } // end if
        
        if (fa != fc && fb != fc)
        {
            // use inverse quadratic interopolation
            s =   ( a * fb * fc / ((fa - fb) * (fa - fc)) )
            + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
            + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
        }
        else
        {
            // secant method
            s = b - fb * (b - a) / (fb - fa);
        }
        
        /*
         Crazy condition statement!:
         -------------------------------------------------------
         (condition 1) s is not between  (3a+b)/4  and b or
         (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
         (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
         (condition 4) (mflag is set and |b−c| < |TOL|) or
         (condition 5) (mflag is false and |c−d| < |TOL|)
         */
        if (    ( (s < (3 * a + b) * 0.25) || (s > b) ) ||
            ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
            ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
            ( mflag && (std::abs(b-c) < tol_x) ) ||
            ( !mflag && (std::abs(c-d) < tol_x))  )
        {
            // bisection method
            s = (a+b)*0.5;
            
            mflag = true;
        }
        else
        {
            mflag = false;
        }
        
        fs = compute_sr(s, mu, epsilon, network, latency);  // calculate fs
        calls++;
        
        d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
        c = b;      // set c equal to upper bound
        fc = fb;    // set f(c) = f(b)
        
        if ( fa * fs < 0)   // fa and fs have opposite signs
        {
            b = s;
            fb = fs;    // set f(b) = f(s)
        }
        else
        {
            a = s;
            fa = fs;    // set f(a) = f(s)
        }
        
        if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
        {
            std::swap(a,b);     // swap a and b
            std::swap(fa,fb);   // make sure f(a) and f(b) are correct after swap
        }
        
    } // end for
    
    assert(false);
    
} // end brent_fun



#endif /* sr_h */



