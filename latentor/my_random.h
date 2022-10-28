//
//  my_random.h
//  latentor
//
//  Created by Eugenio Valdano on 2/13/19.
//  Copyright © 2019 Eugenio Valdano. All rights reserved.
//

#ifndef my_random_h
#define my_random_h

// +++***+++ +++***+++ +++***+++ +++***+++ +++***+++ +++
// +++***+++ CLASS TO HIDE BOOST SBABBIO, AND HAVE A PRACTICAL OBJECT FOR RANDOM GENERATION +++***+++

// REMEMBER: objects of class MyRandom MUST BE PASSED TO FUNCTIONS BY REFERENCE

// relevant BOOST headers (from https://www.boost.org/)
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/exponential_distribution.hpp>

// USAGE
// MyRandom sbabbio(random_seed);
// sbabbio.rand() // random real in [0,1)
// sbabbio.randint(n) // random integer in [0,n)
// sbabbio.randexp(a) // random real from exponential distribution with rate=a


class MyRandom {

private:

    boost::mt19937_64 eng;
    boost::uniform_int<> dist_int;
    boost::uniform_real<> dist_rand;
    boost::exponential_distribution<> dist_exp;

public:

    // constructor
    MyRandom(const unsigned int);

    // provide generators.
    boost::variate_generator<boost::mt19937_64,boost::uniform_int<>> randint;
    boost::variate_generator<boost::mt19937_64,boost::uniform_real<>> rand;
    boost::variate_generator<boost::mt19937_64,boost::exponential_distribution<>> randexp;


};

// Define constructor, using initialization list syntax
MyRandom::MyRandom(const unsigned int rs) : eng(rs), dist_rand(0,1), randint(eng,dist_int), rand(eng,dist_rand), randexp(eng,dist_exp)  {}




#endif /* my_random_h */
