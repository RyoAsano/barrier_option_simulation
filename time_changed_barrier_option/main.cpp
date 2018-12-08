//
//  main.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/19.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include <iostream>
#include <boost/function.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ql/quantlib.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include "vector_field_class.hpp"
#include "conditional_brownian_motion_generator.hpp"
#include <fstream>
#include "function.hpp"


using namespace boost::numeric::ublas;
using namespace QuantLib;


int main() {
    
    BigInteger seed = SeedGenerator::instance().get();
    MersenneTwisterUniformRng unifMt(seed);
    BoxMullerGaussianRng<MersenneTwisterUniformRng> norm_rand_gen(unifMt);
    
    
    
    std::vector<double> init_list = {1.0, 1.0, 0};
    vector<double> init_point(1);
    
    init_point(0) = 0;
    
    VectorFieldsTest V;
    
    boost::function<double(vector<double>)> f = [](vector<double> x){return x(0)*x(0);};
    vector<double> x(1);
    x(0) = 0;
    double t = 1.0;
    double y = 3.0;
    unsigned long int N = 100;
    int j_star = 1;
    
    unsigned long int MonteCarloSimNum = 100;
    
    double running_sum=0;
    
    
    std::cout << running_sum/MonteCarloSimNum << "\n";

    return 0;
}
                                            
/*
 If you want to use exponential r.v. or bernulli r.v.s you can paste the following.
 
 RandomSequenceGenerator<MersenneTwisterUniformRng> unifMtSeq(1,unifMt);
 boost::function<double(double)> ExpFunc = &ExpInverseCumul;
 boost::function<int(double)> Lambda1Func = &Lambda1GenFunc;
 boost::function<int(double)> Lambda2Func = &Lambda2GenFunc;
 InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<double(double)>> exp_rand_gen(unifMtSeq, ExpFunc);
 InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<int(double)>> Lambda1_rand_gen(unifMtSeq, Lambda1Func);
 InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<int(double)>> Lambda2_rand_gen(unifMtSeq, Lambda2Func);
*/

/*
 If you want an output of conditional_BM_paths, paste the following.
 
 std::ofstream outputfile("/Users/asanoryo/Documents/XcodeProjects/time_changed_barrier_option/time_changed_barrier_option/cond_BM_paths.csv");
 for(int i=0; i< MonteCarloSimNum; i++)
 {
 std::vector<double> one_path = conditional_BM_generator(norm_rand_gen, N, t, y);
 for(int i = 0; i<N-1;i++)
 {
 outputfile << one_path[i] << ",";
 }
 outputfile << one_path[N] << "\n";
 }
 outputfile.close();
*/


/*
 If you want to extract maximum values of paths, paste this:
 std::ofstream outputfile("/Users/asanoryo/Documents/XcodeProjects/time_changed_barrier_option/time_changed_barrier_option/max_and_argmax_list.csv");
 outputfile << "argmax, max\n";
 for(int i=0; i<N;i++)
 {
 std::vector<double> path = conditional_BM_generator(norm_rand_gen, N, t, y);
 std::vector<double>::iterator maxIt = std::max_element(path.begin(), --path.end());
 double max = *maxIt;
 size_t argmax = std::distance(path.begin(), maxIt);
 outputfile << argmax << ", " << max <<"\n";
 }
 outputfile.close();
*/




