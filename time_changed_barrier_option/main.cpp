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
#include "generator_Euler_Maruyama_conditional.hpp"
#include <fstream>
#include "function.hpp"
#include "generator_rand_num.hpp"
#include "generator_first_hitting_time.hpp"
#include "generator_Euler_Maruyama.hpp"
#include "generator_Euler_Maruyama_merger.hpp"
#include "BlackScholesTruePrice.hpp"


using namespace boost::numeric::ublas;
using namespace QuantLib;


int main() {
    
    BigInteger seed = SeedGenerator::instance().get();
    MersenneTwisterUniformRng unifMt(seed);
    BoxMullerGaussianRng<MersenneTwisterUniformRng> norm_rand_gen(unifMt);
    
    
    
    vector<double> x(2);

    x(0) = 100;
    x(1) = 0;
    
    vector<double> y(1);
    
    y(0) = 100;
    
    double mu = 0.02;
    double sigma = 0.2;
    double barrier_level = 110;
    double T = 1.0;
    double K = 100;
    
    VectorFieldsTimeChangedBlackScholesWithUpperBarrier V_BS_time_changed(mu, sigma, barrier_level);
    VectorFieldsBlackScholes V_BS_normal(mu, sigma);
    VectorFieldsTimeChangedTest V_time_changed_test;
    VectorFieldsTest V_test;
    
    boost::function<double(vector<double>)> f = [=](vector<double> x){return (K < x(0)) ? (x(0)-K) : 0 ;};
    boost::function<double(vector<double>)> id = [=](vector<double> x){return x(0) ;};
    
    unsigned int j_star = 1;
    
    unsigned long int N_EM_cond = 100;
    unsigned long int N_EM = 10;
    unsigned long int N = N_EM_cond;
    unsigned long int MonteCarloSimNum = 5000000;
    double running_sum=0;
    
    
    for(int i=0;i<MonteCarloSimNum;i++)
    {
        double tau_hat = GeneratorFirstHittingTime(unifMt, abs(V_BS_time_changed.BarrierFunction(x)));
        
        vector<double> running_x = x;
        bool the_process_hits_the_barrier;
        
        running_x = EulerMaruyamaSchemeWithConditionalBMAndStoppingCond(norm_rand_gen, N_EM_cond, V_BS_time_changed, j_star, tau_hat, T, running_x, &the_process_hits_the_barrier);
        
        if(the_process_hits_the_barrier)
        {
            double tau = running_x(V_BS_time_changed.GetBarrierMonitoringIndex());
            
            running_sum += BlackScholesTrueExpEurCall(running_x(0), mu, sigma, T-tau, K);
            
        }
        else
        {
            running_sum += 0;
        }

        //running_sum += IteratedRandomOperatorForEulerMaruyama(norm_rand_gen, N_EM, V_BS_normal, id, T, y);
        //IteratedRandomOperatorForEulerMaruyamaMerger(norm_rand_gen, unifMt, N_EM_cond, N_EM, V_BS_time_changed, V_BS_normal, j_star, T, f, x);
        //IteratedRandomOperatorForEulerMaruyamaMerger(norm_rand_gen, unifMt, N_EM_cond, N_EM, V_BS_time_changed, V_BS_normal, j_star, T, barrier_level, f, x);
        //double tau = GeneratorFirstHittingTime(unifMt, barrier_level-x(0));
        //running_sum += (tau <= T) ? 1 : 0;
        //vector<double> X = EulerMaruyamaSchemeWithConditionalBMAndStoppingCond(norm_rand_gen,N, V_time_changed_test, j_star, tau, T, barrier_level-x(0), x);
        //running_sum += X(V_time_changed_test.GetBarrierMonitoringIndex())> std::numeric_limits<double>::max()?0:f(X);
        //IteratedRandomOperatorForEulerMaruyamaMerger(norm_rand_gen, unifMt, N_EM_cond, N_EM, V_time_changed_test, V_BS_normal, j_star, T, barrier_level, f, x);
        //IteratedRandomOperatorForEulerMaruyama(norm_rand_gen, N_EM, V_time_changed_test, f, T, x);
        //IteratedRandomOperatorForEulerMaruyamaConditional(norm_rand_gen, N_EM, V_BS_time_changed, j_star, f, T, barrier_level, x);
        //IteratedRandomOperatorForEulerMaruyamaMerger(norm_rand_gen, unifMt, N_EM_cond, N_EM, V_BS_time_changed, V_BS_normal, j_star, T, barrier_level, f, x);
    }
    
    std::cout << running_sum/MonteCarloSimNum << "\n";
    //std::cout << BlackScholesTruePrice(x(0), mu, sigma, T, K, barrier_level) << "\n";
    
    std::cout << BlackScholesTrueExpEurCallUpAndIn(x(0), mu, sigma, T, K, barrier_level) << "\n";
    
    return 0;
}

//10  8.80066
//100 8.81281



/*
 //The following results show that EM scheme doesn't work for Black--Scholes.
 vector<double> x(2);
 
 x(0) = 100;
 x(1) = 0;
 
 vector<double> y(1);
 
 y(0) = 100;
 
 double mu = 0.02;
 double sigma = 0.2;
 double barrier_level = 110;
 double T = 1.0;
 double K = 100;
 
 VectorFieldsBlackScholes V_BS_normal(mu, sigma);
 boost::function<double(vector<double>)> id = [=](vector<double> x){return x(0) ;};
 unsigned long int MonteCarloSimNum = 5000000;
 double running_sum=0;
 
 for(int i=0;i<MonteCarloSimNum;i++)
 {
 running_sum += IteratedRandomOperatorForEulerMaruyama(norm_rand_gen, N_EM, V_BS_normal, id, T, y);
 }
 
 std::cout << running_sum/MonteCarloSimNum << "\n";
 
 unsigned long int N_EM = 100;
 //output is 102.02
 
 unsigned long int N_EM = 1000;
 //output is 102.014
*/

/* normal Euler--Maruyama approx.
 vector<double> x(2);
 
 x(0) = 100;
 x(1) = 0;
 
 vector<double> y(1);
 
 y(0) = 100;
 
 double mu = 0.02;
 double sigma = 0.2;
 double barrier_level = 110;
 double T = 1.0;
 double K = 100;
 
 boost::function<double(vector<double>)> f = [=](vector<double> x){return (K < x(0)) ? (x(0)-K) : 0 ;};
 
 unsigned int j_star = 1;
 
 unsigned long int N_EM_cond = 1000;
 unsigned long int N_EM = 1000;
 unsigned long int MonteCarloSimNum = 5000000;
 double running_sum=0;

approx.    9.09156
true value 8.91604
*/



/*
 vector<double> x(2);
 
 x(0) = 100;
 x(1) = 0;
 
 double mu = 0.02;
 double sigma = 0.2;
 double barrier_level = 110;
 double T = 1.0;
 double K = 100;
 unsigned int j_star = 1;
 boost::function<double(vector<double>)> f = [=](vector<double> x){return (K < x(0)) ? (x(0)-K) : 0 ;};

 unsigned long int N_EM_cond = 10;
 unsigned long int N_EM = 10;
 unsigned long int N = N_EM_cond;
 unsigned long int MonteCarloSimNum = 5000000;
 
 8.82656
 8.97456

 
 unsigned long int N_EM_cond = 100;
 unsigned long int N_EM = 100;
 unsigned long int N = N_EM_cond;
 unsigned long int MonteCarloSimNum = 5000000;
 
 approx 8.81195
 true   8.97456
 
 
 unsigned long int N_EM_cond = 1000;
 unsigned long int N_EM = 1000;
 unsigned long int N = N_EM_cond;
 unsigned long int MonteCarloSimNum = 5000000;
 
 approx 8.82543
 true   8.97456

 */

/*
 To use exponential r.v. or bernulli r.v.s you can paste the following.
 
 RandomSequenceGenerator<MersenneTwisterUniformRng> unifMtSeq(1,unifMt);
 boost::function<double(double)> ExpFunc = &ExpInverseCumul;
 boost::function<int(double)> Lambda1Func = &Lambda1GenFunc;
 boost::function<int(double)> Lambda2Func = &Lambda2GenFunc;
 InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<double(double)>> exp_rand_gen(unifMtSeq, ExpFunc);
 InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<int(double)>> Lambda1_rand_gen(unifMtSeq, Lambda1Func);
 InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<int(double)>> Lambda2_rand_gen(unifMtSeq, Lambda2Func);
*/

/*
 If you want an output of conditional_BM_paths, paste this:
 
 std::ofstream outputfile("/Users/asanoryo/Documents/XcodeProjects/time_changed_barrier_option/time_changed_barrier_option/cond_BM_paths.csv");
 for(int i=0; i< MonteCarloSimNum; i++)
 {
 std::vector<double> one_path = GeneratorConditionalBM(norm_rand_gen, N, t, y);
 for(int i = 0; i<N-1;i++)
 {
 outputfile << one_path[i] << ",";
 }
 outputfile << one_path[N] << "\n";
 }
 outputfile.close();
*/


/*
 To extract maximum values of paths, paste this:
 
 std::ofstream outputfile("/Users/asanoryo/Documents/XcodeProjects/time_changed_barrier_option/time_changed_barrier_option/max_and_argmax_list.csv");
 outputfile << "argmax, max\n";
 for(int i=0; i<N;i++)
 {
 std::vector<double> path = GeneratorConditionalBM(norm_rand_gen, N, t, y);
 std::vector<double>::iterator maxIt = std::max_element(path.begin(), --path.end());
 double max = *maxIt;
 size_t argmax = std::distance(path.begin(), maxIt);
 outputfile << argmax << ", " << max <<"\n";
 }
 outputfile.close();
*/


/*
 To test random number generators, paste this:
 
 std::ofstream BetaFile("/Users/asanoryo/Documents/python/Beta.csv");
 std::ofstream GammaFile("/Users/asanoryo/Documents/python/Gamma.csv");
 for(int i=0; i< MonteCarloSimNum; i++)
 {
 BetaFile << GeneratorBetaRandNum(unifMt) << ",";
 GammaFile << GeneratorGammaRandNum(unifMt) << ",";
 }
 BetaFile.close();
 GammaFile.close();
 
 std::ofstream HittingTimeFile("/Users/asanoryo/Documents/python/hitting_time.csv");
 for(int i=0; i< MonteCarloSimNum; i++)
 {
 HittingTimeFile << GeneratorFirstHittingTime(unifMt, 1.0) << ",";
 }
 HittingTimeFile.close();
*/

/*
 std::ofstream BetaFile("/Users/asanoryo/Documents/python/Beta.csv");
 std::ofstream GammaFile("/Users/asanoryo/Documents/python/Gamma.csv");
 std::ofstream ExpFile("/Users/asanoryo/Documents/python/Exp.csv");
 std::ofstream HittingTimeFile("/Users/asanoryo/Documents/python/hitting_time.csv");
 for(int i=0; i< MonteCarloSimNum; i++)
 {
 BetaFile << GeneratorBetaRandNum(unifMt) << "\n";
 GammaFile << GeneratorGammaRandNum(unifMt) << "\n";
 ExpFile << GeneratorExpRandNum(unifMt) << "\n";
 HittingTimeFile << GeneratorFirstHittingTime(unifMt, barrier_level-x(0)) << "\n";
 
 }
 BetaFile.close();
 GammaFile.close();
 ExpFile.close();
 HittingTimeFile.close();
 */

/*
 double t = 0;
 std::ofstream DensityFile("/Users/asanoryo/Documents/python/FirstHittingTimeDeinstyFunc.csv");
 for(int i = 0; i<MonteCarloSimNum;i++)
 {
 t += 1;
 DensityFile << t << "," << DensityFuncOfBMFirstHittingTime(barrier_level, t) << "\n";
 }
 DensityFile.close();
 */

/*tests the CDFOfGeometricBMFirstHittingTime, which should provide the explicit solution of P[tau <= T] for Balck Scholes.
 vector<double> x_BS_test(1);
 x_BS_test(0) = 1.0;
 for(int i=0;i<MonteCarloSimNum;i++)
 {
 bool the_process_hits_the_barrier;
 
 vector<double> X = EulerMaruyamaSchemeWithStoppingCond(norm_rand_gen, N, V_BS_normal, T, x_BS_test, barrier_level, &the_process_hits_the_barrier);
 running_sum += (the_process_hits_the_barrier)?1:0;
 }
 std::cout << running_sum/MonteCarloSimNum << "\n";
 std::cout << CDFOfGeometricBMFirstHittingTime(mu, sigma, T, barrier_level, x_BS_test(0)) << "\n";
 */

/* tests the accuracy of the EulerMaruyamaSchemeWIthConditionalBMAndStoppingCOnd with Black Scholes model in case the test function is binary, which results in the computation of the prob. P[tau <= T].
 for(int i=0;i<MonteCarloSimNum;i++)
 {
 double tau = GeneratorFirstHittingTime(unifMt, barrier_level-x(0));
 bool the_process_hits_the_barrier;
 vector<double> X = EulerMaruyamaSchemeWithConditionalBMAndStoppingCond(norm_rand_gen, N, V_BS_time_changed, j_star, tau, T, barrier_level-x(0), x, &the_process_hits_the_barrier);
 running_sum += (the_process_hits_the_barrier)?1:0;
 }
 
 std::cout << running_sum/MonteCarloSimNum << "\n";
 std::cout << CDFOfGeometricBMFirstHittingTime(mu, sigma, T, barrier_level, x(0)) << "\n";
 */

/*also tests the accuracy of E[1_{tau <= T} f(X(tau))] where X is BS and tau with upper barrier.
 for(int i=0;i<MonteCarloSimNum;i++)
 {
 double tau = GeneratorFirstHittingTime(unifMt, barrier_level-x(0));
 bool the_process_hits_the_barrier;
 vector<double> X = EulerMaruyamaSchemeWithConditionalBMAndStoppingCond(norm_rand_gen, N, V_BS_time_changed, j_star, tau, T, barrier_level-x(0), x, &the_process_hits_the_barrier);
 running_sum += (the_process_hits_the_barrier)?f(X):0;
 }
 
 std::cout << running_sum/MonteCarloSimNum << "\n";
 
 double running_sum2=0;
 vector<double> x_BS_test(1);
 x_BS_test(0) = 1.0;
 for(int i=0;i<MonteCarloSimNum;i++)
 {
 bool the_process_hits_the_barrier;
 
 vector<double> X = EulerMaruyamaSchemeWithStoppingCond(norm_rand_gen, N, V_BS_normal, T, x_BS_test, barrier_level, &the_process_hits_the_barrier);
 running_sum2 += (the_process_hits_the_barrier)?f(X):0;
 }
 
 std::cout << running_sum2/MonteCarloSimNum << "\n";
 */


