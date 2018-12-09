//
//  generator_Euler_Maruyama_conditional.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/08.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "generator_Euler_Maruyama_conditional.hpp"

using namespace QuantLib;
using namespace boost::numeric::ublas;

double lambda(double s, double t, double y, double x)
{
    double inside_absolute = x * sqrt((t-s)*s/t) - y*(t-s)/t;
    double abs_value = (inside_absolute>0) ? inside_absolute:-inside_absolute;
    
    return y-abs_value;
}

vector<double> EulerMaruyamaSchemeWithConditionalBM(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                                                    unsigned long int N, const VectorFields &V, unsigned int j_star, double t, double y, vector<double> x)
{
    vector<double> running_x = x;           //running_x takes on the first argument of the operator Qf that is supposed to update recursively.
    double running_z = 0;                   //running_z takes on the second argument.
    unsigned int d = V.GetNumOfVecFields();
    
    //impliment the recursive argument up to k <= N-1.
    for(unsigned long int i=0; i<=N-2; i++)
    {
        //Since i starts from 0, we regard (i+1) as k and i runs to N-2 so k = i+1 to N-1.
        double s_km = i*t/N;                                     //s_km denotes s_{k-1} (m stands for minus).
        
        double Z_j_star = norm_rand_gen.next().value;
        double arg_for_lambda = lambda(t/N, t-s_km, y-running_z, Z_j_star);
        
        running_z += arg_for_lambda;                            //update conditional brownian motion.
        
        vector<double> running_x_before_update = running_x;
        running_x = running_x_before_update;
        for(unsigned int j=0; j <=d; j++)
        {
            if(j==0)
            {
                running_x += V.GetVal(0, running_x_before_update) * (t/N);
            }
            else if(j==j_star)
            {
                running_x += V.GetVal(j_star, running_x_before_update) * arg_for_lambda;
            }
            else
            {
                double Z = norm_rand_gen.next().value;
                running_x += V.GetVal(j, running_x_before_update) * sqrt(t/N) * Z;
            }
        }
    }
    
    
    //the last substitution for k=N.
    vector<double> running_x_before_update = running_x;
    running_x = running_x_before_update;
    for(unsigned int j=0; j <=d; j++)
    {
        if(j==0)
        {
            running_x += V.GetVal(0, running_x_before_update) * (t/N);
        }
        else if(j==j_star)
        {
            running_x += V.GetVal(j_star, running_x_before_update) * (y-running_z);
        }
        else
        {
            double Z = norm_rand_gen.next().value;
            running_x += V.GetVal(j, running_x_before_update) * sqrt(t/N) * Z;
        }
    }
    
    return running_x;
    
}


double IteratedRandomOperatorForEulerMaruyamaConditional(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                                                         unsigned long int N, const VectorFields &V, unsigned int j_star,
                                                         const boost::function<double (vector<double>)> f, double t, double y, vector<double> x)
{
    vector<double> running_x = EulerMaruyamaSchemeWithConditionalBM(norm_rand_gen,N,V,j_star,t,y,x);
    return f(running_x);
}


std::vector<double> GeneratorConditionalBM(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen, unsigned long int N, double t, double y)
{
    std::vector<double> result;
    
    double running_z = 0;                   //running_z takes on the second argument.
    result.push_back(running_z);
    
    for(unsigned long int i=0; i<=N-2; i++)
    {
        double s_km = i*t/N;                                     //s_km denotes s_{k-1} (m stands for minus).
        double Z_j_star = norm_rand_gen.next().value;
        double arg_for_lambda = lambda(t/N, t-s_km, y-running_z, Z_j_star);
        
        running_z += arg_for_lambda;                            //update conditional brownian motion.
        result.push_back(running_z);
    }
    
    running_z = y;
    result.push_back(running_z);
    
    return result;
    
}
