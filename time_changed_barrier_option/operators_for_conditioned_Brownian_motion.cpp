//
//  operators_for_conditioned_Brownian_motion.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/21.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "operators_for_conditioned_Brownian_motion.hpp"

using namespace QuantLib;

double lambda(double s, double t, double y, double x)
{
    double inside_absolute = x * sqrt((t-s)*s/t) - y*(t-s)/t;
    double abs_value = (inside_absolute>0) ? inside_absolute:-inside_absolute;
    
    return y-abs_value;
}

double IteratedRandomOperators(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                               const InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<double(double)>> &exp_rand_gen,
                               const InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<int(double)>> &Lambda1_rand_gen,
                               const InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<int(double)>> &Lambda2_rand_gen,
                               unsigned long int N, const VectorFields &V, unsigned int j_star,
                               const boost::function<double (vector<double>)> f, double t, double y, vector<double> x)
{
    double pi = M_PI;
    double running_factor = 1.0;            //running_factor denotes the iterated factors before Qg that are supposed to accumulate recursively.
    vector<double> running_x = x;           //running_x takes on the first argument of the operator Qf that is supposed to update recursively.
    double running_z = 0;                   //running_z takes on the second argument.
    unsigned int d = V.GetNumOfVecFields();
    std::ofstream lambda_history("/Users/asanoryo/Documents/XcodeProjects/time_changed_barrier_option/time_changed_barrier_option/lambda_history.csv", std::ios::app);
    
    lambda_history << running_z << ",";
    
    
    
    //impliment the recursive argument up to k <= N-1.
    for(unsigned long int i=0; i<=N-2; i++)
    {
        double s_k = (i+1)*t/N;                                  //Since i starts from 0, we regard (i+1) as k and i runs to N-2 so k = i+1 to N-1.
        double s_km = i*t/N;                                     //s_km denotes s_{k-1} (m stands for minus).
        
        double Z_j_star = norm_rand_gen.next().value;
        
        double arg_for_lambda = lambda(t/N, t-s_km, y-running_z, Z_j_star);
        
        running_factor *= 1.0;
        running_z += arg_for_lambda;
        
        lambda_history << running_z <<",";
        
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
    
    lambda_history << y << "\n";
    
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
    
    lambda_history.close();
    
    return running_factor * f(running_x);
}


/*
double IteratedRandomOperators(const VectorFields &V, double t, double y, unsigned long int N, unsigned int j_star,
                               const matrix<double> Z, const vector<double> E, const vector<double> Lambda1,
                               const vector<double> Lambda2, const boost::function<double (vector<double>)> f, vector<double> x)
{
    double pi = M_PI;
    double running_factor = 1.0;            //running_factor denotes the iterated factors before Qg that are supposed to accumulate recursively.
    vector<double> running_x = x;           //running_x takes on the first argument of the operator Qf that is supposed to update recursively.
    double running_z = 0;                   //running_z takes on the second argument.
    unsigned int d = V.GetNumOfVecFields();

    
    //impliment the recursive argument up to k <= N-1.
    for(unsigned long int i=0; i<=N-2; i++)
    {
        double s_k = (i+1)*t/N;                                  //Since i starts from 0, we regard (i+1) as k and i runs to N-2 so k = i+1 to N-1.
        double s_km = i*t/N;                                     //s_km denotes s_{k-1} (m stands for minus).
        
        double arg_for_lambda = lambda(t/N, t-s_km, y-running_z,
                                       Lambda1(i)*Z(j_star-1,i) + (1.0-Lambda1(i))*Lambda2(i)*sqrt(2.0*E(i)));
        
        running_factor *= 2*Lambda1(i) + 2.0/(y-running_z) * sqrt(2.0*(t-s_km)*t/(pi*(t-s_k)*N)) * (1.0-Lambda1(i)*Lambda2(i));
        running_z += arg_for_lambda;
        
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
                running_x += V.GetVal(j, running_x_before_update) * sqrt(t/N) * Z(j-1,i);
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
            running_x += V.GetVal(j, running_x_before_update) * sqrt(t/N) * Z(j-1,N-1);
        }
    }
    
    return running_factor * f(running_x);
}

 
 {
 double pi = M_PI;
 double running_factor = 1.0;            //running_factor denotes the iterated factors before Qg that are supposed to accumulate recursively.
 vector<double> running_x = x;           //running_x takes on the first argument of the operator Qf that is supposed to update recursively.
 double running_z = 0;                   //running_z takes on the second argument.
 unsigned int d = V.GetNumOfVecFields();
 
 
 //impliment the recursive argument up to k <= N-1.
 for(unsigned long int i=0; i<=N-2; i++)
 {
 double s_k = (i+1)*t/N;                                  //Since i starts from 0, we regard (i+1) as k and i runs to N-2 so k = i+1 to N-1.
 double s_km = i*t/N;                                     //s_km denotes s_{k-1} (m stands for minus).
 
 double Z_j_star = norm_rand_gen.next().value;
 double E = exp_rand_gen.nextSequence().value[0];
 int Lambda1 = Lambda1_rand_gen.nextSequence().value[0];
 int Lambda2 = Lambda2_rand_gen.nextSequence().value[0];
 
 double arg_for_lambda = lambda(t/N, t-s_km, y-running_z,
 Lambda1*Z_j_star + (1.0-Lambda1)*Lambda2*sqrt(2.0*E));
 
 running_factor *= 2*Lambda1 + 2.0/(y-running_z) * sqrt(2.0*(t-s_km)*t/(pi*(t-s_k)*N)) * (1.0-Lambda1)* Lambda2;
 running_z += arg_for_lambda;
 
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
 
 return running_factor * f(running_x);
 }

 
 */
