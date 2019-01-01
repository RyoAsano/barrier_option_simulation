//
//  generator_Euler_Maruyama_conditional.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/08.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "generator_Euler_Maruyama_conditional.hpp"
#include <iostream>
#include <fstream>

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
    unsigned int d = V.GetDimOfDiffusionCoeff();
    
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
                double a = running_x(0);
                double b = running_x(1);
                a;
            }
            else if(j==j_star)
            {
                running_x += V.GetVal(j_star, running_x_before_update) * arg_for_lambda;
                double a = running_x(0);
                double b = running_x(1);
                a;
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


vector<double> EulerMaruyamaSchemeWithConditionalBMAndStoppingCond(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                                                                   unsigned long int N, const VectorFieldsTimeChanged &V, unsigned int j_star, double t, double T,
                                                                   vector<double> x, bool *the_process_hits_the_barrier)
{
    vector<double> running_x = x;           //running_x takes on the first argument of the operator Qf that is supposed to update recursively.
    double running_z = 0;                   //running_z takes on the second argument.
    unsigned int d = V.GetDimOfDiffusionCoeff();
    int i_barrier = V.GetBarrierMonitoringIndex();
    double y =  abs(V.BarrierFunction(x));
    *the_process_hits_the_barrier = true;
    
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
                double a = running_x(0);
                double b = running_x(1);
                a;
            }
            else if(j==j_star)
            {
                running_x += V.GetVal(j_star, running_x_before_update) * arg_for_lambda;
                double a = running_x(0);
                double b = running_x(1);
                a;
            }
            else
            {
                double Z = norm_rand_gen.next().value;
                running_x += V.GetVal(j, running_x_before_update) * sqrt(t/N) * Z;
            }            
        }
        
        if(running_x(i_barrier) > T)
        {
            *the_process_hits_the_barrier = false;
            break;
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
    
    if(running_x(i_barrier) > T)
    {
        *the_process_hits_the_barrier = false;
    }

    
    return running_x;
    
}


double IteratedRandomOperatorEulerMaruyamaSchemeWithConditionalBMAndStoppingCondModified(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                                                                                         unsigned long int multiplier, const VectorFieldsTimeChanged &V, unsigned int j_star, double t, double T,
                                                                                         vector<double> x, const boost::function<double (double, vector<double>)> f)
{
    int d = V.GetDimOfDiffusionCoeff();
    int barrier_index = V.GetBarrierMonitoringIndex();
    double goal = abs(V.BarrierFunction(x));
    int N = int(std::ceil(multiplier*t));
    double aaa = t/N;
    
    vector<double> running_process = x;
    double running_condBM = 0;
    double running_factor = 1.0;
    bool the_process_hits_the_barrier = true;
    
    for(int k=0; k<N; k++){
        vector<double> spot = running_process;
        for(int j=0;j<=d;j++){
            if(j==0){
                running_process += V.GetVal(j, spot)*t/N;
            }else if(j!=j_star){
                double Z = norm_rand_gen.next().value;
                running_process += V.GetVal(j, spot)*Z*sqrt(t/N);
            }else{
                double s_k = (k+1)*t/N;
                double s_km = k*t/N;
                double Z = norm_rand_gen.next().value;
                bool need_to_jump_to_the_goal = (k==N-1);
                double del_condBM = need_to_jump_to_the_goal?
                (goal-running_condBM):
                (goal-running_condBM-abs(Z*sqrt((t-s_k)/(t-s_km) * t/N) - (goal-running_condBM)*(t-s_k)/(t-s_km)));
                running_process += V.GetVal(j, spot)* del_condBM;
                running_factor *= need_to_jump_to_the_goal?1.0:(1.0 + Z/(goal-running_condBM) * sqrt((t-s_km)/(t-s_k) * t/N));
                running_condBM += del_condBM;
            }
        }
        
        /*
        double mu = 0.02;
        double sigma = 0.2;
        double barrier_level = 110;
        double T = 1.0;
        double K = 100;

        double tes = x(0) + running_condBM;
        double tes2 = running_process(0);
        double test = spot(1) + 1.0/(sigma*sigma*spot(0)*spot(0)) * t/N;
        double test2 = running_process(1);
        */
        
        
        if(running_process(barrier_index)>T){
            /*
            std::ofstream output("/Users/asanoryo/Documents/python/col_of_tau_when_process_hits_the_barrier.csv", std::ios::app);
            output << t << "\n";
            output.close();
             */
            the_process_hits_the_barrier=false;
            break;
        }
    }
    
    double result=0;
    if(the_process_hits_the_barrier){
        result = running_factor * f(running_process(barrier_index), running_process);
    }
        
    return result;
}
