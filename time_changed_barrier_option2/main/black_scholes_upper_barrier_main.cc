#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "euler_maruyama_alg.h"
#include "black_scholes_explicit_expectation.h"

int main(int argc, char *argv[]){
    double drift=0;
    double volatility=0.2;
    double maturity=1.0;
    double sde_initial_value=100;
    double strike=100;
    double barrier_level=120;
    unsigned int num_of_subdivisions=0;
    unsigned long num_of_paths=0;

    struct option long_options[]={
        {"num_of_paths", required_argument, 0, 'M'},
        {"num_of_subdivisions", required_argument, 0, 'N'},
        {"drift", optional_argument, 0, 'm'},
        {"volatility", optional_argument, 0, 's'},
        {"maturity", optional_argument, 0, 'T'},
        {"initial_value", optional_argument, 0, 'x'},
        {"strike", optional_argument, 0, 'K'},
        {"barrier_level", optional_argument, 0, 'b'},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    while(1){
        int c=getopt_long(argc, argv, "M:N::m::s::T::x::K::b::h", long_options, NULL);

        if(c==-1){
            break;
        }

        switch(c){
            case 'M':{
                         num_of_paths=atoi(optarg);
                         break;
                     }
            case 'N':{
                         num_of_subdivisions=atoi(optarg);
                         break;
                     }
            case 'm':{
                         drift=(optarg)?atof(optarg):drift;
                         break;
                     }
            case 's':{
                         volatility=(optarg)?atof(optarg):volatility;
                         break;
                     }
            case 'T':{
                         maturity=(optarg)?atof(optarg):maturity;
                         break;
                     }
            case 'x':{
                         sde_initial_value=(optarg)?atof(optarg):sde_initial_value;
                         break;
                     }
            case 'K':{
                         strike=(optarg)?atof(optarg):strike;
                         break;
                     }
            case 'b':{
                         barrier_level=(optarg)?atof(optarg):barrier_level;
                         break;
                     }
            case 'h':{
                         std::cout 
                             << "This is the program for the computation of the up-and-in call vanilla option expectation."
                             << "The following options are available:\n"
                             << "-M or --num_of_paths: the number of paths for Monte-Carlo simulations.\n" 
                             << "-N or --num_of_subdivisions: the number of subdivisions of the time interval for the Euler--Maruyama Scheme.\n"
                             << "-m or --drift: the drift of the Black--Scholes SDE.\n"
                             << "-s or --volatility: the volatility of the Black--Scholes SDE.\n"
                             << "-T or --maturity: the maturity of the option.\n"
                             << "-x or --initial_value: the initial value of the underlying asset price process.\n"
                             << "-K or --strike: the strike of the option.\n"
                             << "-b or --barrier_level: the barrier level of the option.\n";
                         exit(EXIT_SUCCESS);
                         break;
                     }
            default:{
                        std::cerr << "undefined option " << c << "is specified\n";
                        exit(EXIT_FAILURE);
                        break;
                    }
        }
    }

    if(optind<argc){
        std::cerr<<"non-option arguments are reveived.\n";
        exit(EXIT_FAILURE);
    }

    double true_value=black_scholes_explicit_expectation::UpAndInCall(drift,volatility,maturity,strike,sde_initial_value,barrier_level);

    std::cout<<"The computation is starting under the following paramters:\n"
             <<"drift:"<<drift<<std::endl
             <<"volatility:"<<volatility<<std::endl
             <<"maturity:"<<maturity<<std::endl
             <<"initial value of the asset price:"<<sde_initial_value<<std::endl
             <<"strike:"<<strike<<std::endl
             <<"barrier level:"<<barrier_level<<std::endl
             <<"In this case, the true price is given by "<<true_value<<".\n"
             <<std::endl
             <<"Also, we specified:\n"
             <<"number of paths:"<<num_of_paths<<std::endl
             <<"number of subdivisions:"<<num_of_subdivisions<<std::endl
             <<std::endl
             <<"the approximating value is now under computation..."<<std::endl
             <<std::endl;

    time_t start=time(0);
    double approximation_value=euler_maruyama::TimeChangedBlackScholesMonteCarloUpAndInCall(drift, volatility, maturity, 
                                                            sde_initial_value, strike, barrier_level, num_of_subdivisions,num_of_paths);
    time_t end=time(0);

    int time_diff=end-start;
    int hours=(int)(time_diff/3600);
    int minutes=(int)((time_diff-hours*3600)/60);
    int seconds=time_diff-hours*3600-minutes*60;

    std::cout<<"Computation is done.\n"
             <<"The approximating value is:"<<approximation_value<<std::endl
             <<"The true value is:"<<true_value<<std::endl
             <<"The absolute difference is:"<<abs(approximation_value-true_value)<<std::endl
             <<"The computation time is:"
             <<hours<<" h. "
             <<minutes<<" m. "
             <<seconds<<"s.\n";

    return 0;
}
