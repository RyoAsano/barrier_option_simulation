#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "euler_maruyama_alg.h"
#include "statistical_functions.h"

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
        int c=getopt_long(argc, argv, "M:N:m::s::T::x::K::b::h", long_options, NULL);

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
                             << "This is the program for the computation of the up-and-in call vanilla option's variance."
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

    std::cout<<"The computation is starting under the following paramters:\n"
             <<"drift:"<<drift<<std::endl
             <<"volatility:"<<volatility<<std::endl
             <<"maturity:"<<maturity<<std::endl
             <<"initial value of the asset price:"<<sde_initial_value<<std::endl
             <<"strike:"<<strike<<std::endl
             <<"barrier level:"<<barrier_level<<std::endl
             <<std::endl
             <<"Also, we specified:\n"
             <<"number of paths:"<<num_of_paths<<std::endl
             <<"number of subdivisions:"<<num_of_subdivisions<<std::endl
             <<std::endl
             <<"the variance is now  under estimation..."<<std::endl
             <<std::endl;

    double estimated_mean=0;
    double estimated_variance=0;
    time_t start=time(0);
    euler_maruyama::TimeChangedBlackScholesMonteCarloUpAndInCallVariance(drift,volatility,
                                maturity,sde_initial_value,strike,barrier_level,num_of_subdivisions,num_of_paths,&estimated_mean,&estimated_variance);
    time_t end=time(0);

    int time_diff=end-start;
    int hours=(int)(time_diff/3600);
    int minutes=(int)((time_diff-hours*3600)/60);
    int seconds=time_diff-hours*3600-minutes*60;

    std::cout<<"Computation is done.\n"
             <<"The estimated mean is:"<<estimated_mean<<std::endl
             <<"The estimated variance is:"<<estimated_variance<<std::endl
             <<"The computation time is:"
             <<hours<<" h. "
             <<minutes<<" m. "
             <<seconds<<"s.\n";

    char ans_conti;
    while(1){
        std::cout<<"Would you also like to deduce the sample size?(Y/N)\n";
        std::cin >> ans_conti;
        if(ans_conti!='Y' && ans_conti!='N'){
            std::cout<<"Please answer Y or N.\n";
        }else{
            break;
        }
    }

    if(ans_conti=='Y'){
        std::cout << "Henceforth we compute the neccessary sample path size to achieve the following inequality:\n"
                  << "P[|MC-m|<e]>=p\n"
                  << "where MC is the Monte--Carlo estimator m is the target mean, e is an accuracy and p is an acceptance level.\n";
        double accuracy;
        statistical_functions::AcceptanceLevel acceptance_level;
        std::cout<<"Accuracy(e)?\n";
        std::cin >> accuracy;
        while(1){
            char ans_acceptance_level; 
            std::cout<<"Acceptance level?(1:p=95% 2:p=99% 3:p=99.9%)\n";
            std::cin >> ans_acceptance_level;
            if(ans_acceptance_level=='1'){
                acceptance_level=statistical_functions::NinetyFive;
                break;
            }else if(ans_acceptance_level=='2'){
                acceptance_level=statistical_functions::OneOverAHundred;
                break;
            }else if(ans_acceptance_level=='3'){
                acceptance_level=statistical_functions::OneOverAThousand;
                break;
            }else{
                std::cout << "Please answer 1, 2 or 3.\n";
            }
        }
        std::cout<<"The required sample size is:"
                 <<statistical_functions::SampleSizeGenerator(accuracy,estimated_variance,acceptance_level)
                 <<std::endl;
    }

    return 0;
}
