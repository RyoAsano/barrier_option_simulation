#pragma once

namespace euler_maruyama {

void BlackScholesOneStep(double drift, double volatility, double std_gaussian_rand_num, double time_step, double* running_process);

void TimeChangedBlackScholesUpperBarrierOneStep(double drift, double volatility, double stochastic_integrator_increment, 
        double time_increment, double *asset_price_process_ptr, double *barrier_hitting_time_ptr, double *log_density_process_ptr);

void TimeChangedBlackScholesUpperBarrierOnePath(double drift, double volatility, double maturity, double sde_initial_value,
        double barrier_level, double brownian_bridge_goal_time, unsigned int num_of_subdivisions, const double *norm_rand_num_array,
        double *asset_price_process_ptr, double *barrier_hitting_time_ptr, double *log_density_process_ptr, 
        double *factor_ptr, bool *the_process_hits_the_barrier_before_maturity_ptr);

double TimeChangedBlackScholesMonteCarloUpAndInCall(double drift, double volatility, double maturity, double sde_initial_value,
        double strike, double barrier_level, unsigned int num_of_subdivisions, unsigned long long num_of_paths);

void TimeChangedBlackScholesMonteCarloUpAndInCallVariance(double drift, double volatility, double maturity, double sde_initial_value,
        double strike, double barrier_level, unsigned int num_of_subdivisions, unsigned long long num_of_paths, double *mean_ptr, double *variance_ptr);

}//namespace euler_maruyama
