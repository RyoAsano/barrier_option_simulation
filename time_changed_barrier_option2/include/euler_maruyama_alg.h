#pragma once

namespace euler_maruyama {

void BlackScholesOneStep(double drift, double volatility, double std_gaussian_rand_num, double time_step, double* running_process);

void TimeChangedBlackScholesUpperBarrierOneStep(double volatility, double stochastic_integrator_increment, 
        double time_increment, double *asset_price_process_ptr, double *barrier_hitting_time_ptr);

void TimeChangedBlackScholesUpperBarrierSample(double volatility, double maturity, double sde_initial_value,
        double borwnian_bridge_goal_value, double brownian_bridge_goal_time, unsigned int num_of_subdivisions, unsigned long num_of_paths);

void TimeChangedBlackScholesUpperBarrierOnePath(double volatility, double maturity, double sde_initial_value,
        double barrier_level, double brownian_bridge_goal_time, unsigned int num_of_subdivisions, const double *norm_rand_num_array,
        double *asset_price_process_ptr, double *barrier_hitting_time_ptr, bool *the_process_hits_the_barrier_before_maturity_ptr);
}//namespace euler_maruyama