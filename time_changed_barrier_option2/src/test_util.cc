#include "test_util.h"
#include <cassert>
#include "cumulative_dist_func.h"
#include "prob_density_func.h"

namespace test_util {
/*
 * GaussianUnusualMoment computes
 * \mathbb{E}[(\alpha*Z-\beta*y)^{n}],
 * where \alpha, \beta and y are rv_coeff, fixed_coeff and fixed_val, resp.
 */
double GaussianUnusualMoment(double rv_coeff, double fixed_coeff, double fixed_val, unsigned int moment_order){
   
    double result=0;
    switch(moment_order){
        case 1: {
                    result=-fixed_coeff*fixed_val;
                    break;
                }    
        case 2: {
                    result=rv_coeff*rv_coeff+fixed_coeff*fixed_coeff*fixed_val*fixed_val;
                    break;
                }
        case 3: {
                    result=-fixed_coeff*fixed_val*(3.0*rv_coeff*rv_coeff+fixed_coeff*fixed_coeff*fixed_val*fixed_val);
                    break;
                }
        default: {
                     assert(false&&"you have to put moment_order as 1,2 or 3.\n");
                 }
    } 
    return result;
}

/*
 * GaussianTruncatedMoment computes
 * \mathbb{E}[(\alpha*Z-\beta*y)^{n}1_{\alpha*Z-\beta*y>=0}],
 * where \alpha, \beta and y are rv_coeff, fixed_coeff and fixed_val.
 */
double GaussianTruncatedMoment(double rv_coeff, double fixed_coeff, double fixed_val, unsigned int moment_order){
    
    double result=0;
    switch(moment_order){
        case 0: {
                    result=1.0-cumulative_dist_func::StdGaussian(fixed_coeff/rv_coeff*fixed_val);
                    break;
                }
        case 1: {
                    result=rv_coeff*prob_density_func::StdGaussian(fixed_coeff/rv_coeff*fixed_val)
                            -fixed_val*fixed_coeff*(1.0-cumulative_dist_func::StdGaussian(fixed_coeff/rv_coeff*fixed_val));
                    break;
                }
        case 2: {
                   result=(rv_coeff*rv_coeff+fixed_coeff*fixed_coeff*fixed_val*fixed_val)*
                            (1.0-cumulative_dist_func::StdGaussian(fixed_coeff/rv_coeff*fixed_val))
                          -rv_coeff*fixed_coeff*fixed_val*prob_density_func::StdGaussian(fixed_coeff/rv_coeff*fixed_val);
                   break;
                }
        case 3: {
                    result=rv_coeff*(2.0*rv_coeff*rv_coeff+fixed_coeff*fixed_coeff*fixed_val*fixed_val)*
                                    prob_density_func::StdGaussian(fixed_coeff/rv_coeff*fixed_val)
                           -fixed_coeff*fixed_val*(3.0*rv_coeff*rv_coeff+fixed_coeff*fixed_coeff*fixed_val*fixed_val)*
                                    (1.0-cumulative_dist_func::StdGaussian(fixed_coeff/rv_coeff*fixed_val));
                    break;
                }
        default: {
                     assert(false&&"moment_order should be 0, 1, 2 and 3.\n");
                 }
    }

    return result;
}

/*
 * GaussianSignedMoment computes
 * \mathbb{E}[(\alpha*Z-\beta*y)^{n}*\sgn(\alpha*Z-\beta*y)],
 * where \alpha, \beta and y are rv_coeff, fixed_coeff and fixed_val.
 */
double GaussianSignedMoment(double rv_coeff, double fixed_coeff, double fixed_val, unsigned int moment_order){
    assert(moment_order==1||moment_order==2||moment_order==3);
    return -test_util::GaussianUnusualMoment(rv_coeff,fixed_coeff,fixed_val,moment_order)
            +2.0*test_util::GaussianTruncatedMoment(rv_coeff,fixed_coeff,fixed_val,moment_order);
}
}//namespace test_util