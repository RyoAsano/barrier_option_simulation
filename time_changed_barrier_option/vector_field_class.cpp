//
//  vector_field_class.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/19.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "vector_field_class.hpp"

using namespace boost::numeric::ublas;

int VectorFields::GetDimOfDiffusionCoeff() const
{
    return dim_of_diffusion_coeff;
}

int VectorFields::GetDimOfStateSpace() const
{
    return dim_of_state_space;
}

int VectorFieldsTimeChanged::GetBarrierMonitoringIndex() const
{
    return barrier_monitoring_index;
}

VectorFieldsTimeChangedBlackScholesWithUpperBarrier::VectorFieldsTimeChangedBlackScholesWithUpperBarrier(double mu_, double sigma_, double barrier_)
{
    dim_of_diffusion_coeff = 1;
    dim_of_state_space = 2;
    barrier_monitoring_index = 1;
    mu = mu_;
    sigma = sigma_;
    barrier = barrier_;
}


//First coordinate is the price process, second the density process and third the barrier_monitoring process.
vector<double> VectorFieldsTimeChangedBlackScholesWithUpperBarrier::GetVal(int direction, vector<double> current_point) const
{
    vector<double> result(dim_of_state_space);
    
    switch(direction){
        case 0:  //define the vector field of the drift term, i.e. V_0.
            result(0) = 0;
            result(1) = 1.0/(sigma*sigma*current_point[0]*current_point[0]);
            break;
        case 1: //define the vector field w.r.t. the first coordinate of Brownian motion, i.e. V_1.
            result(0) = 1.0;
            result(1) = 0;
            break;
    }
    
    return result;

}

double VectorFieldsTimeChangedBlackScholesWithUpperBarrier::BarrierFunction(boost::numeric::ublas::vector<double> spot) const
{
    return spot(0) - barrier;
}

VectorFieldsTimeChangedTest::VectorFieldsTimeChangedTest()
{
    dim_of_diffusion_coeff = 1;
    dim_of_state_space = 2;
    barrier_monitoring_index = 1;
    
}

vector<double> VectorFieldsTimeChangedTest::GetVal(int direction, boost::numeric::ublas::vector<double> current_point) const
{
    vector<double> result(dim_of_state_space);
    
    switch(direction){
        case 0:  //define the vector field of the drift term, i.e. V_0.
            result(0) = 0;
            result(1) = 1.0;
            break;
        case 1: //define the vector field w.r.t. the first coordinate of Brownian motion, i.e. V_1.
            result(0) = 1.0;
            result(1) = 0;
            break;
    }
    
    return result;
}

double VectorFieldsTimeChangedTest::BarrierFunction(boost::numeric::ublas::vector<double> spot) const
{
    return 1.0;
}

VectorFieldsBlackScholes::VectorFieldsBlackScholes(double mu_, double sigma_)
{
    dim_of_diffusion_coeff = 1;
    dim_of_state_space = 1;
    mu = mu_;
    sigma = sigma_;
}

vector<double> VectorFieldsBlackScholes::GetVal(int direction, boost::numeric::ublas::vector<double> current_point) const
{
    vector<double> result(dim_of_state_space);
    
    switch(direction){
        case 0:  //define the vector field of the drift term, i.e. V_0.
            result(0) = mu * current_point(0);
            break;
        case 1: //define the vector field w.r.t. the first coordinate of Brownian motion, i.e. V_1.
            result(0) = sigma * current_point(0);
            break;
    }
    
    return result;
}

VectorFieldsTest::VectorFieldsTest()
{
    dim_of_diffusion_coeff = 1;
    dim_of_state_space = 1;
}

vector<double> VectorFieldsTest::GetVal(int direction, vector<double> current_point) const
{
    vector<double> result(dim_of_state_space);
    
    switch(direction)
    {
        case 0:
            result(0) = 0;
            break;
        case 1:
            result(0) = 1.0;
            break;
    }
    
    return result;
}

