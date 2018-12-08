//
//  vector_field_class.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/19.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "vector_field_class.hpp"

int VectorFields::GetNumOfVecFields() const{
    return num_of_vector_fields;
}

VectorFieldsTimeChangedBlackScholesWithUpperBarrier::VectorFieldsTimeChangedBlackScholesWithUpperBarrier(double mu_, double sigma_, double barrier_)
{
    num_of_vector_fields = 2;
    dim_of_state_space = 3;
    mu = mu_;
    sigma = sigma_;
    barrier = barrier_;
}

vector<double> VectorFieldsTimeChangedBlackScholesWithUpperBarrier::GetVal(int direction, vector<double> current_point) const
{
    vector<double> result(dim_of_state_space);
    
    switch(direction){
        case 0:  //define the vector field of the drift term, i.e. V_0.
            result(0) = 0;
            result(1) = 0;
            result(2) = 1.0/(sigma*sigma*current_point[0]*current_point[0]);
            break;
        case 1: //define the vector field w.r.t. the first coordinate of Brownian motion, i.e. V_1.
            result(0) = 1.0;
            result(1) = mu * current_point[1] / (sigma*sigma*current_point[0]);
            result(2) = 0;
            break;
    }
    
    return result;
}

VectorFieldsTest::VectorFieldsTest()
{
    num_of_vector_fields = 2;
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

