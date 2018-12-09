//
//  vector_field_class.hpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/19.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#pragma once
#include <boost/numeric/ublas/vector.hpp>

class VectorFields{
public:
    VectorFields(){};
    virtual boost::numeric::ublas::vector<double> GetVal(int direction, boost::numeric::ublas::vector<double> current_point) const=0;
    int GetNumOfVecFields() const;
    
protected:
    int num_of_vector_fields;   //the dimension of Brownian motion, i.e. d.
    int dim_of_state_space;     //the dimension of the state space which coincides with that of the underlying diffusion process X(t), i.e. n.
};


class VectorFieldsTimeChangedBlackScholesWithUpperBarrier : public VectorFields {
public:
    VectorFieldsTimeChangedBlackScholesWithUpperBarrier(double mu_, double sigma_, double barrier_);
    boost::numeric::ublas::vector<double> GetVal(int direction, boost::numeric::ublas::vector<double> current_point) const;     //returns V_i(x) if you put .GetVal(i,x)
    int GetNumOfVecFields() const;      //returns the dimension of Brownian motion.
    
private:
    double mu;
    double sigma;
    double barrier;
};


class VectorFieldsTest : public VectorFields {
public:
    VectorFieldsTest();
    boost::numeric::ublas::vector<double> GetVal(int direction, boost::numeric::ublas::vector<double> current_point) const;
    int GetNumOfVecFields() const;
};

/* vector_field_class_hpp */
