//
//  function.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/26.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "function.hpp"
#include <boost/function.hpp>
#include <cmath>

double ExpInverseCumul(double unif)
{
    return -log(unif);
}

int Lambda1GenFunc(double unif)
{
    if(unif <= 0.5)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int Lambda2GenFunc(double unif)
{
    if(unif <= 0.5)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

