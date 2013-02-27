//
//  main.cpp
//  SplinePython
//
//  Created by David Evans on 2/25/13.
//  Copyright (c) 2013 David Evans. All rights reserved.
//

#include <iostream>
#include <boost/python.hpp>
#include <eigen3/Eigen/Core>
#include "SplineEigen.hpp"

using namespace Eigen;

int main(int argc, const char * argv[])
{
    
    VectorXd x = VectorXd::LinSpaced(100, 0, 3);
    //std::cout<<x<<std::endl;
    VectorXd y(x.rows());
    for (int i=0 ; i< y.rows(); i++) {
        y(i) = exp(x(i));
    }
    std::cout<<y<<std::endl;
    Spline f(x, y, {3});
    
    VectorXd t = VectorXd::LinSpaced( 1000, 0, 3);
    for (int i= 0; i < t.rows(); i++) {
        std::cout<<f(t.row(i))<<","<<exp(t(i))<<std::endl;
    }
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}

