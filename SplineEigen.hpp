//
//  Spline.h
//  SplineInterpolation
//
//  Created by David Evans on 3/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef SplineInterpolation_SplineEigen_hpp
#define SplineInterpolation_SplineEigen_hpp
#include <vector>
#include <map>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
//#include <eigen3/Eigen/UmfPackSupport>
//#include <eigen3/Eigen/CholmodSupport>
#include <fstream>
#include <Python.h>
#include <numpy/arrayobject.h>


class breakpoints {

    //holds the breakpoints
    Eigen::VectorXd v;
    int np;
    friend class SplinePython;
public:
    breakpoints(const Eigen::VectorXd& v_):v(v_){np=v.rows();};
    breakpoints():v(),np(0){};
    
    double operator[](int j) const;
    
    double p() const {return np;} ;
    
    const Eigen::VectorXd& getv(){return v;};
};

class Spline {

    
    //order of the spline fit
    std::vector<int> k;
    //Dimension of the spline fit
    long N;
    
    //vector holding the breakpoints for each set of breakpoints
    std::vector<breakpoints> v;
    
    //vector holding the coefficients for each set of diminsions
    Eigen::VectorXd  c;
    
    
    bool setUpBreakpoints(const Eigen::MatrixXd &X);
    
    void fitcomplete(const Eigen::MatrixXd &X,const Eigen::VectorXd &Y);
    
    void fitIncomplete(const Eigen::MatrixXd &X, const Eigen::VectorXd &Y);
    
    void sortXY(Eigen::MatrixXd &X, Eigen::VectorXd &Y);
    
    double Basis(const breakpoints &v,int k, int j, double x, int der =0) const;
    
    Spline(long _N, const std::vector<int> &_k,const std::vector<breakpoints> &_v,const Eigen::VectorXd &_c);
    
    friend struct SplinePythonPickle;
    
public:
    
    Spline(){};
    
    Spline(Eigen::MatrixXd X, Eigen::VectorXd Y, std::vector<int> k,bool sorted = false);
    
    double operator()(const Eigen::MatrixXd &X) const;
    
    double operator()(const Eigen::MatrixXd &X, const Eigen::VectorXi &d) const;

    
    static Eigen::MatrixXd makeGrid(std::vector<Eigen::VectorXd> &x);
    
    static Spline load(std::ifstream &fin);
    
    const Eigen::VectorXd& getCoeff() const {return c;};
    
    void save(std::ofstream &fout);
    
    const std::vector<breakpoints>& getBreakpoints() const {return v;};
    
    int getN()const {return N;};
    
    const std::vector<int>& getk() const{return k;};
    
};



void saveVector(const Eigen::VectorXd & vec, std::ofstream &fout);
Eigen::VectorXd loadVector(std::ifstream &fin);

    
inline Eigen::MatrixXd kron(const Eigen::MatrixXd &X1, const Eigen::MatrixXd &X2)
    {
        long n1 = X1.rows(), n2 = X2.rows(), m1 = X1.cols(), m2 = X2.cols();
        int n = n1*n2;
        int m = m1*m2;
        
        Eigen::MatrixXd X(n,m);
        
        for(int i1 = 0; i1 < n1; i1++)
        {
            for(int i2 = 0; i2< n2; i2++)
            {
                int i = i1*n2+i2;
                for(int j1 = 0; j1 < m1; j1++)
                {
                    for(int j2 =0; j2<m2; j2++)
                    {
                        int j = j1*m2+j2;
                        X(i,j) = X1(i1,j1)*X2(i2,j2);
                    }
                }
            }
        }
        return X;
    }
    
inline Eigen::SparseMatrix<double,Eigen::RowMajor> kron(Eigen::SparseMatrix<double,Eigen::RowMajor> &X1, Eigen::SparseMatrix<double,Eigen::RowMajor> &X2)
    {
        int n1 = X1.rows(), n2 = X2.rows(), m1 = X1.cols(), m2 = X2.cols();
        int n = n1*n2;
        int m = m1*m2;
        
        std::vector<Eigen::Triplet<double> > list;
        list.reserve(X1.nonZeros()*X2.nonZeros());
        for (int i1 = 0; i1 < n1; i1++) 
        {
            for(int i2 = 0; i2 < n2; i2++)
            {
                int i = i1*n2+i2;
                for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it1(X1,i1); it1; ++it1)
                {
                    for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it2(X2,i2); it2; ++it2)
                    {
                        int j = it1.col()*m2+it2.col();
                        double X_ij = it1.value()*it2.value();
                        list.push_back(Eigen::Triplet<double>(i,j,X_ij));
                    }
                }
            }
        }
        
        Eigen::SparseMatrix<double,Eigen::RowMajor> X(n,m);
        X.setFromTriplets(list.begin(),list.end());
        return X;
    }
    
    

class vecComp {
public:
    bool operator()(const Eigen::RowVectorXd &x1, const Eigen::RowVectorXd &x2) const;
};



inline double breakpoints::operator[](int j) const
{
    if(j >= np)
        return v[np-1];
    if(j < 0)
        return v[0];
    return v[j];
}




#endif
