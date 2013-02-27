//
//  SplineEigen.cpp
//  SplineInterpolation
//
//  Created by David Evans on 9/26/12.
//
//

#include "SplineEigen.hpp"
#include <iostream>

#include <boost/python.hpp>


using namespace Eigen;
using std::endl;
Spline::Spline(long _N, const std::vector<int> &_k,const std::vector<breakpoints> &_v,const Eigen::VectorXd &_c):N(_N),k(_k),v(_v),c(_c)
{
}

Spline::Spline(MatrixXd X, VectorXd Y, std::vector<int> _k,bool sorted): k(_k),N(_k.size()),v(0)
{
    if(!sorted)
        sortXY(X, Y);
    bool complete = setUpBreakpoints(X);
    if(complete)
        fitcomplete(X, Y);
    else
        fitIncomplete(X, Y);
}

void Spline::sortXY(MatrixXd &X, VectorXd &Y)
{
    std::map<RowVectorXd, double, vecComp> temp;
    for (int i=0; i < X.rows(); i++)
        temp[X.row(i)] = Y(i);
    int i = 0;
    X.resize(temp.size(), X.cols());
    Y.resize(temp.size(),1);
    for(std::map<RowVectorXd, double, vecComp>::iterator it = temp.begin(); it!=temp.end(); it++)
    {
        X.row(i) = it->first;
        Y(i) = it->second;
        i++;
    }
}



bool Spline::setUpBreakpoints(const MatrixXd &X)
{
    std::vector<VectorXd> x;
    //for each i store the unique points
    int Nprod = 1;
    int NX = X.rows();
    for(int i =0; i < N; i++)
    {
        x.push_back(X.col(i));
        std::sort(&x[i](0), &x[i](NX-1)+1);
        double *it = std::unique(&x[i](0), &x[i](NX-1)+1);
        x[i].conservativeResize(it-&x[i](0));
        
        Nprod *= x[i].rows();
    }
    double frac;
    bool complete = false;
    if(Nprod != X.rows())
        frac = std::pow(X.rows()/double(Nprod),1.0/N);
    else
    {
        complete = true;
        frac = 1;
    }
    //for now assume model is fuly saturated
    for(int i=0; i < N; i++)
    {
        //to accout for missing values
        int n = frac*x[i].rows();
        //number of breakpoints
        int np = n-k[i]+1;
        VectorXd bp(np);
        bp(0) = x[i](0);
        bp(np-1) = x[i](x[i].rows()-1);
        for(int ip = 1; ip <np-1; ip++)
        {
            //proportion of point < ip
            double P = double(ip)/np;
            int i1 = floor(P*n);
            int i2 = ceil(P*n);
            bp(ip) =x[i](i1) + (P*n-i1)*(x[i](i2)-x[i](i1));
        }
        v.push_back(bp);
    }
    return complete;
}

void Spline::fitcomplete(const MatrixXd &X, const VectorXd &Y)
{
    
    std::vector<VectorXd> x;
    //for each i store the unique points
    int Nprod = 1;
    int NX = X.rows();
    for(int i =0; i < N; i++)
    {
        x.push_back(X.col(i));
        std::sort(&x[i](0), &x[i](NX-1)+1);
        double *it = std::unique(&x[i](0), &x[i](NX-1)+1);
        x[i].conservativeResize(it-&x[i](0));
        
        Nprod *= x[i].rows();
    }
    if(X.cols() < 1100){
        //vector holding the basis matrix for each set of dimensions
        std::vector<MatrixXd> Phi;
        for(int i = 0; i < N; i++)
        {
            int n = x[i].rows();
            Phi.push_back(MatrixXd(n,n));
            //iterate over x's
            for(int ix =0; ix < n; ix++)
            {
                //iterate over basis
                for(int ib=0;ib<n;ib++)
                {
                    Phi[i](ix,ib) = Basis(v[i], k[i], ib, x[i](ix));
                }
            }
        }
        //create inverse matrix
        MatrixXd B = Phi[0].inverse();
        for(int j = 1; j < N; j++)
            B= kron(Phi[j].inverse(),B);
        c = B*Y;
    }else
    {
        SparseMatrix<double,RowMajor> Phi(1,1);
        Phi.coeffRef(0,0) = 1;
        Phi.makeCompressed();
        for (int i =0; i < N; i++) {
            int n = x[i].rows();
            //iterate over x's
            std::vector<Triplet<double> > TripList;
            TripList.reserve(n*k[i]);
            for(int ix =0; ix < n; ix++)
            {
                //iterate over basis
                for(int ib = 0; ib < n; ib++)
                {
                    double temp = Basis(v[i], k[i], ib, x[i](ix));
                    if(temp != 0)
                        TripList.push_back(Triplet<double>(ix,ib,temp));
                }
            }
            SparseMatrix<double,RowMajor> B(n,n);
            B.setFromTriplets(TripList.begin(),TripList.end());
            Phi = kron(B,Phi);
        }
        BiCGSTAB<SparseMatrix<double,RowMajor> >solver;
        //UmfPackLU<SparseMatrix<double,RowMajor> > solver;
        solver.compute(Phi);
        c = solver.solve(Y);
    }
}

void Spline::fitIncomplete(const MatrixXd &X, const VectorXd &Y)
{
    //get number of total number of rows
    int ncols = 1;
    for(int i = 0; i < N; i++)
        ncols *= v[i].p()+k[i]-1;
    if(X.rows() < 1100)
    {
        MatrixXd Phi(X.rows(),ncols);
        //fill each row
        for(int i = 0; i < X.rows(); i++)
        {
            RowVectorXd B(1);
            B(0) = 1;
            for(int j = 0; j < N; j++)
            {
                int n = v[j].p()+k[j]-1;
                RowVectorXd temp(n);
                for(int ib = 0; ib < n; ib++)
                    temp(ib) = Basis(v[j],k[j],ib,X(i,j));
                B = kron(temp,B);
            }
            Phi.row(i) = B;
        }
        //c = Phi.jacobiSvd(ComputeThinU | ComputeThinV).solve(Y);
        MatrixXd A = (Phi.transpose())*Phi;
        VectorXd b = (Phi.transpose())*Y;
        c = A.ldlt().solve(b);
    }else
    {
        SparseMatrix<double,RowMajor> Phi(X.rows(),ncols);
        //fill each row
        for(int i = 0; i < X.rows(); i++)
        {
            SparseMatrix<double,RowMajor> B(1,1);
            B.coeffRef(0,0) = 1;
            B.makeCompressed();
            for(int j =0; j < N; j++)
            {
                int n = v[j].p()+k[j]-1;
                SparseMatrix<double,RowMajor> temp(1,n);
                temp.reserve(2*(k[j]+1));
                for(int ib = 0; ib < n; ib ++)
                {
                    double v_ib = Basis(v[j],k[j],ib,X(i,j));
                    if(v_ib != 0)
                        temp.insert(0, ib) = v_ib;
                }
                B = kron(temp,B);
            }
            Phi.row(i) = B;
        }
        SparseMatrix<double> A = (Phi.transpose())*Phi;
        VectorXd b = ( Phi.transpose() ) * Y;
        
        //CholmodDecomposition<SparseMatrix<double> > solver;
        BiCGSTAB<SparseMatrix<double> > solver;
        //Solve the matrix
        solver.compute(A);
        c = solver.solve(b);
    }
}

double Spline::operator()(const MatrixXd &X) const
{
    //setup
    VectorXd B(1);
    B(0) = 1;
    
    for (int i = 0; i < N; i++) {
        int n =  v[i].p()+k[i]-1;
        VectorXd temp(n);
        for (int ib = 0; ib < n; ib++) {
            temp(ib) = Basis(v[i], k[i], ib, X(i));
        }
        B= kron(temp,B);
    }
    
    return B.dot(c);
}

double Spline::operator()(const MatrixXd &X, const VectorXi &d) const
{
    //setup
    VectorXd B(1);
    B(0) = 1;
    
    for (int i = 0; i < N; i++) {
        int n =  v[i].p()+k[i]-1;
        VectorXd temp(n);
        for (int ib = 0; ib < n; ib++) {
            temp(ib) = Basis(v[i], k[i], ib, X(i),d(i));
        }
        B= kron(temp,B);
    }
    
    return B.dot(c);
}



bool vecComp::operator()(const RowVectorXd &x1, const RowVectorXd &x2) const
{
    for (int i = x1.cols()-1; i >=0; i--) {
        if(x1(i) < x2(i))
            return true;
        else if(x2(i) < x1(i))
            return false;
    }
    return false;
}



MatrixXd Spline::makeGrid(std::vector<VectorXd> &x)
{
    int N = 1;
    std::vector<int> n;
    n.push_back(N);
    for(int i = 0; i < x.size(); i++)
    {
        N *= x[i].rows();
        n.push_back(N);
    }
    MatrixXd X(N,x.size());
    for(int i = 0; i < N; i++)
    {
        int temp = i;
        for(int j = x.size()-1; j>=0; j--)
        {
            X(i,j) = x[j](temp/n[j]);
            temp %=n[j];
        }
    }
    return X;
    
}

double Spline::Basis(const breakpoints &v,int k, int j, double x, int der) const{
    if(der <=0 ){
        if (k <= 0)
        {
            return (x >= v[j])*(x<v[j+1])||(j==v.p()-2)*(x>=v[j])||(j==0)*(x < v[j+1]);//allows for extrapolation
        }
        else {
            //if j == n it will be too large for k-1 so drop second term
            if(j==v.p()+k-2)
                return (x-v[j-k])/(v[j]-v[j-k])*Basis(v,k-1,j-1,x);
            else if(j == 0)//drop first term which whould have j = -1
                return (v[j+1]-x)/(v[j+1]-v[j+1-k])*Basis(v,k-1,j,x);
            else
                return (x-v[j-k])/(v[j]-v[j-k])*Basis(v,k-1,j-1,x) + (v[j+1]-x)/(v[j+1]-v[j+1-k])*Basis(v,k-1,j,x);
        }
    }else{
        if (k <= 0) {
            if ( (x == v[j]) || (x==v[j+1])) {
                return std::numeric_limits<double>::quiet_NaN();
            }else
                return 0.0;
        }else
        {
            if(j==v.p()+k-2)
                return k/(v[j]-v[j-k])*Basis(v, k-1, j-1,x, der-1);
            else if(j==0)
                return -k/(v[j+1]-v[j+1-k])*Basis(v,k-1,j,x,der-1);
            else
                return k/(v[j]-v[j-k])*Basis(v, k-1, j-1,x, der-1)-k/(v[j+1]-v[j+1-k])*Basis(v,k-1,j,x,der-1);
        }
        
    }
}

void saveVector(const VectorXd & vec, std::ofstream &fout)
{
    fout<<vec.rows()<<std::endl;
    for(int i =0; i < vec.rows(); i++)
        fout<<vec(i)<<std::endl;
}

Eigen::VectorXd loadVector(std::ifstream &fin)
{
    int N;
    fin>>N;
    VectorXd vec(N);
    for(int i= 0; i < N; i ++)
        fin>>vec(i);
    return vec;
}

void Spline::save(std::ofstream &fout)
{
    fout<<N<<endl;
    for(int i = 0; i < N; i++)
        fout<<k[i]<<endl;
    for(int i =0; i < N; i++)
        saveVector(v[i].getv(),fout);
    saveVector(c,fout);
}

Spline Spline::load(std::ifstream &fin)
{
    long N;
    fin>>N;
    std::vector<int> k(N);
    for (int i = 0; i< N; i++) {
        fin>>k[i];
    }
    std::vector<breakpoints> v;
    for (int i =0; i < N; i++) {
        v.push_back(breakpoints(loadVector(fin)));
    }
    VectorXd c = loadVector(fin);
    return Spline(N,k,v,c);
    
}

namespace bp = boost::python;

typedef Matrix<double,Dynamic,Dynamic,RowMajor> MatrixXdNP;
typedef Matrix<long,Dynamic,1> VectorXl;

class SplinePython
{
    Spline f;
public:
    SplinePython(PyObject * pyX,PyObject * pyY, bp::list k, bool sorted);
    
    void operator()(PyObject *xpy, PyObject * dpy, PyObject *pyY);
    
    double feval(PyObject *xpy, PyObject *pyd);
    
    void getCoeff(PyObject *pyC);
    
    int getCoeffSize();
    
};

SplinePython::SplinePython(PyObject * pyX,PyObject * pyY, bp::list pyk, bool sorted)
{
    int Xdim = PyArray_NDIM(pyX);
    npy_intp * Xshape = PyArray_DIMS(pyX);
    int nX1 = Xshape[0],nX2;
    if(Xdim==1)
        nX2 = 1;
    else
        nX2 = Xshape[1];
    Map<MatrixXdNP> X((double *)PyArray_DATA(pyX),nX1,nX2);
    npy_intp * Yshape = PyArray_DIMS(pyY);
    Map<VectorXd> Y((double *)PyArray_DATA(pyY),Yshape[0]);
    
    int nk = bp::len(pyk);
    std::vector<int> k(nk);
    for(int i =0; i < nk; i++)
        k[i] = bp::extract<int>(pyk[i]);
    
    f = Spline(X, Y, k,sorted);
}

void SplinePython::operator()(PyObject *pyX, PyObject * pyd, PyObject *pyY)
{
    int Xdim = PyArray_NDIM(pyX);
    npy_intp * Xshape = PyArray_DIMS(pyX);
    int nX1,nX2;
    if(Xdim==1)
    {
        nX1 = Xshape[0];
        nX2 = 1;
    }
    else
    {
        nX1 = Xshape[0];
        nX2 = Xshape[1];
    }
    Map<VectorXl> d((long *)PyArray_DATA(pyd),f.getN());
    Map<MatrixXdNP> X((double *)PyArray_DATA(pyX),nX1,nX2);
    Map<VectorXd> Y((double *)PyArray_DATA(pyY),nX1);
    for(int i =0; i < nX1; i++)
        Y[i] = f(X.row(i),d.cast<int>());
}

double SplinePython::feval(PyObject *pyX,PyObject *pyd)
{
    Map<VectorXd> X((double *)PyArray_DATA(pyX),f.getN());
    Map<VectorXl> d((long *)PyArray_DATA(pyd),f.getN());
    std::cout<<d<<std::endl;
    return f(X,d.cast<int>());
}

void SplinePython::getCoeff(PyObject *pyC)
{
    Map<VectorXd> Y((double *)PyArray_DATA(pyC),f.getCoeff().rows());
    Y = f.getCoeff();
}

int SplinePython::getCoeffSize()
{
    return f.getCoeff().rows();
}


using namespace boost::python;

BOOST_PYTHON_MODULE(Spline_cpp)
{
    class_<SplinePython>("Spline_cpp",init<PyObject*,PyObject *,list,bool>())
    .def("__call__",&SplinePython::operator())
    .def("getCoeff",&SplinePython::getCoeff)
    .def("feval",&SplinePython::feval)
    .def("getCoeffSize",&SplinePython::getCoeffSize);
}



