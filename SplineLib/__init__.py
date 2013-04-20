__author__ = 'dgevans'
from Spline_cpp import Spline_cpp
import numpy as np


class Spline(object):

    def __init__(self,X,y,k=None):
        self.f = Spline_cpp()
        self.fit(X,y,k)

    @staticmethod
    def makeGrid(x):
        N = 1
        n = []
        n.append(N)
        for i in range(0,len(x)):
            N *= x[i].shape[0]
            n.append(N)
        X = np.zeros((N,len(x)))
        for i in range(0,N):
            temp = i
            for j in range(len(x)-1,-1,-1):
                X[i,j] = x[j][temp/n[j]]
                temp %= n[j]
        return X

    def fit(self,X,y,k=None,sorted=False):
        """


        :param sorted:
        :param X:
        :param y:
        :param k:
        """
        X = np.asarray(X)
        y = np.asarray(y)
        if not X.flags['C_CONTIGUOUS']:
            X = np.ascontiguousarray(X)
        if not y.flags['C_CONTIGUOUS']:
            y = np.ascontiguousarray(y)
        assert y.ndim == 1
        if X.ndim > 1:
            assert X.ndim==2 and (X.shape[0] == y.shape[0] or X.shape[1] == y.shape[0])

            if(X.shape[0] != y.shape[0]):
                X = X.transpose()
            self.N = X.shape[1]
        else:
            assert X.shape[0] == y.shape[0]
            self.N = 1
        if k == None:
            k = []
            for i in range(0,self.N):
                k.append(3)
        k = list(k)
        assert len(k) == self.N

        self.f.fit(X,y,k,sorted)


    def __call__(self,X,d = None):
        if d == None:
            d = np.zeros(self.N,dtype=np.int)
        d = np.asarray(d,dtype=np.int)
        X = np.asarray(X)
        if not X.flags['C_CONTIGUOUS']:
            X = np.ascontiguousarray(X)
        if not d.flags['C_CONTIGUOUS']:
            d = np.ascontiguousarray(d)
            #make sure arrays are of the right shape
        X = np.atleast_2d(X)
        d = np.atleast_1d(d)
        assert X.ndim == 2 and (X.shape[0] == self.N or X.shape[1] == self.N) #X cannot be more than 2d, one dimension must be length N
        assert d.ndim == 1 and len(d) == self.N
        if X.shape[0] == self.N:
            X = X.transpose()
        y = np.zeros(X.shape[0])
        self.f(X,d,y)
        return y

    def getCoeffs(self):
        c = np.zeros(self.f.getCoeffSize())
        self.f.getCoeff(c)
        return c

