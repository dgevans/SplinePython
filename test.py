import numpy as np
from Spline.Spline import Spline

x = np.linspace(0,1,10)
X = Spline.makeGrid([x,x])
y = 1/(0.5 +(X[:,0]+X[:,1])**2)

fV = Spline(X,y)

xhat = [0.2,0.2]
print fV(xhat)

print fV(xhat,[1,0])
print fV(xhat,[0,1])