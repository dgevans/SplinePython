import numpy as np
from Spline import Spline
import matplotlib.pyplot as plt

x = np.linspace(0,1,10)
y = np.exp(x)
f = Spline(x,y)
g = Spline(x,x)
X = Spline.makeGrid([x,x])
def g(x,y):
    return 1/(0.5 +(x+y)**2)
y = 1/(0.5 +(X[:,0]+X[:,1])**2)

fV = Spline(X,y)

xhat = 0.5*np.ones((1000,2))
xhat[:,0] = np.linspace(0,1,1000)

plt.plot(xhat[:,0],map(g,xhat[:,0],xhat[:,1]),'b')
plt.plot(xhat[:,0],fV(xhat),'g')
plt.show()
#print fV(xhat)

#print fV(xhat,[1,0])
#print fV(xhat,[0,1])