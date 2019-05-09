import numpy as np
import matplotlib.pyplot as plt




def f(y,t=None):
    return -y*t


#===================
#euler method

def euler(y0,t0,h):
    ynext=y0+h*f(y0,t0)
    return ynext

y=np.zeros(11)
t=y.copy()

ycurr=10.0 #initial condition
tcurr=0.0 #integrate from zero
h=0.2     #step size

for i in range(11):
    y[i],t[i]=ycurr,tcurr
    ynew=euler(ycurr,tcurr,h)
    ycurr=ynew
    tcurr+=h

plt.plot(t,y,label='Euler')


#===================
#fourth order RK

y=np.zeros(11)
t=y.copy()

ycurr=10.0 #initial condition
tcurr=0.0 #integrate from zero
h=0.2     #step size


def RK4(y0,t0,h):
    k1=h*f(y0,t0)
    k2=h*f(y0+0.5*k1,t0+0.5*h)
    k3=h*f(y0+0.5*k2,t0+0.5*h)
    k4=h*f(y0+k3,t0+h)
    ynext=y0+(1./6.)*k1+(1./3.)*k2+(1./3.)*k3+(1./6.)*k4
    return ynext

for i in range(11):
    y[i],t[i]=ycurr,tcurr
    ynew=RK4(ycurr,tcurr,h)
    ycurr=ynew
    tcurr+=h

plt.plot(t,y,label='RK4')

plt.legend()
plt.show()
