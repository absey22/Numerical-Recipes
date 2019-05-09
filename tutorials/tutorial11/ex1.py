import numpy as np

Mstar=1*Msolar
Ma=1*Mj;Pa=12yr
Mb=0.011*Mj;Pb=1yr

def euler_r(r0,t0,h):
    rnext=r0+h*euler_v(r0,t0)
    return rnext

def euler_v(r0,t0,h):
    def grav(M,m):
        return G*(M+m)/r**3.*euler_r
    vnext=v0+h*grav(v0,t0)
    return vnext



r1,r2=np.zeros((100,2))
t=np.zeros(100)

r1curr=(G*(mj+solar)/(4.*np.pi**2.))**1./3.*(12yrs); r2curr=r1curr #initial condition
tcurr=0.0 #integrate from zero
h=0.2     #step size

for i in range(100):
    r1[i],r2[i]=r1curr,r2curr
    ynew=euler(ycurr,tcurr,h)
    ycurr=ynew
    tcurr+=h
