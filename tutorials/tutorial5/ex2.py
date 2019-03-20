import numpy as np
import matplotlib.pyplot as plt


#functions to use-----------------------
def function1(x):
    f=x**3-6.*x**2+11.*x-6
    fprime=3.*x**2-12.*x+11
    return f,fprime

def function2(x):
    f=np.tan(np.pi*x)-6.
    fprime=np.pi*(np.cos(np.pi*x))**-2
    return f,fprime

def function3(x):
    f=x**3-2.*x+2
    fprime=3.*x**2-2
    return f,fprime

def function4(x):
    f=np.exp(10.*(x-1))-0.1
    fprime=10.*np.exp(10.*(x-1))
    return f,fprime
#---------------------------------------


#CONTROL PANEL: choose bracket and function
bounds=(2.5,4.0)
funx=function1



def secant(fx,x1,x2,niter):
    xprev=x1
    xcurr=x2
    for n in range(1,niter+1):
        #print(n,":",xprev,"-->",xcurr)
        try: #when at root, denom=0, so div by 0 means convergence at root
            xint=xprev-(xcurr-xprev)*fx(xprev)[0]/(fx(xcurr)[0]-fx(xprev)[0]) #[0] is f, [1] in f'
        except:
            print("Secant converged at i =",n)
            break
        xprev=xcurr #update xvalues for next iteration
        xcurr=xint
        if n==100:
            print("Secant reached ",niter,"iterations.")
    return xcurr

print(secant(funx,bounds[0],bounds[1],100))



def newtraph(fx,x1,x2,niter):
    #xlow=x1
    xprev=np.nan
    xcurr=x2
    for n in range(1,niter+1):
        #print(n,":",xprev,"-->",xcurr)
        xint=xcurr-(fx(xcurr)[0]/fx(xcurr)[1]) #[0] is f, [1] in f'
        xprev=xcurr
        xcurr=xint
        if abs(xcurr-xprev)<=1e-12:  #accuracy to converge to
            print("Newt-Raph converged at i =",n)
            break
        if n==100:
            print("Newt-Raph reached ",niter,"iterations.")
        #if xcurr<xlow:
        #    print("no root in interval")
        #    break
    return xcurr

print(newtraph(funx,bounds[0],bounds[1],100))


def falseposition(fx,x1,x2,niter):
    xcurr=x1
    for n in range(1,niter+1):
        #print(n,":",x2,"<--",xcurr)
        xint=x2-(fx(x2)[0]*(xcurr-x2)/(fx(xcurr)[0]-fx(x2)[0]))
        #print("num",fx(x1)[0]*(xcurr-x1))
        #print("denom1",fx(xcurr)[0])
        #print("denom2",-fx(x1)[0])
        if abs(xcurr-xint)<=1e-12: #accuracy to converge to
            print("False Pos. converged at i =",n)
            break
        xcurr=xint
        if n==100:
            print("False Position reached ",niter,"iterations.")
    return xcurr

print(falseposition(funx,bounds[0],bounds[1],10000))


def bisection(fx,x1,x2,niter):
    xl=x1;xr=x2
    if fx(xl)[0]*fx(xr)[0]>0.:
        print("not bracketing a root")        
    for n in range(1,niter+1):
        xguess=(xl+xr)/2.
        cond=fx(xl)[0]*fx(xguess)[0] #check function values at these points
        if cond<0.:
            xr=xguess #cut off upper region
        if cond>0.:
            xl=xguess #cut off lower bracket region
        if abs(xr-xl)<=1e-12:  #accuracy to converge to
            print("Bisection converged at i =",n)
            break
        if n==100:
            print("Bisection. reached ",niter,"iterations.")
    return xguess

print(bisection(funx,bounds[0],bounds[1],100))


def brents():
    return None

print(brents())
