import numpy as np


def function1(x):
    return x**3.-6.*x**2.+11.*x-6.

def function2(x):
    return np.tan(np.pi*x)-6.




def secant(fx,x1,x2,niter):
    xprev=x1
    xcurr=x2
    for n in range(1,niter+1):
        print(n,":",xprev,"-->",xcurr)
        try:
            xint=xprev-(xcurr-xprev)*fx(xprev)/(fx(xcurr)-fx(xprev))
        except:
            break
        xprev=xcurr
        xcurr=xint
    return xcurr


print(secant(function1,2.5,4.0,100))



    
    
