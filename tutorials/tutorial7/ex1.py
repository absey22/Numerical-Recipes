import numpy as np
import matplotlib.pyplot as plt

xspace=np.linspace(0,4,20)

def f(x,a=2,b=4,c=2):
    return a*x+b*x**0.5+c


realization=np.tile(f(xspace)[:,np.newaxis],10000)

print(realization[:,:3])

noisyrealization=realization+np.random.normal(0,1,(len(xspace),10000))
print(noisyrealization[:,:3])


plt.plot(xspace,realization[:,0])

plt.plot(xspace,noisyrealization[:,0])

plt.show()
