import numpy as np
import matplotlib.pyplot as plt
from myfunctions import poisson,I_0,rng



# ==========================  1(a)   ==========================
print(" ")
print(" ------- ")
print("Exercise 1(a):")
print(" ")
print("Poisson: P_lambda(k)")
print("P_1(0) =",poisson(1,0))
print("P_5(10) =",poisson(5,10))
print("P_3(20) =",poisson(3,21))
print("P_2.6(40) =",poisson(2.6,40))





# ==========================  1(b)   ==========================

# SEED FOUND IN "myfunctions.py"



print(" ")
print(" ------- ")
print("Exercise 1(b):")
print(" ")


#generate 1e6 numbers from mLCG/XORshift iterations for binning
x,y=rng(I_0)

#take the LAST value of one of the sequential list:
rng_out=y[-1]

print("The seed", I_0,"generates a random number", rng_out)



#plot the resulting first 1000 floats generated:

plt.figure()
plt.subplot(1,2,1)
plt.plot(x[10:1010],y[10:1010],'ro')
plt.xlabel("$I_j$")
plt.ylabel("$I_{j+1}$")
plt.title("n = 1,000")

plt.subplot(1,2,2)
plt.hist(y,bins=20,edgecolor='black', linewidth=1.2)
plt.xlabel("binned pseudorandom value")
plt.ylabel("frequency")
plt.ylim(48500,51500)
#plt.yticks([49600,49700,49800,49900,50000,50100,50200,50300])

plt.title("n = 1,000,000")

plt.show()


