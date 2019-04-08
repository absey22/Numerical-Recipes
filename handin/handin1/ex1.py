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
print("P_3(21) =",poisson(3,21))
print("P_2.6(40) =",poisson(2.6,40))
#print("P_101(200) =",poisson(101,200))



# ==========================  1(b)   ==========================

# SEED FOUND IN "myfunctions.py"



print(" ")
print(" ------- ")
print("Exercise 1(b):")
print(" ")


#generate 1e6 numbers from mLCG/XORshift iterations for binning
x,y=rng(I_0,niter=int(1e6))

#take the LAST value of one of the sequential list:
rng_out=y[-1]

print("The seed", I_0,"generates a random number", rng_out)



#plot the resulting 1000 floats generated:

plt.figure()
plt.subplot(1,2,1)
plt.plot(x[10:1010],y[10:1010],'ro')
plt.xlabel("$I_j$")
plt.ylabel("$I_{j+1}$")
plt.title("Sequential Random Numbers (n = 1,000)")

#plot 1,000,000 random numbers in a histogram

plt.subplot(1,2,2)
freq,bins,patches=plt.hist(y,bins=20,edgecolor='black', linewidth=1.2)
binmidpts=0.5*(bins[1:] + bins[:-1])                  #find mid point of each bin
plt.errorbar(binmidpts, freq, yerr=(freq)**0.5, fmt='none',label='Poissonian Std. dev.')
plt.hlines(50000,0,1,linestyles=':',label='Ideal Uniform Dist.')
plt.title("Binning random generations (n = 1,000,000)")
plt.xlabel("pseudorandom value (0.05 wide bins)")
plt.ylabel("frequency")
plt.ylim(45000,52000)
#plt.yticks([49600,49700,49800,49900,50000,50100,50200,50300])

plt.legend()
plt.savefig("./plots/rngquality.png")


