import numpy as np
import matplotlib.pyplot as plt
from myfunctions import poisson,mlcg,xorshift



# ==========================  1(a)   ==========================
print(" ")
print(" ------- ")
print("Exercise 1(a):")
print(" ")
print("Poisson: P_lambda(k)")
print("P_1(0) =",poisson(1,0))
print("P_5(10) =",poisson(5,10))
print("P_3(20) =",poisson(3,20))
print("P_2.6(40) =",poisson(2.6,40))





# ==========================  1(b)   ==========================
#seed should be 1 < I_0 < m-1, prime, and an unsigned 64-bit
# FOUND IN "myfunctions.py"


#64-bit XOR-shift: the input must be 64-bit and not equal to 0
print(" ")
print(" ------- ")
print("Exercise 1(b):")
print(" ")


#take the LAST value from MLCG:
mlcg_out=mlcg(I_0)[1][-1]
#plug it into the XOR shift algorithm
rng_out=xorshift( mlcg_out )

print("The seed chosen:", I_0)
print("generated the random number:", rng_out)



#plot the resulting first 1000 floats generated:
x=xorshift(mlcg(I_0)[0])
y=xorshift(mlcg(I_0)[1])

plt.figure()
plt.subplot(1,2,1)
plt.plot(x[:1000],y[:1000],'ro')
plt.xlabel("$I_j$")
plt.ylabel("$I_{j+1}$")
plt.title("n = 1,000")

plt.subplot(1,2,2)
plt.hist(y,bins=20,edgecolor='black', linewidth=1.2)
plt.xlabel("binned value")
plt.ylabel("frequency")
plt.ylim(49000,50600)
#plt.yticks([49600,49700,49800,49900,50000,50100,50200,50300])

plt.title("n = 1,000,000")

plt.show()


