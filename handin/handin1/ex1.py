import numpy as np
import matplotlib.pyplot as plt


def poisson(lam,k):
    k=int(k)           #makes the fctn "continuous"
    if k<0:            #k cant be negative
        k=-1*k
    lam=int(lam)
    if lam<0:
        lam=-1*lam
    if k==0:           
        kfact=1.0
    else:
        kfact=k        #calculate the factorial of k
        for i in range(k-1,0,-1):
            kfact*=1.*i
    P=(lam**(1.*k)*np.exp(-lam))/kfact
    return P


#M(CLG) with c=0: linear method for recursive number generation

#all the possible integers between 0 and m-1 will occur eventually.
#Starting from I_0, sequence takes off so that any successive I_j are random.
# modulus, and parameter values from Press et al. pp.349 Method E
#will have period modulus-1
def mlcg(seed,a=10014146,modulus=549755813881,niter=int(1e6)):
    a=np.int64(a)                     #convert to 64-bit
    modulus=np.int64(modulus)
    I_j=seed                            #supply the seed
    pI_j=[]
    pI_j1=[]
    for i in range(niter):
        I_j1=(a*I_j) % modulus          #update via the m(LCG) method
        pI_j.append(I_j)          #save each iteration for plotting
        pI_j1.append(I_j1)
        I_j=I_j1                        #reassign the value to change it for the next iteration
    return np.asarray(pI_j,dtype=np.uint64), np.asarray(pI_j1,dtype=np.uint64) #return the two lists for plotting


#this has a period of up to 2^64 - 1

#using parameter a1,a2,a3 from the lecture slides
def xorshift(val,a=(21,35,4)):
    a=np.uint64(a)                   #recast the shifting parameters to avoid
    val=val^(val>>a[0])              #  stuffing val into a smaller data type
    val=val^(val<<a[1])
    val=val^(val>>a[2])
    val=val/((2**64)-1)              #normalize this to a float btwn 0.0 and 1.0
    return val

#===============================================================================




#display results:
print(" ")
print(" ------- ")
print("Exercise 1(a):")
print(" ")
print("Poisson: P_lambda(k)")
print("P_1(0) =",poisson(1,0))
print("P_5(10) =",poisson(5,10))
print("P_3(20) =",poisson(3,20))
print("P_2.6(40) =",poisson(2.6,40))






#seed should be 1 < I_0 < m-1, prime, and an unsigned 64-bit
I_0=np.int64(1014674096)


#64-bit XOR-shift: the input must be 64-bit and not equal to 0
print(" ")
print(" ------- ")
print("Exercise 1(b):")
print(" ")

#display results:
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


