import numpy as np


#FUNCTIONS to be called by the other scripts of this hand in assignment

# ==========================  1(a)   ==========================
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


# ==========================  1(b)   ==========================

I_0=np.int64(1014674096)

#M(CLG) with c=0: linear method for recursive number generation

#all the possible integers between 0 and m-1 will occur eventually.
#Starting from I_0, sequence takes off so that any successive I_j are random.
# modulus, and parameter values from Press et al. pp.349 Method E
#will have period modulus-1
def mlcg(seed=I_0,a=10014146,modulus=549755813881,niter=int(1e6)):
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




# ==========================  2(a)   ==========================

def rng_normalize(low,high):
    rnglist=xorshift(mlcg(I_0)[1])  #generate random floats 0 to 1.
    rescaled=high+((rnglist-np.amax(rnglist))*(high-low)/np.ptp(rnglist))
    return rescaled[-1] #return the last item in list of rescaled random floats
