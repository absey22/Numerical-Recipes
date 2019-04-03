import numpy as np

#FUNCTIONS to be called by the other scripts of this hand in assignment



# ==========================  1(a)   ==========================
def poisson(lam,k):
    lam=np.float64(lam)
    k=np.int64(k)           #makes the fctn "continuous"
    if k<0:            #k cant be negative
        k=-1*k
    if lam<0:
        lam=-1*lam
    if k==0:           
        kfact=1.0
    else:
        kfact=k       #calculate the factorial of k
        for i in range(k-1,0,-1):
            kfact*=1.*i
    P=(lam**(1.*k)*np.exp(-lam))/kfact
    return P


# ==========================  1(b)   ==========================
#seed should be 1 < I_0 < m-1, prime, and an unsigned 64-bit
I_0=np.uint64(7919)

#m(CLG) with c=0: linear method for recursive number generation

#all the possible integers between 0 and m-1 will occur eventually.
#Starting from I_0, sequence takes off so that any successive I_j are random.
# modulus, and parameter values from Press et al. pp.349 Method E
#will have period modulus-1
def mlcg(val,a=np.int64(10014146),modulus=np.int64(549755813881)):
    return (a*val) % modulus        


#64-bit XOR-shift: the input must be 64-bit and not equal to 0

#this has a period of up to 2^64 - 1
#using parameter a1,a2,a3 from the lecture slides
def xorshift(val,a=(21,35,4)):
    a=np.uint64(a)
    val=np.uint64(val)                   #recast the shifting parameters to avoid
    val=val^(val>>a[0])                  #  stuffing val into a smaller data type with XOR shift
    val=val^(val<<a[1])
    val=val^(val>>a[2])             
    return val                           #normalize this to a float btwn 0.0 and 1.0


def rng(seed,niter=int(1e6)):
    I_j=seed #supply the seed
    pI_j=[]
    pI_j1=[]
    for i in range(niter): 
        I_j1=xorshift(mlcg(I_j))         #apply m(LCG) method and XORshift that result
        pI_j.append(I_j)
        pI_j1.append(I_j1)
        I_j=I_j1                         #store this iteration
    return np.asarray(pI_j,dtype=np.uint64)/((2**64)-1) , np.asarray(pI_j1,dtype=np.uint64)/((2**64)-1)  #return the two lists for plotting



# ==========================  2(a)   ==========================

def rng_normalize(low,high):
    rnglist=rng(I_0,1000)[1]  #generate random floats 0 to 1.
    rescaled=high+((rnglist-np.amax(rnglist))*(high-low)/np.ptp(rnglist))
    return rescaled #return the last item in list of rescaled random floats

#return the 3D integrand integrated via trapezoid rule, but
# including factor of 4*pi for integration over azimuth and elevation
# (the average number of satellites cancels)
def trapezoidrule_3Dnorm(func,panels,x1,x2):
    h=(x2-x1)/panels
    integral=0.5*(func(x1)+func(x2))
    for i in range(1,panels):
        integral+=func(x1+i*h)
    return 1./(4*np.pi*(h*integral))


# ==========================  2(b)   ==========================

def neville(interppts,sample_pt,P):
    for order in range(len(interppts)):
        H=[]
        if order==0:
            H=P
            continue
        for i in range(len(P)-1):
            j=i+order
            num1=(sample_pt - interppts[j]) * P[i]
            num2=(sample_pt - interppts[i]) * P[i+1]
            den=interppts[i]-interppts[j]
            result=( num1 - num2 )  / (den)
            H.append( result )
        Pfinal=H
        P=Pfinal
    return Pfinal[0]

def lininterp(interppts,sample_pt,matchpts,ydata):
    interp=[]
    for i in range(len(interppts)-1):
        #print(" ")
        #print("datapt",xinterp[i],ydata[i])
        #print("samppts",sample_points[j_low[i]:j_low[i+1]])
        y=ydata[i] + (sample_pt[matchpts[i]:matchpts[i+1]]-interppts[i]) *    \
                     ( (ydata[i+1]-ydata[i]) / (interppts[i+1]-interppts[i]) )
        interp.append(y)
    return np.hstack(np.asarray(interp))



# ==========================  2(c)   ==========================


def centdiffapprox(func,evalpt,h,normconst):
    return (func(evalpt+h,normconst)-func(evalpt-h,normconst))/(2*h)


def ridders(func,normconst,evalpt,h_init=0.1,m=20,d=2.):
    approx=np.empty(m)
    h=h_init
    for i in range(m):
        approx[i]=centdiffapprox(func,evalpt,h,normconst) #calc the approx
        h/=d #divide h by factor d
    
    while len(approx)>1:  #combine approximations
        newapprox=np.empty(len(approx)-1)
        for i in range(len(approx)-1):
            D= (d**(2.*(i+1))*approx[i+1]-approx[i]) / (d**(2.*(i+1))-1)
            newapprox[i]=D
        approx=newapprox
    return approx[0]
