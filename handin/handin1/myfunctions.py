import numpy as np

#FUNCTIONS to be called by the other scripts of this hand in assignment



# ==========================  1(a)   ==========================
def poisson(lam,k):
    lam=np.float64(lam)
    k=np.int64(k)       #makes the fctn "continuous"
    if k<0:            #k cant be negative
        k=-1*k
    if lam<0:
        lam=-1.*lam
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


def rng(seed,niter=2000):
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

def rng_normalize(rnglist,low,high):
    rescaled=high+((rnglist-np.amax(rnglist))*(high-low)/np.ptp(rnglist))
    return rescaled #return list of rescaled random floats

#return the 3D integrand integrated via trapezoid rule, but
# including factor of 4*pi for integration over azimuth and elevation
# (the average number of satellites cancels)
def trapezoidrule_3Dnorm(func,panels,x1,x2,A,B,C):
    h=(x2-x1)/panels
    integral=0.5*(func(x1,A,B,C)+func(x2,A,B,C))
    for i in range(1,panels):
        integral+=func(x1+i*h,A,B,C)
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



# ==========================  2(d)   ==========================


def rejectionsample(rnglist,func,normconst,N_sat=100.):
    xsamples=[]
    x,y=rnglist #generate points 0,1
    x=rng_normalize(x,0.0,5.0) #normalize them to size of p(x)
    y=rng_normalize(y,0.0,np.nanmax(4*np.pi*normconst*func(x)/N_sat))
    evalpts=4*np.pi*normconst*func(x)/N_sat
    for i in range(len(x)):
        if y[i]<=evalpts[i]: #keep x if y falls below p(x)
            xsamples.append(x[i])
        if len(xsamples)==100: #only keep first 100 drawn
            break
    return xsamples



# ==========================  2(e)   ==========================


def createhaloes(func,normconst,N):
    haloes=np.empty((N,100,3))
    rnglist=rng(I_0,niter=2000)[1] #generate list of random numbers used for seeds
    cnt=0 #counter for supplying fresh seeds to each halo from rnglist
    for i in range(N):
        t=rng_normalize(rng(rnglist[cnt],niter=100)[1],0.0,np.pi) #elevation angles
        p=rng_normalize(rng(rnglist[cnt+1],niter=100)[1],0.0,2.*np.pi) #azimuthal angles
        
        r=rejectionsample(rng(rnglist[cnt+2],niter=1000),func,normconst)
        
        haloes[i]=np.stack((r,t,p),axis=1) #append generated halo
        cnt=i+3 #increment in the seed counter
    return  haloes



# ==========================  2(f)   ==========================


def bisection(fx,x1,x2,niter):
    xl=x1;xr=x2
    if fx(xl)*fx(xr)>0.:
        print("not bracketing a root")
        return np.nan
    else:
        for n in range(niter):
            xguess=(xl+xr)/2.
            cond=fx(xl)*fx(xguess) #check function values at these points
            if cond<0.:
                xr=xguess #cut off upper region
            if cond>0.:
                xl=xguess #cut off lower bracket region
            if abs(xr-xl)<=1e-12:  #accuracy to converge to
                print("Bisection converged to within 1e-12 after ",n,"iterations")
                break
            if n==niter:
                print("Bisection reached ",niter,"iterations without a root.")
        return xguess


    
# ==========================  2(g)   ==========================


def selectionsort(data): #perform Selection Sort on data
    sorteddata=data
    for i in range(0,len(data)-1):
        i_min=i
        for j in range(i+1,len(data)):
            if sorteddata[j]<sorteddata[i_min]:
                i_min=j
        if i_min!=i:
            sorteddata[i_min],sorteddata[i]=sorteddata[i],sorteddata[i_min]
    return sorteddata

def calcpercentile(data,perc): #calculate the request percentile perc in the data
    if perc>1:
        perc/=100.
    galcnt=len(data)
    return data[int(perc*galcnt)]



# ==========================  2(h)   ==========================


def lininterp2D_onedim(data,func,aa,bb,cc,ax,density):
    cubeinterp=[]
    for z in range(data.shape[0]):
        sheet=data[z,:,:] #define sheet in which to "2D interpolate"
        sheetinterp=[]
        for i in range(sheet.shape[0]): #iterate over one edge of the sheet
            if ax=="0": #"ax" keyword for which dimension is being calculated
                basedata=sheet[i,:]
                ydata=trapezoidrule_3Dnorm(func,panels=int(1e3),x1=0,x2=5,A=aa[z],B=bb[i],C=cc[:])
            elif ax=="1":
                basedata=sheet[int(i/density),:]
                ydata=trapezoidrule_3Dnorm(func,panels=int(1e3),x1=0,x2=5,A=aa[:],B=bb[i],C=cc[z])
            elif ax=="2":
                basedata=sheet[int(i/density),:]
                ydata=trapezoidrule_3Dnorm(func,panels=int(1e3),x1=0,x2=5,A=aa[z],B=bb[:],C=cc[i])
            sample_points=np.linspace(min(basedata),max(basedata),density*len(basedata))
            j_low=[] #determine matching indices for make interpolation points
            for bd in basedata:
                interp_point=np.argmin(abs(bd-sample_points))
                j_low.append(  interp_point  )
            j_low=np.asarray(j_low)
            xinterp=sample_points[ j_low[:].astype(int) ] #get those points in the sample_points
            sheetinterp.append(lininterp(xinterp,sample_points,j_low,ydata)) #sheet of 2D interpolation
        cubeinterp.append(sheetinterp) #add each sheet to the new interpolated axis of Acube
    return np.asarray(cubeinterp)



# ==========================  3(a)   ==========================

#loglikelihood to maximize
def loglikelihood(normconst,A,B,C,xpt):
    return np.log(normconst)+ (A-3.)*np.log(xpt/B)-(xpt/B)**C

def INTEGRANDdensityprofile(x,A,B,C):
    #x=np.float64(x)
    return (x/B)**(A-1.)*np.exp(-(x/B)**C)
