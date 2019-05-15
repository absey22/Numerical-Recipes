import numpy as np
import matplotlib.pyplot as plt

"""

# ==========================  1(a)   ==========================
from myfunctions2 import Urng

#make a random realization of 1e6 uniformly drawn scalars
uniformdraw=Urng(size=int(1e6))

#(FORMATTING PARTIALLY BORROWED FROM MY EXERCISE 1)
#plot the resulting 1000 floats generated:
plt.figure()
plt.subplot(2,2,1)
plt.plot(uniformdraw[0:1000:2],uniformdraw[1:1001:2],'ro')
plt.xlabel("$x_i$")
plt.ylabel("$x_{i+1}$")
plt.title("$n_{draws}$ = 1,000")

#plot index vs. magnitude
plt.subplot(2,2,2)
plt.xlabel("index i")
plt.ylabel("$x_i$")
plt.title("$n_{draws}$ = 1,000")
plt.bar(np.arange(len(uniformdraw[0:1000])),uniformdraw[0:1000])

#histogram the 1,000,000 random numbers
plt.subplot(2,1,2)
freq,bins,_=plt.hist(uniformdraw,bins=20,edgecolor='black', linewidth=1.2)
binmidpts=0.5*(bins[1:] + bins[:-1])                  #find mid point of each bin
plt.errorbar(binmidpts, freq, yerr=(freq)**0.5, fmt='none',label='Poissonian Std. dev.')
plt.hlines(50000,0,1,linestyles=':',label='Ideal Uniform')
plt.title("MWC & Xor Shift Uniform Variate ($n_{draws}$ = 1,000,000)")
plt.xlabel("True values (0.05 wide bins)")
plt.ylabel("Frequency")
plt.ylim(46000,52000)
plt.legend(loc=8)

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig("./plots/uniformityanalysis.png")
plt.clf()



# ==========================  1(b)   ==========================
from myfunctions2 import BoxMuller,GaussianPDF

#perform Box Muller tranform from these uniform variates to a standard normal variate
# This function draws two lists of uniform variates from within
# (transformed again to given mean=3.0 and stdev=2.4)
mean=3.0
stddev=2.4
normdraw=BoxMuller(N=1000,mean=mean,stdv=stddev)

#plot the result
freq,bins,_=plt.hist(normdraw,bins=int(9*stddev),edgecolor='black', linewidth=1.5,density=True)
#list of stddev locations on x-axis
sigmabars=np.concatenate((mean-np.linspace(5*stddev,stddev,5),mean+np.linspace(stddev,5*stddev,5)))
#plot stddev locations and mean
plt.vlines(mean,0,GaussianPDF(mean,mean,stddev),linewidth=1.5,linestyles='--')
plt.vlines(sigmabars,0,GaussianPDF(sigmabars,mean,stddev),linewidth=1.8,linestyles=':',color='r')
plt.title("Box-Muller Normal Variate ($n_{draws}$ = 1,000,000)")
plt.xlabel("True values")
plt.ylabel("Pseudorandom (pdf-normalized) Frequency")
#plot expected shape of ideal normal
plt.plot(bins, GaussianPDF(bins,mean,stddev),linewidth=2.2,linestyle="--",color='k',label='Ideal Normal ($\mu$='+str(mean)+', $\sigma$='+str(stddev)+')')
plt.legend()

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig("./plots/gaussiancomparison.png")
plt.clf()



# ==========================  1(c)   ==========================
from myfunctions2 import mykstest,BoxMuller
from scipy.stats import kstest

#Prove: null hypothesis is that actual data follows a normal distribution
# accomplish by finding max distance between Gaussian CDF and data CDF

#The KS statistic D can be calculated as shown in Section 6.14 of Press et al.
# a p-value of the signficance of the D of this observation can therefore be found
# from that distribution.

N=int(1e5) # number of points to reach up to in testing
dex=0.1 # to increment by ("n_points=10^dex")
start,stop=1,5 # to start at and stop at (stop=log_10(N))
cnt=np.ones(int((stop-start)/dex)+1) # +1 to reach up to N in slicing randomdata
pval=np.zeros((len(cnt),2)) # 2 for accepting a D AND a p-value
Dval=pval.copy() # the same data shape will be required to store the D statistic
data=BoxMuller(N,mean=0,stdv=1.0)
for i in range(len(cnt)):
    cnt[i]+=round(i*dex,2)   # increment the cnt by dex each loop for slicing & plotting
    #normaldraw=BoxMuller(10**cnt[i],mean=0,stdv=1.0) #create a realization of data drawn from my Box-Muller algorithm (normal variates)
    dataslice=data[:int(10**cnt[i])]  # sucessively larger slices into randomdata to analyze
    Dval[i,0],pval[i,0]=mykstest(dataslice) # take D,pval from mykstest() (default: normal cdf)
    Dval[i,1],pval[i,1]=kstest(dataslice,'norm') # take D,pval from scipy.stats.kstest


plt.subplot(2,1,1)
plt.title("K-S Test")
plt.plot(cnt,pval[:,0],label="My p-value")
plt.plot(cnt,pval[:,1],":",label="SciPy p-value")
plt.ylabel("p-value")
plt.legend()

plt.subplot(2,1,2)
plt.plot(cnt,Dval[:,0],label="My D statistic")
plt.plot(cnt,Dval[:,1],":",label="SciPy D statistic")
plt.xlabel("dex ($N_{points} = 10^{dex}$)")
plt.ylabel("D statistic")
plt.legend()

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig("./plots/kstest.png")
plt.clf()

"""

# ==========================  1(d)   ==========================
from myfunctions2 import kuiperstest,BoxMuller
from scipy.stats import kstest

#See 1(c) for comments explaining these steps here:
N=int(1e5) 
dex=0.1 
start,stop=1,5 
cnt=np.ones(int((stop-start)/dex)+1)
mykuiper=np.zeros((len(cnt),2))
scipykuiper=mykuiper.copy()
randomdata=BoxMuller(N,mean=0,stdv=1.0)
for i in range(len(cnt)):
    cnt[i]+=round(i*dex,2) 
    #normaldraw=BoxMuller(10**cnt[i],mean=0,stdv=1.0) #create a realization of data drawn from my Box-Muller algorithm (normal variates)
    dataslice=randomdata[:int(10**cnt[i])]
    mykuiper[i]=kuiperstest(dataslice) 
    scipykuiper[i,0]=kstest(dataslice,'norm',alternative='greater')[0] 
    scipykuiper[i,1]=kstest(dataslice,'norm',alternative='less')[0]

    
#make the same plot as 1(c):
plt.figure()
plt.subplot(2,2,1)
plt.plot(cnt,mykuiper[:,0],label="My Kuiper test, D+")
plt.plot(cnt,scipykuiper[:,0],label="Scipy Kuiper test, D+")
plt.xlabel("dex ($N_{points} = 10^{dex}$)")
plt.ylabel("D+ statistic")
plt.legend()

plt.subplot(2,2,2)
plt.plot(cnt,mykuiper[:,1],label="My Kuiper D- statistic")
plt.plot(cnt,scipykuiper[:,1],label="Scipy D- statistic")
plt.xlabel("dex ($N_{points} = 10^{dex}$)")
plt.ylabel("D- statistic")
plt.legend()

plt.subplot(2,1,2)
plt.title("Kuiper Test   (V $\equiv$ D- + D+)")
plt.plot(cnt,np.add(mykuiper[:,0],mykuiper[:,1]),label="My Kuiper V statistic")
plt.plot(cnt,np.add(scipykuiper[:,0],scipykuiper[:,1]),label="Scipy V statistic")
plt.xlabel("dex ($N_{points} = 10^{dex}$)")
plt.ylabel("V statistic")
plt.legend()

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig("./plots/kuipertest.png")
plt.clf()


"""
# ==========================  1(e)   ==========================
from myfunctions2 import mykstest

data=np.genfromtxt("randomnumbers.txt").T #transpose it for convenience

dex=0.1 
start,stop=1,5 
cnt=np.ones(int((stop-start)/dex)+1)
mykuiper=np.zeros((len(cnt),2))
scipykuiper=mykuiper.copy()
randomdata=BoxMuller(N,mean=0,stdv=1.0)
for i in range(len(cnt)):
    cnt[i]+=round(i*dex,2) 
    #normaldraw=BoxMuller(10**cnt[i],mean=0,stdv=1.0) #create a realization of data drawn from my Box-Muller algorithm (normal variates)
    dataslice=randomdata[:int(10**cnt[i])]
    mykuiper[i]=kuiperstest(dataslice) 
    scipykuiper[i,0]=kstest(dataslice,'norm',alternative='greater')[0] 
    scipykuiper[i,1]=kstest(dataslice,'norm',alternative='less')[0]

    
#make the same plot as 1(c):
plt.subplot(2,1,1)
plt.title("Kuiper Test")
plt.plot(cnt,mykuiper[:,0],label="My Kuiper test, D+")
plt.plot(cnt,scipykuiper[:,0],label="Scipy Kuiper test, D+")
plt.ylabel("statistic")
plt.legend()

plt.subplot(2,1,2)
plt.plot(cnt,mykuiper[:,1],label="My Kuiper test, D-")
plt.plot(cnt,scipykuiper[:,1],label="Scipy Kuiper test, D-")
plt.xlabel("dex ($N_{points} = 10^{dex}$)")
plt.ylabel("statistic")
plt.legend()

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig("./plots/kuipertest.png")
plt.clf()

"""
