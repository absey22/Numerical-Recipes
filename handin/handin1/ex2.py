import numpy as np
import matplotlib.pyplot as plt


#RNG seed found in "myfunctions.py"




# ==========================  2(a)   ==========================
from myfunctions import I_0,rng,rng_normalize,trapezoidrule_3Dnorm,neville,lininterp

#generate the free parameters controlling of the exp drop-off:
#   1.1 < a < 2.5
#   0.5 < b < 2.0
#   1.5 < c < 4.0

a=rng(I_0,100)[1] #seed a list of random numbers
a=rng_normalize(a,1.1,2.5)[53] #take an arbitrary item from this list after renormalizing
b=rng(I_0,100)[1]
b=rng_normalize(b,0.5,2.0)[11]
c=rng(I_0,100)[1]
c=rng_normalize(c,1.5,4.0)[78]


#define the profile function
def densityprofile(x,normalization,N_sat=100.):
    return normalization*N_sat*(x/b)**(a-3.)*np.exp(-(x/b)**c)


#Assuming densityprofile is 3D symmetric about origin then
# volume integral reduces to 4pi*INT(densityprofile*x^2*dx)
# (change the power to which the x prefactor is raised to in
#  densityprofile() in order to form the integrand)
def INTEGRANDdensityprofile(x,A=a,B=b,C=c):
    #x=np.float64(x)
    return (x/B)**(A-1.)*np.exp(-(x/B)**C)


# Via implementing the trapezoid rule, calculate the normalization constant
# (the average number of satellites cancels)
A = trapezoidrule_3Dnorm(INTEGRANDdensityprofile,panels=int(1e5),x1=0,x2=5,A=a,B=b,C=c)

#display result:
print(" ")
print(" ------- ")
print("Exercise 2(a):")
print(" ")
print("Using generated  a = %f" % a,", b = %f" % b, ", and c = %f" % c)
print("  The profile normalization constant is A = %f" % A)



# ==========================  2(b)   ==========================
from myfunctions import neville,lininterp
plt.clf()
plt.tight_layout()

#make the initial data points
xpts=np.asarray([1e-4,1e-2,1e-1,1.,5.])
xdata=np.log(xpts)
ydata=np.log(densityprofile(xpts,A))

#plot base data and profile behavior
plt.figure()
plt.subplot(2,1,1)
plt.title("Log Profile Interpolations")
plt.plot(xdata,ydata,'r^',markersize=10,label='"Data" values')
tempspace=np.linspace(xdata[0],xdata[-1],25)
plt.plot(tempspace,np.log(densityprofile(np.exp(tempspace),A)),'b+',label='Profile behavior')



#interpolating:

#Ringing of polynomial interpolation at edges due to high order polynomial
# (large number of points).

sample_points=np.linspace(xdata[0],xdata[-1],40)

#find the interpolation values in sample_points
j_low=[]
for x in xdata:
    interp_point=np.argmin(abs(x-sample_points))
    j_low.append(  interp_point  )
j_low=np.asarray(j_low)

#define the matching points in sample_points to interpolate at
xinterp=sample_points[ j_low[:].astype(int) ]

#initilize the polynomial points for Nevilles (explicitly)
P_init=ydata
P=P_init

#interpolate and plot
poly_interp=neville(xinterp,sample_points,P)
plt.plot(sample_points,poly_interp,'bo',label="Neville's (4th degree polynomial)")

plt.legend()


#Spline interpolation would overcome polynomial oscillations,
# but linear interpolation captures similar behavior with the
# least amount of required computation.

#interpolate and plot
plt.subplot(2,1,2)
plt.plot(xdata,ydata,'r^',markersize=10,label='"Data" values')
tempspace=np.linspace(xdata[0],xdata[-1],25)
plt.plot(tempspace,np.log(densityprofile(np.exp(tempspace),A)),'b+',label='Profile behavior')
lin_interp=lininterp(xinterp,sample_points,j_low,ydata)
plt.plot(sample_points[:len(lin_interp)],lin_interp,'go',label="Linear Interpolation")

#display the results:
print(" ")
print(" ------- ")
print("Exercise 2(b): (plotted)")
print(" ")

plt.xlabel("log(x)")
plt.ylabel("log( n(x) )")
plt.xlim(xdata[0]-0.1,xdata[-1]+0.1)


plt.legend()
plt.savefig("./plots/interpolationcomparison.png")




# ==========================  2(c)   ==========================
from myfunctions import ridders
plt.clf()

print(" ")
print(" ------- ")
print("Exercise 2(c):")
print(" ")


print("Via Ridder's method in 20 approximations:")
print("  dn(x)/dx (@ x=b) =",ridders(densityprofile,normconst=A,evalpt=b))

from sympy import Symbol,exp,diff
x = Symbol('x')
f = 100.*A*(x/b)**(a-3.)*exp(-(x/b)**c)
fprime = f.diff(x)

print("Via Sympy:")
print("  dn(x)/dx (@ x=b) =", fprime.subs(x,b))



# ==========================  2(d)   ==========================
from myfunctions import rejectionsample
plt.clf()
plt.tight_layout()

print(" ")
print(" ------- ")
print("Exercise 2(d):")
print(" ")


theta,phi=rng(I_0,100)
theta=rng_normalize(theta,0.0,np.pi) #elevation angles
phi=rng_normalize(phi,0.0,2.*np.pi) #azimuthal angles


r=rejectionsample(rng(I_0,niter=1000),INTEGRANDdensityprofile,normconst=A)

#display the result:
print("Via rejection sampling using uniform generators to create radii given the known profile:")
print(" (   r(x/xvir)   ,  theta(rad)   ,  phi(rad)    )")
coords=np.stack((r,theta,phi),axis=1)
print(coords)



# ==========================  2(e)   ==========================
from myfunctions import createhaloes
plt.clf()
plt.tight_layout()

print(" ")
print(" ------- ")
print("Exercise 2(e): (plotted)")
print(" ")


haloescube=createhaloes(INTEGRANDdensityprofile,normconst=A,N=1000)

#get the radial component of all N haloes and histogram plot them
radialcomponent=haloescube[:,:,0]

#20 evenly log-spaced bins AND "density=True" ensures the area (or integral) under the histogram will sum to 1 by dividing the count by the number of observations times the bin width and not dividing by the total number of observations.

evenlogbins=10.**np.linspace(-4,np.log(5.)/np.log(10.),20)
#binwidths=logbins[1:]-logbins[:-1] #calculate width of each bin in log space
#make histogram AND evenly space it by normalizing by the binwidths
#radialhist,binedges=np.histogram(radialcomponent.ravel(),logbins)
#radialhistnorm=radialhist/binwidths
plt.hist(radialcomponent.ravel(),bins=evenlogbins,edgecolor='black',linewidth=1.2,density=True,label= '(prob. density normalized counts)')
#plt.bar(logbins[:-1],radialhistnorm,binwidths,label='Avg. satellites in 1000 haloes (normalized)')

#plot the profile itself N(x)=n(x)4pi*x^2
xspace=np.linspace(1e-4,5,10000)
plt.plot(xspace,4*np.pi*A*INTEGRANDdensityprofile(xspace),":",lw=3,label='Profile of Dist.')

plt.title('Satellites in 1000 haloes')
plt.ylabel("log( N(x) )")
plt.xlabel("log(x)")
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.9e-4,5.1)
plt.legend(loc=3)

plt.savefig("./plots/satellitesof1000haloes.png")



# ==========================  2(f)   ==========================
from myfunctions import bisection
plt.clf()
plt.tight_layout()

print(" ")
print(" ------- ")
print("Exercise 2(f): ")
print(" ")


#find roots as intersection(s) of y/2 and N(x)

#First get the maximum and its location
y=np.amax(4*np.pi*A*INTEGRANDdensityprofile(xspace))
locmax=np.argmax(4*np.pi*A*INTEGRANDdensityprofile(xspace))

#By inspection, the profile of the distribution has one roots on either
# side of the maximum at y. There, call bisection once for each of them
# by specifying the bounds of the root finding on either side of y
region1=(1e-4,xspace[locmax])
region2=(xspace[locmax],5.)

#by shifting the profile down by y/2, we can call bisection on it.
def shiftedprofile(xpt):
    return 4*np.pi*A*INTEGRANDdensityprofile(xpt)-(y/2.)

#call bisection method on each region
root1=bisection(shiftedprofile,region1[0],region1[1],1000)
root2=bisection(shiftedprofile,region2[0],region2[1],1000)


#display the result:
print("Two roots were found:")
print("1: x,y =",root1,",",4*np.pi*A*INTEGRANDdensityprofile(root1))
print("2: x,y =",root2,",",4*np.pi*A*INTEGRANDdensityprofile(root2))

plt.plot(root1,INTEGRANDdensityprofile(root1),'y*',label='roots')
plt.plot(root2,INTEGRANDdensityprofile(root2),'y*')
plt.plot(region1[0],0,"r<",label='region 1')
plt.plot(region1[1]-0.1,0,"r>")
plt.plot(region2[0]+0.1,0,"b<",label='region 2')
plt.plot(region2[1],0,"b>")
plt.plot(xspace,INTEGRANDdensityprofile(xspace))
plt.hlines(y/2,xspace[0],xspace[-1],linestyle=":")
plt.legend()


plt.savefig("./plots/bi-intersection.png")




# ==========================  2(g)   ==========================
from myfunctions import selectionsort,calcpercentile,poisson
plt.clf()
plt.tight_layout()

print(" ")
print(" ------- ")
print("Exercise 2(g): ")
print(" ")


#---Sorting galaxies:

#using radial data from (e) 
radialhist,binedges=np.histogram(radialcomponent.ravel(),20)

#the largest bin occurs at
ind=np.argmax(radialhist)
#print(binedges[ind],binedges[ind+1])
#print(radialhist[ind])

#crop the haloes array to this bin
haloes=[]
for i in range(radialcomponent.shape[0]):
    halo=radialcomponent[i,:] #clip each halo to that largest bin
    clippedhalo=halo[(halo>binedges[ind])&(halo<binedges[ind+1])]
    haloes.append(clippedhalo)


galaxies=[gal for gal in haloes for gal in gal] #sort them
sortedgalaxies=selectionsort(galaxies.copy())

#16th, 84th percentile (-1, +1 sigma around mean) and median of sortedgalaxies:
print("16th percentile:",calcpercentile(sortedgalaxies,16))
print("50th percentile (median):",calcpercentile(sortedgalaxies,50))
print("84th percentile:",calcpercentile(sortedgalaxies,84))



#---Poisson comparison:

heights=[len(galaxies) for galaxies in haloes]

plt.hist(heights,bins=24,edgecolor='black',linewidth=1.2,density=True,label='(prob. density normalized, max = '+str(np.amax(heights))+')')
k=np.arange(np.amin(heights),np.amax(heights))     #measurement k's for plotting
lamb=np.mean(heights) #mean of poisson distribution to plot
poissondist=[poisson(lamb,k) for k in k]

plt.title("No. of Galaxies in Largest Radial Bin")
plt.ylabel("Normalized Counts")
plt.xlabel("")
plt.plot(k,poissondist,":",lw=3,label='Poisson dist., mean = '+str(lamb))

plt.legend()
plt.savefig("./plots/poissoncomparison.png")



# ==========================  2(h)   ==========================
from myfunctions import lininterp2D_onedim
plt.clf()
plt.tight_layout()

print(" ")
print(" ------- ")
print("Exercise 2(h): ")
print(" ")

aa=np.arange(1.1,2.6,0.1) #15
bb=np.arange(0.5,2.1,0.1) #16
cc=np.arange(1.5,4.1,0.1) #26


#Evaluate the normalization constant over a grid of a, b, and c parameters to create a cub
Acube=trapezoidrule_3Dnorm(INTEGRANDdensityprofile,panels=int(1e5),x1=0,x2=5,A=aa[:,None,None],B=bb[None,:,None],C=cc[None,None,:])


print("3D Interpolation between",Acube.shape[0]*Acube.shape[1]*Acube.shape[2],"original values.")

print("Starting from",Acube.shape[0],"sheets of dimension",Acube.shape[1],"by",Acube.shape[2])




density=2 #change the density parameter to determine the stretching factor of all dimensions

adense=np.arange(1.1,2.6,0.1/density)
bdense=np.arange(0.5,2.1,0.1/density)
cdense=np.arange(1.5,4.1,0.1/density)
#Do Linear interpolation over ONE DIMENSION of Acube. "density parameter determines
# how many multiples of data points to interpolate into Acube with.
# (e.g. density=2 will double the size of the datacube)
# "ax" determines which axis you are stretching

#first dimension is over "c" parameter in Acube
Acubeinterp1=lininterp2D_onedim(Acube,INTEGRANDdensityprofile,adense,bdense,cdense,ax="0",density=density)


print("Interpolating in c leaves",Acubeinterp1.shape[0],"sheets of dimension",Acubeinterp1.shape[1],"by",Acubeinterp1.shape[2])


plt.figure()
plt.subplot(3,2,1)
plt.title("Original Cube")
plt.imshow(Acube[7,:,:],origin='lower')
plt.subplot(3,2,2)
plt.title("Interpolated Along 1st Dimension")
plt.xlabel("(here)")
plt.imshow(Acubeinterp1[7,:,:],origin='lower')

plt.subplot(3,2,3)

plt.imshow(Acube[:,8,:],origin='lower')
plt.subplot(3,2,4)
plt.xlabel("(here)")
plt.imshow(Acubeinterp1[:,8,:],origin='lower')

plt.subplot(3,2,5)

plt.imshow(Acube[:,:,15],origin='lower')
plt.subplot(3,2,6)

plt.imshow(Acubeinterp1[:,:,25],origin='lower')

plt.savefig("./plots/interpolationfirstdim.png")


#Must transpose the result in order to work over next axis of the array Acube(interp)
Acubeinterp1_flipped=np.transpose(Acubeinterp1,(2,1,0))

#interpolate over "a" parameter
Acubeinterp2=lininterp2D_onedim(Acubeinterp1_flipped,INTEGRANDdensityprofile,adense,bdense,cdense,ax="1",density=density)

print("Interpolating in a leaves",Acubeinterp2.shape[2],"sheets of dimension",Acubeinterp2.shape[1],"by",Acubeinterp2.shape[0])



plt.figure()
plt.subplot(3,2,1)
plt.title("Original Cube")
plt.imshow(Acube[7,:,:],origin='lower')
plt.subplot(3,2,2)
plt.title("Interpolated Along 1st AND 2nd Dimension")
plt.xlabel("(here)")
plt.imshow(np.transpose(Acubeinterp2,(2,1,0))[7,:,:],origin='lower')

plt.subplot(3,2,3)

plt.imshow(Acube[:,8,:],origin='lower')
plt.subplot(3,2,4)
plt.xlabel("(here)")
plt.imshow(np.transpose(Acubeinterp2,(2,1,0))[:,8,:],origin='lower')

plt.subplot(3,2,5)
plt.imshow(Acube[:,:,15],origin='lower')
plt.subplot(3,2,6)
plt.ylabel("(here)")
plt.imshow(np.transpose(Acubeinterp2,(2,1,0))[:,:,25],origin='lower')


plt.savefig("./plots/interpolationseconddim.png")


#transpose to a new face for the interpolation
Acubeinterp2_flipped=np.transpose(Acubeinterp2,(2,0,1))

#interpolate over "b" parameter
Acubeinterp3=lininterp2D_onedim(Acubeinterp2_flipped,INTEGRANDdensityprofile,adense,bdense,cdense,ax="2",density=density)

print("Interpolating in b leaves",Acubeinterp2.shape[1],"sheets of dimension",Acubeinterp2.shape[0],"which doesn't make sense. It should still be 3D.")

"""
plt.figure()
plt.subplot(3,2,1)
plt.title("Original Cube")
plt.imshow(Acube[7,:,:])
plt.subplot(3,2,2)
plt.title("Interpolated Along 1st, 2nd, AND 3rd Dimensions")
plt.imshow(np.transpose(Acubeinterp3,(2,0,1))[7,:,:])

plt.subplot(3,2,3)

plt.imshow(Acube[:,8,:])
plt.subplot(3,2,4)

plt.imshow(np.transpose(Acubeinterp3,(2,0,1))[:,8,:])

plt.subplot(3,2,5)

plt.imshow(Acube[:,:,15])
plt.subplot(3,2,6)

plt.imshow(np.transpose(Acubeinterp3,(2,0,1))[:,:,25])

plt.show()
"""

#**************************************************************************************
#WEIRD BUG 3D Interpolation: I couldnt figure out why my interpolation fails for the
# final dimension. The line fitting returns empty arrays from lininterpolation within
# lininterp2D_onedim (therefore final dimension is not displayed).
#**************************************************************************************
