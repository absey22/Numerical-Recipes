import numpy as np
import matplotlib.pyplot as plt


#RNG seed found in "myfunctions.py"




# ==========================  2(a)   ==========================
from myfunctions import I_0,rng,rng_normalize,trapezoidrule_3Dnorm,neville,lininterp

#generate the free parameters controlling of the exp drop-off:
#   1.1 < a < 2.5
#   0.5 < b < 2.0
#   1.5 < c < 4.0

a=rng(I_0,100)[1]
a=rng_normalize(a,1.1,2.5)[52]
b=rng(I_0,100)[1]
b=rng_normalize(b,0.5,2.0)[10]
c=rng(I_0,100)[1]
c=rng_normalize(c,1.5,4.0)[78]


#define the profile function
def densityprofile(x,normalization,N_sat=100.):
    return normalization*N_sat*(x/b)**(a-3.)*np.exp(-(x/b)**c)


#Assuming densityprofile is 3D symmetric about origin then
# volume integral reduces to 4pi*INT(densityprofile*x^2*dx)
# (change the power to which the x prefactor is raised to in
#  densityprofile() in order to form the integrand)
def INTEGRANDdensityprofile(x):
    #x=np.float64(x)
    return (x/b)**(a-1.)*np.exp(-(x/b)**c)


# Via implementing the trapezoid rule, calculate the normalization constant
# (the average number of satellites cancels)
A = trapezoidrule_3Dnorm(INTEGRANDdensityprofile,panels=int(1e5),x1=0,x2=5)

#display result:
print(" ")
print(" ------- ")
print("Exercise 2(a):")
print(" ")
print("Using generated  a = %f" % a,", b = %f" % b, ", and c = %f" % c)
print("  The profile normalization constant is A = %f" % A)



# ==========================  2(b)   ==========================
from myfunctions import neville,lininterp

#make the initial data points
xpts=np.asarray([1e-4,1e-2,1e-1,1.,5.])
xdata=np.log(xpts)
ydata=np.log(densityprofile(xpts,A))

#plot base data and profile behavior
plt.figure()
plt.subplot(2,1,1)
plt.title("Log Profile Interpolations")
plt.plot(xdata,ydata,'r^',markersize=10,label='Interpolation values')
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
plt.plot(xdata,ydata,'r^',markersize=10,label='Interpolation values')
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

plt.show()




# ==========================  2(c)   ==========================
from myfunctions import ridders

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

print(" ")
print(" ------- ")
print("Exercise 2(d):")
print(" ")


theta,phi=rng(I_0,100)
theta=rng_normalize(theta,0.0,180.0) #elevation angles
phi=rng_normalize(phi,0.0,360.0) #azimuthal angles


r=rejectionsample(rng(I_0,niter=2000),INTEGRANDdensityprofile,normconst=A)

#display the result:
print(" (      r       ,     theta    ,      phi     )")
coords=np.stack((r,theta,phi),axis=1)
print(coords)



# ==========================  2(e)   ==========================
from myfunctions import createhaloes

print(" ")
print(" ------- ")
print("Exercise 2(e): (plotted)")
print(" ")


haloescube=createhaloes(INTEGRANDdensityprofile,normconst=A,N=1000)



#plt.plot(haloescube[100,:,0],4*np.pi*A*INTEGRANDdensityprofile(haloescube[100,:,0]),'ro')
#plt.plot(haloescube[210,:,0],4*np.pi*A*INTEGRANDdensityprofile(haloescube[210,:,0]),'go')
#plt.plot(haloescube[920,:,0],4*np.pi*A*INTEGRANDdensityprofile(haloescube[920,:,0]),'bo')

plt.hist(haloescube[:,:,0].ravel(),bins=10.**np.linspace(np.log10(1e-4),np.log10(5.),20))

xspace=np.linspace(0,6,1000)
plt.plot(xspace,90000*np.pi*A*INTEGRANDdensityprofile(xspace),":")

plt.ylabel("log( N(x) )")
plt.xlabel("log(x)")
plt.xscale("log")
plt.yscale("log")
#plt.xlim(-4-0.1,np.log(6.)+0.1)
plt.show()
