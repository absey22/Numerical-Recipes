import numpy as np
import matplotlib.pyplot as plt

basedata=np.linspace(1,20,7)


def fx(data):
    return 2*data*np.sin(0.8*data)



sample_points=np.linspace(1,25,100)

#find the points in the sample points mesh that lie adjacent to each of the data
# points. I think in this case that will be M=3 points per data point in the
# sample point space.
# One to the left of the data point, the data point itself, and the one
# right of it. "j_lo" is that the first of those: the left one.


j_low=[]
for x in basedata:
    interp_point=np.argmin(abs(x-sample_points))
    j_low.append(  interp_point  )
j_low=np.asarray(j_low)


xinterp=sample_points[ j_low[:].astype(int) ]
ydata=fx(xinterp)

plt.plot(basedata,np.zeros(len(basedata)),'ro')
plt.plot(xinterp,np.zeros(len(xinterp)),'bo')
plt.plot(xinterp,ydata,'yo')


#linear interpolate on these values:
lin_interpolation=[]
for i in range(len(xinterp)):
    print "section",i+1,"/",len(xinterp)
    if i<len(xinterp)-1:
        print "from:",min(sample_points[j_low[i]:j_low[i+1]]),"to:",max(sample_points[j_low[i]:j_low[i+1]])
        y=ydata[i]+(sample_points[j_low[i]:j_low[i+1]]-xinterp[i])*( (ydata[i+1]-ydata[i]) / (xinterp[i+1]-xinterp[i]) )
        lin_interpolation.append(y)
    else:
        print "(must extrapolate)"
        print i
        #y=ydata[i-1]+(sample_points[j_low[i]:]-xinterp[i-1])*( (ydata[i]-ydata[i-1]) / (xinterp[i]-xinterp[i-1]) )
lin_interpolation=np.asarray(lin_interpolation)
lin_interpolation=lin_interpolation.ravel()

print sample_points.shape,len(lin_interpolation)

plt.plot(sample_points[:len(lin_interpolation)],lin_interpolation,'g^')



Neville=[]
init_P=ydata
P=init_P
Htemp=[]
for i in range(len(xinterp)):
    H=( (xinterp[i+1]-sample_points)*P[i+1] + (sample_points-xinterp[i])*P[i] ) / (xinterp[i]-xinterp[i+1])
    Htemp.append(H)
    P=Htemp
       










plt.show()

"""  
# do this via bisection algorithm:

    #loop over each data point  
    #set edge points in the sample_points initially to the min and max
    edgelow=min(sample_points)
    edgehigh=max(sample_points)
    #find the mid point of this by len/2
    midpoint=(edgelow+edgehigh)/2
    #check if the data point is < or > that mid point
    if d<midpoint:
        edgehigh=midpoint
    #reset edge point(s) accordingly (move min to half, or max to half)
    #set new mid point, and check where data point is again.
    #repeat until data point either is a sample point OR is in btwn two sample points.
"""



