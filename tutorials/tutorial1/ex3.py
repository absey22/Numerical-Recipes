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

"""

#plt.plot(basedata,np.zeros(len(basedata)),'ro')
#plt.plot(xinterp,np.zeros(len(xinterp)),'bo')
#plt.plot(xinterp,ydata,'yo')


#LINEAR INTERPOLATION on these values:
lin_interpolation=[]
for i in range(len(xinterp)):
    print("section",i+1,"/",len(xinterp))
    if i<len(xinterp)-1:
        print("from:",min(sample_points[j_low[i]:j_low[i+1]]),"to:",max(sample_points[j_low[i]:j_low[i+1]]))
        y=ydata[i]+(sample_points[j_low[i]:j_low[i+1]]-xinterp[i])*( (ydata[i+1]-ydata[i]) / (xinterp[i+1]-xinterp[i]) )
        lin_interpolation.append(y)
    else:
        print("(must extrapolate)")
        
        #y=ydata[i-1]+(sample_points[j_low[i]:]-xinterp[i-1])*( (ydata[i]-ydata[i-1]) / (xinterp[i]-xinterp[i-1]) )
lin_interpolation=np.asarray(lin_interpolation)
lin_interpolation=lin_interpolation.ravel()



#plt.plot(sample_points[:len(lin_interpolation)],lin_interpolation,'g^')

"""


P_init=ydata
P=P_init



#NEVILLES ALGORITHM AT SINGLE POINTS

#def nevillealg(x0,x1,p1,p2):
#    P=( (x-x0)*p1 - (x-x1)*p2 )/(x1-x0)
#    return P
poly=[]
for s in sample_points:
    print "s:",s
    print "P:",P
    for order in range(len(xinterp)-1):
        print "order",order
        H=[]
        for i in range(len(P)-1):
            j=i+(order+1)
            print "(i,j):",(i,j)
            if order==0:
                result=P
            else:
                num1=(s - xinterp[j]) * P[i]
                num2=(s - xinterp[i]) * P[i+1]
                den=xinterp[i]-xinterp[j]
                result=( num1 - num2 )  / (den)
        H.append( result )
        print H
    P=H
poly.append(P)
"""
#NEVILLES ALGORITHM ON SET OF POINTS
for order in range(len(xinterp)-1):
    print
    print(":::FINDING ORDER",order+1," terms:::  (using order",order,"with",len(P),"entries)")
    for p in P:
        print(p)
    H=[]
    for i in range(len(P)-1):
        j=i+(order+1)
        print
        print("--==--the p-th item:","P_"+str((i,j)))
        if order==0:
            num1=(sample_points[j_low[i]:j_low[j]] - xinterp[j]) *P[i]
            num2=(sample_points[j_low[i]:j_low[j]] - xinterp[i]) *P[i+1]
            den=xinterp[j]-xinterp[i]
        else:
            print("1st term:","P_"+str((i,j-1)))
            print("---interpolate on:",min(sample_points[j_low[i]:j_low[j-1]]),"to",max(sample_points[j_low[i]:j_low[j-1]]),"btwn",xinterp[i],"and",xinterp[j-1],"(",len(sample_points[j_low[i]:j_low[j-1]]),"pts.)")
            print("---subtract",xinterp[j],"from it")
            print("---multiply that result by",P[i])
            
            print("2nd term:","P_"+str((i+1,j)))
            print("---interpolate on:",min(sample_points[j_low[i+1]:j_low[j]]),"to",max(sample_points[j_low[i+1]:j_low[j]]),"btwn",xinterp[i+1],"and",xinterp[j],"(",len(sample_points[j_low[i+1]:j_low[j]]),"pts.)")
            print("---subtract",xinterp[i],"from it")
            print("---multiply that result by",P[i+1])

            num1=(sample_points[j_low[i]:j_low[j-1]] - xinterp[j]) * P[i]
            num2=(sample_points[j_low[i+1]:j_low[j]] - xinterp[i]) * P[i+1]
            den=xinterp[i]-xinterp[j]
            
        result=( num1 - num2 )  / (den)
        print("RESULT =",result," (",len(result),"items)")
        H.append( result )
    #H=np.asarray(H)
    #P=H.ravel()
    P=H
print(len(P))

#plt.plot(sample_points[:len(P)],P,'p*')
#plt.show()
"""



"""

Neville=[]
init_P=ydata
P=init_P

for s in range(len(sample_points)):
    H=( (xinterp[i+1]-sample_points[])*P[i+1] + \
        (sample_points[]-xinterp[i]) * P[i]  )  \
       /  (xinterp[i]-xinterp[i+1])
    H.append()
    P=H
       
"""









#plt.show()

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



