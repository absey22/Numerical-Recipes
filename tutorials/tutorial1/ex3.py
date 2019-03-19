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






print("LINEAR INTERPOLATION")
lin_interpolation=[]
for i in range(len(xinterp)):
    #print("section",i+1,"/",len(xinterp))
    if i<len(xinterp)-1:
        #print("from:",min(sample_points[j_low[i]:j_low[i+1]]),"to:",max(sample_points[j_low[i]:j_low[i+1]]))
        y=ydata[i]+(sample_points[j_low[i]:j_low[i+1]]-xinterp[i])*( (ydata[i+1]-ydata[i]) / (xinterp[i+1]-xinterp[i]) )
        lin_interpolation.append(y)
    else:
        print("(must extrapolate)")
        
        #y=ydata[i-1]+(sample_points[j_low[i]:]-xinterp[i-1])*( (ydata[i]-ydata[i-1]) / (xinterp[i]-xinterp[i-1]) )
lin_interpolation=np.asarray(lin_interpolation)
lin_interpolation=lin_interpolation.ravel()

plt.figure()
plt.subplot(3,1,1)
plt.plot(basedata,np.zeros(len(basedata)),'ro',label="data abscissae")
plt.plot(xinterp,ydata,'yo',label="interp function values")
plt.plot(sample_points[:len(lin_interpolation)],lin_interpolation,'g^')
plt.legend()
plt.title("Linear Interpolation")






print("POLYNOMIAL INTERPOLATION (Neville's)")
P_init=ydata
P=P_init
def neville(interppts,sample_pt,P):
    for order in range(len(interppts)):
        H=[]
        #print "order",order
        if order==0:
            H=P
            #print H
            continue
        #print "BEFORE: P:",P, "len(P)",len(P)
        for i in range(len(P)-1):
            j=i+order
            #print("--==--the p-th item:","P_"+str((i,j)))
            #print P[i]
            num1=(sample_pt - interppts[j]) * P[i]
            num2=(sample_pt - interppts[i]) * P[i+1]
            den=interppts[i]-interppts[j]
            result=( num1 - num2 )  / (den)
            H.append( result )
            #print H
        Pfinal=H
        P=Pfinal
        #print "AFTER: P:",Pfinal, "len(P)",len(Pfinal)
    return Pfinal[0]

plt.subplot(3,1,2)
plt.plot(sample_points[:],neville(xinterp,sample_points[:],P),'bo')
#plt.plot(basedata,fx(basedata),'ro')
plt.ylim(-30,30)
plt.plot(basedata,np.zeros(len(basedata)),'ro',label="data abscissae")
plt.plot(xinterp,ydata,'yo',label="interp function values")
plt.legend()
plt.title("Neville's Algorithm")


"""
#NEVILLES ALGORITHM ON SIMULTANEOUS SET OF POINTS (doesnt work)
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


print("POLYNOMIAL INTERPOLATION (Natural Cubic Splines)")

print("# of abssicae: "+str(len(xinterp)))

a=ydata
print a, len(a)

h=[]
for i in range(len(xinterp)-1):
    h.append(xinterp[i+1]-xinterp[i])

print "h length", len(h)
beta=[]
for i in range(1,len(xinterp)-1):
    
    term1=(3./h[i])*(a[i+1]-a[i])
    term2=(3./h[i-1])*(a[i]-a[i-1])
    beta.append(term1-term2)
print "beta length", len(beta)

l=[];mu=[];z=[]
l.append(1.);mu.append(0.);z.append(0.)
for i in range(1,len(xinterp)-1):
    l.append(2*(xinterp[i+1]-xinterp[i-1])-(h[i-1]*mu[i-1]))
    mu.append(h[i]/l[i])
    z.append((beta[i-1]-h[i-1]*z[i-1])/l[i])
    

#boundary condition
l.append(1.);z.append(0.)

b=np.zeros(len(xinterp))
c=np.zeros(len(xinterp))
d=np.zeros(len(xinterp))
#c[-1]=0.

print "len(z)",len(z)
print "len(mu)",len(mu)
print "len(c)",len(c)
for i in range(len(xinterp)-2,0,-1):
    print i
    c[i]=z[i]-mu[i]*c[i-1]
    b[i]=((a[i+1]-a[i])/h[i])-(h[i]*(c[i+1]+2.*c[i])/3.)
    d[i]=(c[i+1]-c[i])/(3.*h[i])

print a
print b
print c
print d



plt.subplot(3,1,3)
plt.plot(basedata,np.zeros(len(basedata)),'ro',label="data abscissae")
plt.plot(xinterp,ydata,'yo',label="interp function values")
plt.title("Natural Cubic Splines")
plt.legend()


#plt.show()
