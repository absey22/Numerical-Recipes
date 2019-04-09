import numpy as np


# ==========================  3(a)   ==========================
from myfunctions import trapezoidrule_3Dnorm,loglikelihood,INTEGRANDdensityprofile

print(" ")
print(" ------- ")
print("Exercise 3(a): ")
print(" ")

fnames=["./satgals_m11.txt","./satgals_m11.txt","./satgals_m11.txt","./satgals_m11.txt","./satgals_m11.txt","./satgals_m11.txt"]

haloes=[]
for i in range(len(fnames)):
    haloes.append(np.genfromtxt(fnames[i],dtype=float,skip_header=5))




#Attempt to calculate the maximum would involve evaluating the likelihood under the parameters A,b,c and in each galaxy with xi position. This would give a 5D array which you could locate the maximum of, and find out what the values of A,b,c are there.

#find value of likelihood there: (use numpy generation for sake of time)
#this has to be done in each mass bin:

larray=[]
for i in range(len(haloes)):
    halo=haloes[i]
    x=halo[:,0]
    #generate random a,b,c at this point (*use numpy package)
    a=np.random.uniform(1.1,2.5)
    b=np.random.uniform(0.5,2.0)
    c=np.random.uniform(1.5,4.0)
    #calculating the normalization constant there
    A = trapezoidrule_3Dnorm(INTEGRANDdensityprofile,panels=int(1e5),x1=0,x2=5,A=a,B=b,C=c)
    calculation=[loglikelihood(A,a,b,c,x),[a,b,c]]
    larray.append(calculation) #make list of likelihood and a,b,c there

for l in larray: #loop through each mass bin
    massbin=np.asarray(l)
    litem=np.asarray(l[0])
    print("likelihood max:",massbin[np.argmax(litem[0])]) #take a,b,c corresponding to max likelihood

    #this should just be one value and output the corresponding a,b,ccd plot
    



