import numpy as np

A=1.*np.asarray([[3,8,1,-12,-4],[1,0,0,-1,-0],[4,4,3,-40,-3],[0,2,1,-3,-2],[0,1,0,-12,-0]])

b=np.asarray([[2,0,1,0,0]])


#at each column, find the pivot (argmax) (if none, exit as singular matrix)
Atemp=np.zeros(A.shape)
acnt=0

for i in range(A.shape[1]):
    print(A[:,i])
    #pivloc=np.where(np.any(A[:,i])==1.,np.argmin(A[:,i]-1.0),np.argmax(A[i:,i])+acnt)
    pivloc=int(A[:,i][A[:,i]==1.0][0])
    print("    i =",i,"  pivot loc:",pivloc)
    print(A[i:,i])
    rowholder=A[i,:].copy()
    pivotrow=A[pivloc,:].copy()
    pivotrow/=pivotrow[i]
    print("swap the pivot row:")
    print(pivotrow)
    print("and the current row:")
    print(rowholder)
    print(" ")
    A[i,:]=pivotrow
    A[pivloc,:]=rowholder
    Atemp[i,:]=pivotrow
   
    #print(A)
    acnt+=1

A=Atemp
print(A)


#scale the row so the pivot =1.0 (scale b too)

#swap that row to the right place in A (and b)

