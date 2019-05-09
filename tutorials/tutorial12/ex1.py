import numpy as np
import matplotlib.pyplot as plt

# has m=442 sample stars with X and Y coordinates given n=2 labels, or features
dataset=np.genfromtxt("dataset_LATTE.txt")


features=dataset[:,:2]
m=features.shape[0]
n=features[0].shape[0]



#feature scaling
def rescale(x):
        return (x-np.mean(x))/np.std(x)
rescaledfeatures=np.zeros(features.shape)
for i in range(len(features[0,:])):
    rescaledfeatures[:,i]=rescale(features[:,i])

#plt.figure()
#plt.subplot(2,1,1)
#plt.plot(dataset[:,0],dataset[:,1],'ro')
#plt.subplot(2,1,2)
#plt.plot(rescaledfeatures[:,0],rescaledfeatures[:,1],'bo')
#plt.show()



#overwrite the features
features=rescaledfeatures


#gradient descent


#logistic regression


def hypothesis(data,parameters):
    p=parameters
    bestestimate=p[0]+1. + p[1]*data[:,0] + p[2]*data[:,1]
    return 1./( 1.+np.exp(-1.*bestestimate) )
    #return bestestimate

def lossfunction(bestestimate,label):#label=dataset[:,2]
    #bestestimate is the hypothesis
    return -1.* ( label * np.log(bestestimate) + (1.-label)*np.log(1.-bestestimate) )

def costfunction(lossfunction):
    return (1./m)*np.sum(lossfunction)

def gradientdescent(parameter,hypothesis,label,sample,learningrate=0.1):
    temp=np.empty((sample.shape[0],3))
    temp[:,0]=np.ones(sample.shape[0])
    temp[:,1:]=sample
    sample=temp
    return parameter - (learningrate/m)*np.sum((hypothesis-label)[:,None]*sample,axis=0)
    #return parameter - (learningrate/m)*np.sum(np.matmul((hypothesis-label),sample),axis=0)

#initial parameters
parameters=np.asarray([0.,0.,0.])

 #get labels of binary classification
labels=dataset[:,2]


#create initial best estimate yhat:
h_i=hypothesis(features,parameters)
#compute the loss function:
loss_i=lossfunction(h_i,labels)
#calculate the resulting cost:
cost_i=costfunction(loss_i)

err_th=1e-6
err=cost_i
errlist=[]
#update parameters via GD
while err>err_th:
    #compute new parameters:
    parameters=gradientdescent(parameters,h_i,labels,features)  
    #create new estimate yhat:
    h_new=hypothesis(features,parameters)
    #compute the loss function:
    loss_new=lossfunction(h_new,labels)
    #calculate the resulting cost:
    cost_new=costfunction(loss_new)
    
    err=abs(cost_i-cost_new)
    cost_i=cost_new
    errlist.append(err)

print(parameters)
plt.plot(np.arange(len(errlist)),errlist)
plt.show()
