import numpy as np

#----1(a)----

#create list
numbers=[]
i=0
while i<100:
   i+=1
   numbers.append(i)

#calculate
avg=sum(numbers)/len(numbers)
stddev=np.std(numbers)

print("avg",avg)
print("stddev",stddev)


#----1(b)----


#for n in numbers:
#    if n%2==0:
        
even=[n%2==0 for n in numbers]
odd=np.logical_not(even)

if even:
    
if odd:
        

