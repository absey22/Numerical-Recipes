import numpy as np
from myfunctions import mlcg,xorshift,rng_normalize
#RNG seed found in "myfunctions.py"

#generate the free parameters controlling of the exp drop-off:
#   1.1 < a < 2.5
#   0.5 < b < 2
#   1.5 < c < 4

a=rng_normalize(1.1,2.5)
b=rng_normalize(0.5,2.0)
c=rng_normalize(1.5,4.0)

