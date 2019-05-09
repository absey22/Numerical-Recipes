import numpy as np
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt

from scipy.fftpack import fft2

wolfimage=plt.imread("wolf.jpg", format='jpg')

wolfimage=wolfimage.astype(float)

print(wolfimage.shape)

wolfft=fft2(wolfimage)
print(wolfimage.dtype)
plt.figure()
plt.subplot(1,2,1)
plt.imshow(wolfimage)
plt.subplot(1,2,2)
plt.imshow(wolfft)
plt.show()

wolftile=np.tile(wolfimage,(2,2))
wolftileft=fft2(wolftile)


plt.figure()
plt.subplot(1,2,1)
plt.imshow(wolftile)
plt.subplot(1,2,2)
plt.imshow(wolftileft)
plt.show()

