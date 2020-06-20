import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

mat = np.array([[0,33,21,33,8],[0,56,51,53,26],[23,120,149,73,18],[55,101,116,50,16],[11,78,26,2,10]])

print(mat)

def x_centroid(matrix):
  num = 0
  for i in range(matrix.shape[1]):
    num = num + (i * matrix[:,i:i+1].sum())
  denom = matrix.sum()
  centroid = num / denom
  return centroid

def y_centroid(matrix):
  num = 0
  for i in range(matrix.shape[0]):
    num = num + (i * matrix[i:i+1,:].sum())
  denom = matrix.sum()
  centroid = num / denom
  return centroid

def centroid(matrix):
  centroid = x_centroid(matrix),y_centroid(matrix)
  return centroid

print(centroid(mat))
print(np.std(mat,axis = 0))

def x_sd(matrix):
  nums = []
  for i in range(matrix.shape[1]):
    nums.append((matrix[:,i:i+1].sum()) * (i - x_centroid(matrix))**2)
  num = sum(nums)
  denom = matrix.sum() - 1
  x_sd = (num / denom)**(.5)
  x_sd_mean = x_sd / (matrix.sum()**(.5))
  return x_sd_mean

def y_sd(matrix):
  nums = []
  for i in range(matrix.shape[0]):
    nums.append((matrix[i:i+1,:].sum()) * (i - y_centroid(matrix))**2)
  num = sum(nums)
  denom = matrix.sum() - 1
  y_sd = (num / denom)**(.5)
  y_sd_mean = y_sd / (matrix.sum()**(.5))
  return y_sd_mean

def uncertainty(matrix):
  uncertainty = x_sd(matrix),y_sd(matrix)
  return uncertainty

print(uncertainty(mat))

im = fits.getdata(input("FITS file name:"))
plt.imshow(im)
plt.show()
print(im)
x = int(input("Specify x coordinate:"))
y = int(input("Specify y coordinate:"))
points = im[x-1:x+2,y-1:y+2]
print(points)
print(centroid(points))
print(uncertainty(points))
