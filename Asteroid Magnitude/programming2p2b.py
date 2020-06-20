def centroidAndUncertainty(image,xcord,ycord,length,width):

    import numpy as np
    from astropy.io import fits
    import matplotlib.pyplot as plt

    def x_centroid(matrix,y):
      num = 0
      for i in range(matrix.shape[1]):
        num = num + ((i+y-1) * matrix[:,i:i+1].sum())
      denom = matrix.sum()
      centroid = num / denom
      return centroid

    def y_centroid(matrix,x):
      num = 0
      for i in range(matrix.shape[0]):
        num = num + ((i+x-1) * matrix[i:i+1,:].sum())
      denom = matrix.sum()
      centroid = num / denom
      return centroid

    def x_sd(matrix,y):
      nums = []
      for i in range(matrix.shape[1]):
        nums.append((matrix[:,i:i+1].sum()) * ((i+y-1) - x_centroid(matrix,y))**2)
      num = sum(nums)
      denom = matrix.sum() - 1
      x_sd = (num / denom)**(.5)
      x_sd_mean = x_sd / (matrix.sum()**(.5))
      return x_sd_mean

    def y_sd(matrix,x):
      nums = []
      for i in range(matrix.shape[0]):
        nums.append((matrix[i:i+1,:].sum()) * ((i+x-1) - y_centroid(matrix,x))**2)
      num = sum(nums)
      denom = matrix.sum() - 1
      y_sd = (num / denom)**(.5)
      y_sd_mean = y_sd / (matrix.sum()**(.5))
      return y_sd_mean

    im = fits.getdata(image)
    plt.imshow(im)
    plt.show()
    x = int(xcord)
    y = int(ycord)
    l = int(length)
    w = int(width)
    #points = im[y-int((w-1)/2):y+int((w-1)/2)+1,x-int((l-1)/2):x+int((l-1)/2)+1]
    points = im[y - int(np.floor(w/2)):y + int(np.floor(w/2))+1, x - int(np.floor(l/2)):x + int(np.floor(l/2))+1]
    return ["x centroid: " + str(x_centroid(points,y)),"y_centroid: " + str(y_centroid(points,x)),"x_sd: " + str(x_sd(points,y)),"y_sd: " + str(y_sd(points,x))]

print(centroidAndUncertainty("2002Q15.00000001.Entered Coordinates.fit",371,741,8,8))
