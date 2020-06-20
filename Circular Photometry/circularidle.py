#Circular Photometry

#Libraries
import numpy as np
from astropy.io import fits
from math import sin,cos,tan,atan2,sqrt,pi,radians as rad, degrees as deg, log10

#Inputs
im = "aptest.FIT" #input("Image file name:")
im = fits.getdata(im)
x = 490 #int(input("x:"))
y = 293 #int(input("y:"))
rAperture = 5 #int(input("Aperture radius:"))
rAnnulusInner = 8 #int(input("Annulus radius inner:"))
rAnnulusOuter = 13 #int(input("Annulus radius outer:"))

#Functions
def aperturef(im,x,y,r):
    aperture = 0
    count = 0
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            mat = np.full((10,10),im[i,j])
            for m in range(mat.shape[0]):
                for n in range(mat.shape[1]):
                    if (x - n)**2 + (y - m)**2 <= r**2:
                        aperture += (im[i,j]) * (1 / 10**2)
                        count += 1
    return aperture,count

aperture,count = aperturef(im,x,y,rAperture)
print(aperture,count)
