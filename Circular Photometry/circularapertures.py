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
t = input("reject ,accept, or frac:")

#Functions

#Circular Photometry

#Libraries
import numpy as np
from astropy.io import fits
from math import sin,cos,tan,atan2,sqrt,pi,radians as rad, degrees as deg, log10

#Inputs
im0 = "aptest.FIT" #input("Image file name:")
im0 = fits.getdata(im0)
x = 490 #int(input("x:"))
y = 293 #int(input("y:"))
im = im0[y-24:y+24,x-24:x+24]
rAperture = 5 #int(input("Aperture radius:"))
rAnnulusInner = 8 #int(input("Annulus radius inner:"))
rAnnulusOuter = 13 #int(input("Annulus radius outer:"))

#Functions
def apertureFracf(im,x,y,r):
    aperture = 0
    count = 0
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            for m in np.arange(i-.5,i+.5,.1):
                for n in np.arange(j-.5,j+.5,.1):
                    if (x - n)**2 + (y - m)**2 <= r**2:
                        aperture += (im[i,j]) * (1 / 10**2)
                        count += 1
    return aperture,count

aperture,count = aperturef(im,x,y,rAperture)
print(aperture,count)

def apertureRejectf(im,x,y,r,t):
    aperture = 0
    count = 0
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            for k in np.arange(-.5,.5,.1):
                m0 = ((j+.5)-x)**2 + ((i+k)-y)**2
                m1 = ((j-.5)-x)**2 + ((i+k)-y)**2
                m2 = ((j+k)-x)**2 + ((i-.5)-y)**2
                m3 = ((j+k)-x)**2 + ((i+.5)-y)**2
                if m0 <= r**2 and m1 <= r**2 and m2 <= r**2 and m3 <= r**2:
                    aperture += im[i,j]
                    count += 1
    return aperture,count


aperture,count = apertureFracf(im,x,y,rAperture,t)
print(aperture,count)
