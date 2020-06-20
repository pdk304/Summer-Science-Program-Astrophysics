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
def aperturef(im,x,y,r,t):
    aperture = 0
    count = 0
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            mat = np.full((10,10),im[i,j])
            if t == "frac":
                for m in range(mat.shape[0]):
                    for n in range(mat.shape[1]):
                        if (x - n)**2 + (y - m)**2 <= r**2:
                            aperture += (im[i,j]) * (1 / 10**2)
                            count += 1
            elif t == "reject":
                for k in np.arange(-.5,.5,.1):
                    m0 = ((j+.5)-x)**2 + ((i+k)-y)**2
                    m1 = ((j-.5)-x)**2 + ((i+k)-y)**2
                    m2 = ((j+k)-x)**2 + ((i-.5)-y)**2
                    m3 = ((j+k)-x)**2 + ((i+.5)-y)**2
                    if m0 <= r**2 and m1 <= r**2 and m2 <= r**2 and m3 <= r**2:
                        aperture += im[i,j]
                        count += 1
            else:
                for k in np.arange(-.5,.5,.1):
                    m0 = ((j+.5)-x)**2 + ((i+k)-y)**2
                    m1 = ((j-.5)-x)**2 + ((i+k)-y)**2
                    m2 = ((j+k)-x)**2 + ((i-.5)-y)**2
                    m3 = ((j+k)-x)**2 + ((i+.5)-y)**2
                    if m0 <= r**2 or m1 <= r**2 or m2 <= r**2 or m3 <= r**2:
                        aperture += im[i,j]
                        count += 1
    return aperture,count


aperture,count = aperturef(im,x,y,rAperture,t)
print(aperture,count)

annulus = aperturef(im,x,y,rAnnulusOuter,t)[0] - aperturef(im,x,y,rAnnulusInner,t)[0]
print(annulus)

def avgAnnulusf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t):
    annulus = aperturef(im,x,y,rAnnulusOuter,t)[0] - aperturef(im,x,y,rAnnulusInner,t)[0]
    annulusPixels = aperturef(im,x,y,rAnnulusOuter,t)[1] - aperturef(im,x,y,rAnnulusInner,t)[1]
    avgAnnulus = annulus / annulusPixels
    return avgAnnulus

def signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t):
    objectPlusSky = aperturef(im,x,y,rAperture,t)[0]
    avgSky = avgAnnulusf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
    print(avgSky)
    signal = objectPlusSky - avgSky * aperturef(im,x,y,rAperture)[1]
    return signal


signal = signalf(im,x,y,rAperture,annulus,rAnnulusInner,rAnnulusOuter,t)
print(signal)

def signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t):
    objectPlusSky = aperturef(im,x,y,rAperture,t)[0]
    avgSky = avgAnnulusf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
    print(avgSky)
    signal = objectPlusSky - avgSky * aperturef(im,x,y,rAperture,t)[1]
    return signal

signal = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
print(signal)

def magf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t):
    signalObject = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
    mag = -2.5 * log10(signalObject)
    return mag

mag = magf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
print(mag)
