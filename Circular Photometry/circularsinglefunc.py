#Circular Photometry

#Libraries
import numpy as np
from astropy.io import fits
from math import sin,cos,tan,atan2,sqrt,pi,radians as rad, degrees as deg, log10

#Inputs
im = "2002QF15.00000048.Entered Coordinates.fit" #input("Image file name:")
im = fits.getdata(im)
x = np.rint(621.9572212065814) #int(input("x:"))
y = np.rint(595.8402193784278) #int(input("y:"))
rAperture = np.rint(9 / 2) #int(input("Aperture radius:"))
rAnnulusInner = np.rint((8 / 5) * rAperture) #rAperture + 3 #int(input("Annulus radius inner:")) +3 rAperture
rAnnulusOuter = np.rint((13 / 8) * rAnnulusInner) #rAnnulusInner + 5 #int(input("Annulus radius outer:")) +5 rAnnulusInner
D = 10
g = 2
readNoise = 15
t = 0 #0 for frac, 1 for border reject, 2 for border accept
#Note: atom was crashing when using user inputs, so I hard-coded the values instead

#Functions
def aperturef(im,x,y,r,t):
    aperture = 0
    count = 0
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            if (j - x)**2 + (i - y)**2 <= (r+1)**2:
                m0 = ((j+.5)-x)**2 + ((i+.5)-y)**2
                m1 = ((j+.5)-x)**2 + ((i-.5)-y)**2
                m2 = ((j-.5)-x)**2 + ((i+.5)-y)**2
                m3 = ((j-.5)-x)**2 + ((i-.5)-y)**2
                if t == 0: #Fractional photometry; divides each pixel into 100 subpixels
                    for m in np.arange(i-.5,i+.5,.1):
                        for n in np.arange(j-.5,j+.5,.1):
                            if (x - n)**2 + (y - m)**2 <= r**2:
                                aperture += (im[i,j]) * (1 / 10**2)
                                count += 1/100
                if t == 1: #Border reject; checks each corner of each pixel; if all corners are inside the circle, accepts
                    if m0 <= r**2 and m1 <= r**2 and m2 <= r**2 and m3 <= r**2:
                        aperture += im[i,j]
                        count += 1
                if t == 2: #Border accept; checks each corner of each pixel; if one or more of them are inside the circle, accepts
                    if m0 <= r**2 or m1 <= r**2 or m2 <= r**2 or m3 <= r**2:
                        aperture += im[i,j]
                        count += 1
    return aperture,count

def avgAnnulusf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t):
    annulus = aperturef(im,x,y,rAnnulusOuter,t)[0] - aperturef(im,x,y,rAnnulusInner,t)[0]
    annulusPixels = aperturef(im,x,y,rAnnulusOuter,t)[1] - aperturef(im,x,y,rAnnulusInner,t)[1]
    avgAnnulus = annulus / annulusPixels
    return avgAnnulus

def signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t):
    objectPlusSky = aperturef(im,x,y,rAperture,t)[0]
    avgSky = avgAnnulusf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
    signal = objectPlusSky - avgSky * aperturef(im,x,y,rAperture,t)[1]
    return signal

def magf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t):
    signalObject = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
    mag = -2.5 * log10(signalObject)
    return mag

def SNRf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t):
    S = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
    Sky = avgAnnulusf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
    n_ap = aperturef(im,x,y,rAperture,t)[1]
    n_an = aperturef(im,x,y,rAnnulusOuter,t)[1] - aperturef(im,x,y,rAnnulusInner,t)[1]
    rho = sqrt((readNoise)**2 + (g / sqrt(12))**2)
    SNR = sqrt(g * S) / sqrt(1 + n_ap * (1 + n_ap / n_an) * ((g * Sky + g * D + rho**2) / (g * S)))
    return SNR

def noisef(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t):
    S = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
    SNR = SNRf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t)
    noise = S / SNR
    return noise

def uncertaintyf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t):
    SNR = SNRf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t)
    uncertainty = 1.0857 / SNR
    return uncertainty

def rangef(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t):
    S = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
    noise = noisef(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t)
    inf = S - noise
    sup = S + noise
    return [inf,sup]

#Main
signal = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
print("signal:",signal)

SNR = SNRf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t)
print("SNR:",SNR)

mag = magf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
print("mag:",mag)

uncertainty = uncertaintyf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t)
print("uncertainty:",uncertainty)
