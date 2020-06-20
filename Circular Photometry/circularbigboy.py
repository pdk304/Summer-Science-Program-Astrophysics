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
D = 10
g = 2
readNoise = 15

#Functions
def apertureFracf(im,x,y,r):
    aperture = 0
    count = 0
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            if (j - x)**2 + (i - y)**2 <= (r+2)**2:
                for m in np.arange(i-.5,i+.5,.1):
                    for n in np.arange(j-.5,j+.5,.1):
                        if (x - n)**2 + (y - m)**2 <= r**2:
                            aperture += (im[i,j]) * (1 / 10**2)
                            count += 1/100
    return aperture,count

def aperturef(im,x,y,r):
    aperture = 0
    count = 0
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            if (j - x)**2 + (i - y)**2 <= (r+2)**2:
                m0 = ((j+.5)-x)**2 + ((i+.5)-y)**2
                m1 = ((j+.5)-x)**2 + ((i-.5)-y)**2
                m2 = ((j-.5)-x)**2 + ((i+.5)-y)**2
                m3 = ((j-.5)-x)**2 + ((i-.5)-y)**2
                if m0 <= r**2 and m1 <= r**2 and m2 <= r**2 and m3 <= r**2:
                    aperture += im[i,j]
                    count += 1
    return aperture,count

def apertureAcceptf(im,x,y,r):
    aperture = 0
    count = 0
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            if (j - x)**2 + (i - y)**2 <= (r+2)**2:
                m0 = ((j+.5)-x)**2 + ((i+.5)-y)**2
                m1 = ((j+.5)-x)**2 + ((i-.5)-y)**2
                m2 = ((j-.5)-x)**2 + ((i+.5)-y)**2
                m3 = ((j-.5)-x)**2 + ((i-.5)-y)**2
                if m0 <= r**2 or m1 <= r**2 or m2 <= r**2 or m3 <= r**2:
                    aperture += im[i,j]
                    count += 1
    return aperture,count

def avgAnnulusf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter):
    annulus = aperturef(im,x,y,rAnnulusOuter)[0] - aperturef(im,x,y,rAnnulusInner)[0]
    annulusPixels = aperturef(im,x,y,rAnnulusOuter)[1] - aperturef(im,x,y,rAnnulusInner)[1]
    avgAnnulus = annulus / annulusPixels
    return avgAnnulus

def signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter):
    objectPlusSky = aperturef(im,x,y,rAperture)[0]
    avgSky = avgAnnulusf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter)
    signal = objectPlusSky - avgSky * aperturef(im,x,y,rAperture)[1]
    return signal

def magf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter):
    signalObject = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter)
    mag = -2.5 * log10(signalObject)
    return mag

def SNRf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise):
    S = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter)
    n_ap = aperturef(im,x,y,rAperture)[1]
    n_an = aperturef(im,x,y,rAnnulusOuter)[1] - aperturef(im,x,y,rAnnulusInner)[1]
    rho = sqrt((readNoise)**2 + (g / sqrt(12))**2)
    SNR = sqrt((g * S) / sqrt(1 + n_ap * (1 + n_ap / n_an) * ((g * S + g * D + rho**2) / (g * S))))
    return SNR

def noisef(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise):
    S = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter)
    SNR = SNRf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise)
    noise = S / SNR
    return noise

def uncertaintyf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise):
    SNR = SNRf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise)
    uncertainty = 1.0857 / SNR
    return uncertainty

def rangef(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise):
    S = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter)
    noise = noisef(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise)
    inf = S - noise
    sup = S + noise
    return [inf,sup]

#Main
aperture,count = aperturef(im,x,y,rAperture)
print("aperture:",aperture,"aperture count:",count)

annulus = aperturef(im,x,y,rAnnulusOuter)[0] - aperturef(im,x,y,rAnnulusInner)[0]
print("annulus:",annulus)

signal = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter)
print("signal:",signal)

mag = magf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter)
print("mag:",mag)

SNR = SNRf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise)
print("SNR:",SNR)

noise = noisef(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise)
print("noise:",noise)

uncertainty = uncertaintyf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise)
print("uncertainty:",uncertainty)

range = rangef(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise)
print("range:",range)
