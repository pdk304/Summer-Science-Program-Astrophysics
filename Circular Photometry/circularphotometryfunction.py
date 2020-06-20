#Circular Photometry

#Libraries
import numpy as np
from astropy.io import fits
from math import sin,cos,tan,atan2,sqrt,pi,radians as rad, degrees as deg, log10

#Inputs
im = "2002Q15.00000030.ENTERED COORDINATES.REDUCED.FIT" #input("Image file name:")
im = fits.getdata(im)
x = np.rint(499.0895265423242) #int(input("x:"))
y = np.rint(687.0763271162124) #int(input("y:"))
rAperture = 6 #np.rint(9 / 2) #int(input("Aperture radius:"))
rAnnulusInner = 9 #np.rint((8 / 5) * rAperture) #rAperture + 3 #int(input("Annulus radius inner:")) +3 rAperture
rAnnulusOuter = 10 #np.rint((13 / 8) * rAnnulusInner) #rAnnulusInner + 5 #int(input("Annulus radius outer:")) +5 rAnnulusInner
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
'''
signal = signalf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
print("signal:",signal)

SNR = SNRf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t)
print("SNR:",SNR)

mag = magf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,t)
print("mag:",mag)

uncertainty = uncertaintyf(im,x,y,rAperture,rAnnulusInner,rAnnulusOuter,D,g,readNoise,t)
print("uncertainty:",uncertainty)
'''
def signalsf(im,x_list,y_list,rAperture,rAnnulusInner,rAnnulusOuter,t):
    signals = []
    for i in range(len(x_list)):
        signals.append(signalf(im,x_list[i],y_list[i],rAperture,rAnnulusInner,rAnnulusOuter,t))
    return signals

x_list4 = 166.1018483,347.8688243,930.1476355,519.2293351,866.2144675
y_list4 = 627.668265,881.7914428,754.1228374,491.0830784,1015.259147

x_list3 = 521.3751743,596.3180284,732.861394,743.7307128,734.4137686,748.2197092
y_list3 = 695.1496314,712.3581195,760.8003862,718.5251572,521.3150287,505.1943712

x_list2 = 398.5837535,302.1589271,501.40925,780.3292657,918.4169246,743.0968666
y_list2 = 955.5780907,634.8886794,474.4724623,362.3354559,732.9526193,836.8442667

x_list = 684.2154073,424.8698794,660.1534239,349.7381557,296.8370833
y_list = 437.1670219,849.1703208,677.0409993,557.9992269,456.2862673

signals = signalsf(im,x_list,y_list,rAperture,rAnnulusInner,rAnnulusOuter,t)
print(signals)
