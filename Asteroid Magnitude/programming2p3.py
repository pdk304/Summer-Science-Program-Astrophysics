import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import math

def magAsteroid(image,xStar,yStar,wStar,magStar,xSky,ySky,xAsteroid,yAsteroid,wAsteroid):
    
    im = fits.getdata(image)
    
    starPlusSky = im[yStar-int(np.floor(wStar/2))-2:yStar+int(np.floor(wStar/2))+3,xStar-int(np.floor(wStar/2))-2:xStar+int(np.floor(wStar/2))+3]

    theSky = im[int(ySky)-1:int(ySky)+2,int(xSky)-1:int(xSky)+2]
    print("theSky",theSky)

    avgSky = theSky.sum() / theSky.size
    print("avgSky",avgSky)

    signalStar = starPlusSky.sum() - avgSky * starPlusSky.size
    print("signalStar",signalStar)

    const = magStar + 2.5 * math.log10(signalStar)
    print("const",const)

    asteroidPlusSky = im[yAsteroid-int(np.floor(wAsteroid/2))-2:yAsteroid+int(np.floor(wAsteroid/2))+3,xAsteroid-int(np.floor(wAsteroid/2))-2:xAsteroid+int(np.floor(wAsteroid/2))+3]

    signalAsteroid = asteroidPlusSky.sum() - avgSky * asteroidPlusSky.size
    print("signalAsteroid",signalAsteroid)

    magAsteroid = -2.5 * math.log10(signalAsteroid) + const

    return magAsteroid

print(magAsteroid("sampleimage.fits",173,342,5,15.26,200,200,351,154,3))

print(magAsteroid("sampleimage.fits",355,285,5,16.11,200,200,351,154,3))
