#LSPR code

import math
from math import sin, cos, atan2, radians as rad, degrees as deg, sqrt, atan
import numpy as np

x_centroid,y_centroid = 617.0468466545726,598.0791804503255
L = (3911*10**-3)/(24*10**-6)

#opening file
file = np.loadtxt("Obs2.txt")

#Declare lists
deltas = file[:,3]
D = rad(np.mean(deltas))
deltas.tolist()
alphas = file[:,2]
A = rad(np.mean(alphas))
alphas.tolist()

#number of finding stars
N = alphas.size

#Centroids
x_centroids = file[:,0]
print(x_centroids)
y_centroids = file[:,1]
print(y_centroids)

def flatten_alpha(delta,alpha,A,D,L,x_centroid,y_centroid):
    delta = rad(delta)
    alpha = rad(alpha)
    H = sin(delta)*sin(D)+cos(delta)*cos(D)*cos(alpha-A)
    flat_alpha = (cos(delta)*sin(alpha - A))/H - x_centroid/L
    return flat_alpha

def flatten_delta(delta,alpha,A,D,L,x_centroid,y_centroid):
    delta = rad(delta)
    alpha = rad(alpha)
    H = sin(delta)*sin(D)+cos(delta)*cos(D)*cos(alpha-A)
    flat_delta = (sin(delta)*cos(D) - cos(delta)*sin(D)*cos(alpha - A))/H - y_centroid/L
    return flat_delta

flat_alphas = []
for i in range(N):
    flat_alphas.append(flatten_alpha(deltas[i],alphas[i],A,D,L,x_centroid,y_centroid))
np.asarray(flat_alphas)
sum_alphas = np.sum(flat_alphas)
print(flat_alphas)

flat_deltas = []
for i in range(N):
    flat_deltas.append(flatten_delta(deltas[i],alphas[i],A,D,L,x_centroid,y_centroid))
np.asarray(flat_deltas)
sum_deltas = np.sum(flat_deltas)
print(flat_deltas)

def plate_elements(alphas,deltas,x_centroids,y_centroids):
    sum_alpha = np.sum(alphas)
    sum_x = np.sum(x_centroids)
    sum_y = np.sum(y_centroids)
    sum_alpha_x = np.dot(alphas,x_centroids)
    sum_alpha_y = np.dot(alphas,y_centroids)
    sum_x_squared = np.dot(x_centroids,x_centroids)
    sum_y_squared = np.dot(y_centroids,y_centroids)
    sum_x_y = np.dot(x_centroids,y_centroids)

    A = np.array([[N,sum_x,sum_y],[sum_x,sum_x_squared,sum_x_y],[sum_y,sum_x_y,sum_y_squared]])
    ra_vector = np.array([sum_alpha,sum_alpha_x,sum_alpha_y])
    ra_constants = np.dot(np.linalg.inv(A),ra_vector)

    sum_delta = np.sum(deltas)
    sum_delta_x = np.dot(deltas,x_centroids)
    sum_delta_y = np.dot(deltas,y_centroids)
    delta_vector = np.array([sum_delta,sum_delta_x,sum_delta_y])
    dec_constants = np.dot(np.linalg.inv(A),delta_vector)

    return ra_constants[0],dec_constants[0],ra_constants[1],dec_constants[1],ra_constants[2],dec_constants[2]

b1 = plate_elements(flat_alphas,flat_deltas,x_centroids,y_centroids)[0]
b2 = plate_elements(flat_alphas,flat_deltas,x_centroids,y_centroids)[1]
a11 = plate_elements(flat_alphas,flat_deltas,x_centroids,y_centroids)[2]
a21 = plate_elements(flat_alphas,flat_deltas,x_centroids,y_centroids)[3]
a12 = plate_elements(flat_alphas,flat_deltas,x_centroids,y_centroids)[4]
a22 = plate_elements(flat_alphas,flat_deltas,x_centroids,y_centroids)[5]

def ra_and_dec(b1,b2,a11,a21,a12,a22,x_centroid,y_centroid):
    estimated_RA = (b1 + a11 * x_centroid + a12 * y_centroid)/15
    estimated_dec = b2 + a21 * x_centroid + a22 * y_centroid
    return estimated_RA,estimated_dec

alpha_flat = ra_and_dec(b1,b2,a11,a21,a12,a22,x_centroid,y_centroid)[0] + x_centroid / L
delta_flat = ra_and_dec(b1,b2,a11,a21,a12,a22,x_centroid,y_centroid)[1] + y_centroid / L

def unflatten(A,D,alpha_flat,delta_flat):
    Delta = cos(D) - delta_flat * sin(D)
    Gamma = sqrt(alpha_flat**2 + Delta**2)
    alpha = (deg(A) + deg(atan2(alpha_flat,Delta)))/15
    delta = deg(atan2(sin(D) + delta_flat * cos(D),Gamma))
    return alpha,delta

alpha = unflatten(A,D,alpha_flat,delta_flat)[0]
delta = unflatten(A,D,alpha_flat,delta_flat)[1]

def uncertainty(alphas,deltas,b1,b2,a11,a21,a12,a22,x_centroids,y_centroids):
    fitSumAlpha = 0
    for i in range(0,N):
        fitSumAlpha += ((alphas[i]) - (b1) - (a11) * (x_centroids[i]) - (a12) * (y_centroids[i]))**2
    fitSumDelta = 0
    for j in range(0,N):
        fitSumDelta += ((deltas[j]) - (b2) - (a21) * (x_centroids[j]) - (a22) * (y_centroids[j]))**2
    sigma_ra = math.sqrt((1/(N-3))*fitSumAlpha)*3600
    sigma_dec = math.sqrt((1/(N-3))*fitSumDelta)*3600
    return sigma_ra,sigma_dec

sigma_ra = uncertainty(alphas,deltas,b1,b2,a11,a21,a12,a22,x_centroids,y_centroids)[0]
sigma_dec = uncertainty(alphas,deltas,b1,b2,a11,a21,a12,a22,x_centroids,y_centroids)[1]

#estimated RA
alpha_hr = int(alpha)
alpha_min = int((alpha - alpha_hr)*60)
alpha_sec = (alpha - alpha_hr - alpha_min / 60)*3600

#Estimated Dec
delta_deg = int(delta)
delta_min = int((delta - delta_deg)*60)
delta_sec = (delta - delta_deg - delta_min / 60)*3600


#Print everything
print("****************")
print("plate constants")
print("****************")
print("b1:",b1,"deg")
print("b2:",b2,"deg")
print("a11:",a11,"deg/pix")
print("a12:",a12,"deg/pix")
print("a21:",a21,"deg/pix")
print("a22:",a22,"deg/pix")
print("****************")
print("uncertainty")
print("****************")
print("RA:",sigma_ra,"arcsec")
print("Dec:",sigma_dec,"arcsec")
print("****************")
print("astrometry for (x,y) =",(x_centroid,y_centroid))
print("****************")
print("RA =",alpha_hr,"hr",alpha_min,"min",alpha_sec,"sec")
print("Dec =",delta_deg,"hr",delta_min,"min",delta_sec,"sec")
