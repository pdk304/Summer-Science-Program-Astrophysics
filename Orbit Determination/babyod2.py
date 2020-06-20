from math import sqrt, sin, cos, atan2, asin, acos, radians as rad, degrees as deg,pi
import numpy as np

k = .01720209895

epsilon = rad(-23.4352)

#Dummy test values for debugging purposes
r0 = np.array([.13216218,-1.23287955,-.32258829])
r0dot = np.array([.96835101, .13676487, -.15293932])

epsilonrot = np.array([[1,0,0],[0,cos(epsilon),-sin(epsilon)],[0,sin(epsilon),cos(epsilon)]])
r = np.dot(epsilonrot,r0)
rdot = np.dot(epsilonrot,r0dot)

rnorm = np.linalg.norm(r)
rdotnorm = np.linalg.norm(rdot)

#Calculate a
def get_a(rnorm,rdot):
    a = (2*rnorm**-1 - np.dot(rdot,rdot))**-1
    return a

a = get_a(rnorm,rdot)

n = sqrt(k**2 / a**3)

#Calculate e
def get_e(r,rdot,a):
    e = sqrt(1 - (np.linalg.norm(np.cross(r,rdot)))**2 / a)
    return e

e = get_e(r,rdot,a)

#Calculate h
h = k * np.cross(r,rdot)
hnorm = np.linalg.norm(h)

#Calculate I
I = acos(h[2] / hnorm)

#Calculate Omega
def get_Omega(h,I):
    sin_Omega = h[0] / (hnorm * sin(I))
    cos_Omega = - h[1] / (hnorm * sin(I))
    Omega = atan2(sin_Omega,cos_Omega)
    return Omega

Omega = get_Omega(h,I)

#Calculate omegaf
def get_omegaf(r,rnorm,I,Omega):
    sin_omegaf = r[2] / (rnorm * sin(I))
    cos_omegaf = (1 / cos(Omega)) * (r[0]/rnorm + cos(I)*sin_omegaf * sin(Omega))
    omegaf = atan2(sin_omegaf,cos_omegaf)
    return omegaf

omegaf = get_omegaf(r,rnorm,I,Omega)

#Calculate f
def get_f(a,e,r,rdot,rnorm):
    cos_f = (1 / e) * ((a*(1 - e**2))/rnorm - 1)
    sin_f = (np.dot(r,rdot) / (e * rnorm)) * sqrt(a*(1-e**2))
    f = atan2(sin_f,cos_f)
    return f

f = get_f(a,e,r,rdot,rnorm)

#Calculate omega
def get_omega(omegaf,f):
    omega = omegaf - f
    return omega

omega = get_omega(omegaf,f)

#Calculate M
def get_M(a,e,rnorm,f,t_2,t_0):
    n = sqrt(k**2 / a**3)
    E = acos((1/e)*(1 - rnorm / a))
    M_t2 = E - e * sin(E)
    M_t0 = M_t2 + n*(t_0 - t_2)
    return M_t0

#encapsulate all functions into a single orbital elements function
def get_orbital_elements(r,rdot,t_2,t_0):
    rnorm = np.linalg.norm(r)
    rdotnorm = np.linalg.norm(rdot)
    a = get_a(rnorm,rdot)
    n = sqrt(k**2 / a**3)
    e = get_e(r,rdot,a)
    h = k * np.cross(r,rdot)
    hnorm = np.linalg.norm(h)
    I = acos(h[2] / hnorm)
    Omega = get_Omega(h,I)
    omegaf = get_omegaf(r,rnorm,I,Omega)
    f = get_f(a,e,r,rdot,rnorm)
    omega = get_omega(omegaf,f)
    M = get_M(a,e,rnorm,f,t_2,t_0)
    I = deg(I) % 360
    Omega = deg(Omega) % 360
    omega = deg(omega) % 360
    M = deg(M) % 360
    return a,e,I,Omega,omega,M
