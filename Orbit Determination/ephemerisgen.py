#Ephemeris generation

import numpy as np
from math import sin,cos,sqrt,radians as rad,degrees as deg,acos,asin,atan2
from time import sleep

#Step 1
e = 0.4516752702182729
a = 1.703044972857831 #AU
i = rad(4.644780120140178)
Omega = rad(144.4119981532961 )
omega = rad(197.831071215316)
M_0 = rad(282.0228443667625)

k = .01720209894
mu = k**2

T = 2456800.5 #JHD
t = 2456842.5 #JHD

n = sqrt(mu / a**3)
M = n*(t - T) + M_0

#Step 2
def solve_kep(M):
    E = M
    Mguess = E - (e*sin(E))
    while abs(Mguess - M) > 10e-6:
        f = M - (E - e*sin(E))
        fprime = (e*cos(E)) - 1
        E = E - (f / fprime)
        Mguess = E - (e*sin(E))
    return E

E = solve_kep(M)
print("E:",E)

#Step 3
r = np.array([a*cos(E)-a*e,a*sqrt(1-e**2)*sin(E),0])
print("r:",r)

#Step 4
omegarot = np.array([[cos(omega),-sin(omega),0],[sin(omega),cos(omega),0],[0,0,1]])
irot = np.array([[1,0,0],[0,cos(i),-sin(i)],[0,sin(i),cos(i)]])
Omegarot = np.array([[cos(Omega),-sin(Omega),0],[sin(Omega),cos(Omega),0],[0,0,1]])

rot = np.dot(Omegarot,np.dot(irot,omegarot))
print("rot:",rot)
r1 = np.dot(rot,r)
print("r1:",r1)

#Step 5
epsilon = rad(23.43688)
epsilonrot = np.array([[1,0,0],[0,cos(epsilon),-sin(epsilon)],[0,sin(epsilon),cos(epsilon)]])
r2 = np.dot(epsilonrot,r1)

#Step 6
R = np.array([-2.604100023131745e-01,9.827797286991941e-01,-4.385819099838496e-5])
print("R:",R)

#Step 7
rho = r1 + R
rho = np.dot(epsilonrot,rho)
print("rho:",rho)
norm_rho = sqrt(rho[0]**2 + rho[1]**2 + rho[2]**2)
rho_hat_0 = rho[0] / norm_rho
rho_hat_1 = rho[1] / norm_rho
rho_hat_2 = rho[2] / norm_rho
rho_hat = np.array([rho_hat_0,rho_hat_1,rho_hat_2])

delta = asin(rho_hat[2])

def get_alpha(vec):
    sin_alpha = vec[1]
    cos_alpha = vec[0] / cos(delta)
    alpha = atan2(sin_alpha,cos_alpha)
    return alpha

alpha = deg(get_alpha(rho_hat))/15
alpha_hr = int(alpha)
alpha_min = int((alpha - alpha_hr)*60)
alpha_sec = (alpha - alpha_hr - alpha_min / 60)*3600

delta = deg(delta)
delta_deg = int(delta)
delta_min = int((delta - delta_deg)*60)
delta_sec = (delta - delta_deg - delta_min / 60)*3600

#Print
print("delta =",delta_deg,":",delta_min,":",delta_sec)
print("alpha =",alpha_hr,":",alpha_min,":",alpha_sec)
