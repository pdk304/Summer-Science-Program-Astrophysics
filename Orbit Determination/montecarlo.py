import numpy as np
import pandas as pd
import matplotlib as plt
from mogfunction import *
from functions import *
from babyod2 import *

file = np.loadtxt("2002QF15moginputs.txt")

years,months,days = file[:,0],file[:,1],file[:,2]
hours,minutes,seconds = file[:,3],file[:,4],file[:,5]
alpha_hours,alpha_minutes,alpha_seconds = file[:,6],file[:,7],file[:,8]
delta_deg,delta_minutes,delta_seconds = file[:,9],file[:,10],file[:,11]
R_x,R_y,R_z = file[:,12],file[:,13],file[:,14]

t = JDf(years,months,days,hours,minutes,seconds)
t_0 = t_0 = JDmini(2018,7,22,6,0,0)

alphas = get_alphasf(alpha_hours,alpha_minutes,alpha_seconds)
deltas = get_deltasf(delta_deg,delta_minutes,delta_seconds)

sigma_alpha_arcsecs = 0.890,1.514,1.30
sigma_delta_arcsecs = 0.539,0.5530,.56

sigma_alphas = arcsec_to_rad(sigma_alpha_arcsecs)
sigma_deltas = arcsec_to_rad(sigma_delta_arcsecs)

alpha_1 = alphas[0]
alpha_2 = alphas[1]
alpha_3 = alphas[2]

sigma_alpha_1 = sigma_alphas[0]
sigma_alpha_2 = sigma_alphas[1]
sigma_alpha_3 = sigma_alphas[2]

delta_1 = deltas[0]
delta_2 = deltas[1]
delta_3 = deltas[2]

sigma_delta_1 = sigma_deltas[0]
sigma_delta_2 = sigma_deltas[1]
sigma_delta_3 = sigma_deltas[2]

N = 10000

sample_alphas_1 = np.random.normal(alpha_1,sigma_alpha_1,N)
sample_alphas_2 = np.random.normal(alpha_2,sigma_alpha_2,N)
sample_alphas_3 = np.random.normal(alpha_3,sigma_alpha_3,N)
sample_alphas = []
for i in range(N):
    sample_alphas.append((sample_alphas_1[i],sample_alphas_2[i],sample_alphas_3[i]))

print(sample_alphas)

sample_deltas_1 = np.random.normal(delta_1,sigma_delta_1,N)
sample_deltas_2 = np.random.normal(delta_2,sigma_delta_2,N)
sample_deltas_3 = np.random.normal(delta_3,sigma_delta_3,N)
sample_deltas = []

for i in range(N):
    sample_deltas.append((sample_deltas_1[i],sample_deltas_2[i],sample_deltas_3[i]))

print(sample_deltas)

list = []
[ list.append(mog(sample_alphas[i],sample_deltas[i],t,R_x,R_y,R_z,t_0)) for i in range(N) ]
print(list)

sample_elements_arr = np.asarray(list)
sample_elements_df = pd.DataFrame(sample_elements_arr,columns = ['a','e','I','Omega','omega','M'])
sample_elements_df.to_csv('monte data 10000.csv',sep = ' ')
