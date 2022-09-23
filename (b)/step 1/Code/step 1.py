
'''Calling library'''
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

'''Functions'''
# Susceptible
def f_sus(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd):
    x = - beta[i] * sus_old * inf_old - v_speed_1st[i]
    return x

# Vaccination1
def f_vac1(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd):
    x = v_speed_1st[i] - v_speed_2nd[i] - beta[i] * vac1_old * inf_old * (1 - effi_1st)
    return x

# Vaccination2
def f_vac2(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd):
    x = v_speed_2nd[i] - beta[i] * vac2_old * inf_old * (1 - effi_2nd)
    return x

# Exposed
def f_exp(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd):
    x = beta[i] * sus_old * inf_old + beta[i] * vac1_old * inf_old * (1 - effi_1st) + beta[i] * vac2_old * inf_old * (1 - effi_2nd) - alpha * exp_old
    return x

# Infectious
def f_inf(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd): 
    x = alpha * exp_old - gamma * inf_old
    return x

# Recovered + Death
def f_rem(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd):
    x = gamma * inf_old
    return x

'''Main'''

# Upload csv files
Spreading_data = pd.read_csv("Spreading_data.csv")

# Time and the number of days
days = len(Spreading_data) - 2
days_divided = days * 10
dt = 1.0 / 10.0 # [day]

# Matrix for output
sus = np.zeros(days_divided+1)
vac1 = np.zeros(days_divided+1)
vac2 = np.zeros(days_divided+1)
exp = np.zeros(days_divided+1)
inf = np.zeros(days_divided+1)
inf_new = np.zeros(days_divided+1)
rem = np.zeros(days_divided+1)

sus_cal = np.zeros(days+1)
vac1_cal = np.zeros(days+1)
vac2_cal = np.zeros(days+1)
exp_cal = np.zeros(days+1)
inf_cal = np.zeros(days+1)
rem_cal = np.zeros(days+1)
inf_tot_cal = np.zeros(days+1)
inf_new_cal = np.zeros(days+1)

# Initial conditions about SEIR
n_total = 1.263 * 1.0e+8
exp_ini = 0.1786

inf_ini = Spreading_data.iloc[0]['I_now']
rem_ini = Spreading_data.iloc[0]['R']
vac1_ini = Spreading_data.iloc[0]['v_speed_1st']
vac2_ini = Spreading_data.iloc[0]['v_speed_2nd']
sus_ini = n_total - exp_ini - inf_ini - rem_ini - vac1_ini - vac2_ini

sus[0] = sus_ini
vac1[0] = vac1_ini
vac2[0] = vac2_ini
exp[0] = exp_ini
inf[0] = inf_ini
rem[0] = rem_ini

# Temporary recording
sus_old = sus_ini
vac1_old = vac1_ini
vac2_old = vac2_ini
exp_old = exp_ini
inf_old = inf_ini
rem_old = rem_ini

# Nature of Disease
alpha = 0.4
gamma = 0.0691712976652005
beta = np.zeros(days_divided)
beta_daily = np.zeros(days)

# Conditions of Vaccine
effi_1st = 0.53
effi_2nd = 0.87
v_speed_1st_daily = Spreading_data.iloc[0:]['v_speed_1st'] 
v_speed_2nd_daily = Spreading_data.iloc[0:]['v_speed_2nd']

v_speed_1st = np.zeros(days_divided)
v_speed_2nd = np.zeros(days_divided)

# Matrix for provisional E
exp_fm_exInew = np.zeros(days+1)
exp_fm_exInew[0] = 0.0

# Counter
counter_day = 0
k = 0

# Numerical solution: Runge-Kutta method
for i in range(days_divided):
    if i % int(1.0 / dt) == 0:
        exp_fm_exInew[counter_day+1] = (Spreading_data.iloc[counter_day + 1]['I_new'] + Spreading_data.iloc[counter_day + 2]['I_new'])/(2*alpha)
        dE = exp_fm_exInew[counter_day+1] - exp_old
        beta_daily[counter_day] = (dE + alpha * exp_old)/(inf_old * (sus_old + vac1_old * effi_1st + vac2_old * effi_2nd))
        counter_day += 1

    beta[i] = beta_daily[counter_day-1]
    v_speed_1st[i] = v_speed_1st_daily[counter_day-1]
    v_speed_2nd[i] = v_speed_2nd_daily[counter_day-1]
    
    k_sus_1 = f_sus(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_vac1_1 = f_vac1(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_vac2_1 = f_vac2(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_exp_1 = f_exp(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_inf_1 = f_inf(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_rem_1 = f_rem(sus_old, exp_old, inf_old, vac1_old, vac2_old, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    
    k_sus_2 = f_sus(sus_old + k_sus_1 * 0.5 * dt, exp_old + k_exp_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, vac1_old + k_vac1_1 * 0.5 * dt, vac2_old + k_vac2_1 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_vac1_2 = f_vac1(sus_old + k_sus_1 * 0.5 * dt, exp_old + k_exp_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, vac1_old + k_vac1_1 * 0.5 * dt, vac2_old + k_vac2_1 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_vac2_2 = f_vac2(sus_old + k_sus_1 * 0.5 * dt, exp_old + k_exp_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, vac1_old + k_vac1_1 * 0.5 * dt, vac2_old + k_vac2_1 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_exp_2 = f_exp(sus_old + k_sus_1 * 0.5 * dt, exp_old + k_exp_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, vac1_old + k_vac1_1 * 0.5 * dt, vac2_old + k_vac2_1 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_inf_2 = f_inf(sus_old + k_sus_1 * 0.5 * dt, exp_old + k_exp_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, vac1_old + k_vac1_1 * 0.5 * dt, vac2_old + k_vac2_1 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_rem_2 = f_rem(sus_old + k_sus_1 * 0.5 * dt, exp_old + k_exp_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, vac1_old + k_vac1_1 * 0.5 * dt, vac2_old + k_vac2_1 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    
    k_sus_3 = f_sus(sus_old + k_sus_2 * 0.5 * dt, exp_old + k_exp_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, vac1_old + k_vac1_2 * 0.5 * dt, vac2_old + k_vac2_2 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_vac1_3 = f_vac1(sus_old + k_sus_2 * 0.5 * dt, exp_old + k_exp_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, vac1_old + k_vac1_2 * 0.5 * dt, vac2_old + k_vac2_2 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_vac2_3 = f_vac2(sus_old + k_sus_2 * 0.5 * dt, exp_old + k_exp_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, vac1_old + k_vac1_2 * 0.5 * dt, vac2_old + k_vac2_2 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_exp_3 = f_exp(sus_old + k_sus_2 * 0.5 * dt, exp_old + k_exp_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, vac1_old + k_vac1_2 * 0.5 * dt, vac2_old + k_vac2_2 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_inf_3 = f_inf(sus_old + k_sus_2 * 0.5 * dt, exp_old + k_exp_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, vac1_old + k_vac1_2 * 0.5 * dt, vac2_old + k_vac2_2 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_rem_3 = f_rem(sus_old + k_sus_2 * 0.5 * dt, exp_old + k_exp_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, vac1_old + k_vac1_2 * 0.5 * dt, vac2_old + k_vac2_2 * 0.5 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
            
    k_sus_4 = f_sus(sus_old + k_sus_3 * dt, exp_old + k_exp_3 * dt, inf_old + k_inf_3 * dt, vac1_old + k_vac1_3 * dt, vac2_old + k_vac2_3 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_vac1_4 = f_vac1(sus_old + k_sus_3 * dt, exp_old + k_exp_3 * dt, inf_old + k_inf_3 * dt, vac1_old + k_vac1_3 * dt, vac2_old + k_vac2_3 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_vac2_4 = f_vac2(sus_old + k_sus_3 * dt, exp_old + k_exp_3 * dt, inf_old + k_inf_3 * dt, vac1_old + k_vac1_3 * dt, vac2_old + k_vac2_3 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_exp_4 = f_exp(sus_old + k_sus_3 * dt, exp_old + k_exp_3 * dt, inf_old + k_inf_3 * dt, vac1_old + k_vac1_3 * dt, vac2_old + k_vac2_3 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_inf_4 = f_inf(sus_old + k_sus_3 * dt, exp_old + k_exp_3 * dt, inf_old + k_inf_3 * dt, vac1_old + k_vac1_3 * dt, vac2_old + k_vac2_3 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
    k_rem_4 = f_rem(sus_old + k_sus_3 * dt, exp_old + k_exp_3 * dt, inf_old + k_inf_3 * dt, vac1_old + k_vac1_3 * dt, vac2_old + k_vac2_3 * dt, beta, alpha, gamma, v_speed_1st, v_speed_2nd, effi_1st, effi_2nd)
            
    sus_new = sus_old + dt / 6.0 * (k_sus_1 + 2.0 * k_sus_2 + 2.0 * k_sus_3 + k_sus_4)
    vac1_new = vac1_old + dt / 6.0 * (k_vac1_1 + 2.0 * k_vac1_2 + 2.0 * k_vac1_3 + k_vac1_4)
    vac2_new = vac2_old + dt / 6.0 * (k_vac2_1 + 2.0 * k_vac2_2 + 2.0 * k_vac2_3 + k_vac2_4)
    exp_new = exp_old + dt / 6.0 * (k_exp_1 + 2.0 * k_exp_2 + 2.0 * k_exp_3 + k_exp_4)
    inf_new = inf_old + dt / 6.0 * (k_inf_1 + 2.0 * k_inf_2 + 2.0 * k_inf_3 + k_inf_4)
    rem_new = rem_old + dt / 6.0 * (k_rem_1 + 2.0 * k_rem_2 + 2.0 * k_rem_3 + k_rem_4)
    
    sus_old = sus_new
    vac1_old = vac1_new
    vac2_old = vac2_new
    exp_old = exp_new
    inf_old = inf_new
    rem_old = rem_new
            
    sus[i+1] = sus_old
    vac1[i+1] = vac1_old
    vac2[i+1] = vac2_old
    exp[i+1] = exp_old
    inf[i+1] = inf_old
    rem[i+1] = rem_old 

    i += 1
            
for j in range(0, days_divided+1):
    if j % int(1.0 / dt) == 0:
        sus_cal[k] = sus[j] 
        vac1_cal[k] = vac1[j]
        vac2_cal[k] = vac2[j]
        exp_cal[k] = exp[j] 
        inf_cal[k] = inf[j] 
        rem_cal[k] = rem[j] 
        inf_tot_cal[k] = (inf[j] + rem[j]) 
        k += 1

inf_new_cal[0] = inf_tot_cal[0]
for ij in range(1, days):
    inf_new_cal[ij] = inf_tot_cal[ij] - inf_tot_cal[ij - 1]


beta_data = pd.DataFrame(beta_daily * n_total, columns= ['beta'])
date_data = pd.DataFrame(Spreading_data.iloc[0:len(Spreading_data) - 2]['date'])

R2_I = r2_score(Spreading_data.iloc[0:days + 1]['I_now'], inf_cal)
R2_I_new = r2_score(Spreading_data.iloc[0:days + 1]['I_new'], inf_new_cal)
R2_R = r2_score(Spreading_data.iloc[0:days + 1]['R'], rem_cal)
print(R2_I)
print(R2_I_new)
print(R2_R)

outputdata = pd.concat([date_data, beta_data], axis=1)
outputdata.to_csv("beta_fm_exIandR.csv", index=False)

