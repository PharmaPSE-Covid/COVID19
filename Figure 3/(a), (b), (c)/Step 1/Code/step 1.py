
'''Calling library'''
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

'''Functions'''
def f_s_e(sus_old, inf_old, inf_v1_old, inf_v2_old, beta):
    x = beta[i] * sus_old * (inf_old + inf_v1_old + inf_v2_old)
    return x

def f_v1_e1(vac1_old, inf_old, inf_v1_old, inf_v2_old, beta, effi_1st):
    x = beta[i] * vac1_old * (inf_old + inf_v1_old + inf_v2_old) * (1 - effi_1st)
    return x

def f_v2_e2(vac2_old, inf_old, inf_v1_old, inf_v2_old, beta, effi_2nd):
    x = beta[i] * vac2_old * (inf_old + inf_v1_old + inf_v2_old) * (1 - effi_2nd)
    return x

def f_e_i(exp_old, alpha):
    x = alpha * exp_old
    return x

def f_e1_i1(exp_v1_old, alpha_v1):
    x = alpha_v1 * exp_v1_old
    return x

def f_e2_i2(exp_v2_old, alpha_v2):
    x = alpha_v2 * exp_v2_old
    return x

def f_i_r(inf_old, gamma):
    x = gamma * inf_old
    return x

def f_i1_r(inf_v1_old, gamma_v1):
    x = gamma_v1 * inf_v1_old
    return x

def f_i2_r(inf_v2_old, gamma_v2):
    x = gamma_v2 * inf_v2_old
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
exp_v1 = np.zeros(days_divided+1)
exp_v2 = np.zeros(days_divided+1)
inf = np.zeros(days_divided+1)
inf_v1 = np.zeros(days_divided+1)
inf_v2 = np.zeros(days_divided+1)
rem = np.zeros(days_divided+1)

v1_to_e1 = np.zeros(days_divided)
v2_to_e2 = np.zeros(days_divided)
e_to_i = np.zeros(days_divided)
e1_to_i1 = np.zeros(days_divided)
e2_to_i2 = np.zeros(days_divided)
i_to_r = np.zeros(days_divided)
i1_to_r = np.zeros(days_divided)
i2_to_r = np.zeros(days_divided)


sus_cal = np.zeros(days+1)
vac1_cal = np.zeros(days+1)
vac2_cal = np.zeros(days+1)
exp_cal = np.zeros(days+1)
exp_v1_cal = np.zeros(days+1)
exp_v2_cal = np.zeros(days+1)
inf_cal = np.zeros(days+1)
inf_v1_cal = np.zeros(days+1)
inf_v2_cal = np.zeros(days+1)
rem_cal = np.zeros(days+1)

inf_new_cal = np.zeros(days+1)
inf_v1_new_cal = np.zeros(days+1)
inf_v2_new_cal = np.zeros(days+1)

# Initial conditions about SEIR
n_total = 1.263 * 1.0e+8
exp_ini = 0.1786
exp_v1_ini = 0.0000
exp_v2_ini = 0.0000

inf_ini = Spreading_data.iloc[0]['I_now']
inf_v1_ini = 0.0000
inf_v2_ini = 0.0000
rem_ini = Spreading_data.iloc[0]['Rem']
vac1_ini = Spreading_data.iloc[0]['v_speed_1st']
vac2_ini = Spreading_data.iloc[0]['v_speed_2nd']
sus_ini = n_total - vac1_ini - vac2_ini - exp_ini - exp_v1_ini- exp_v2_ini - inf_ini - inf_v1_ini - inf_v2_ini - rem_ini

sus[0] = sus_ini
vac1[0] = vac1_ini
vac2[0] = vac2_ini
exp[0] = exp_ini
exp_v1[0] = exp_v1_ini
exp_v2[0] = exp_v2_ini
inf[0] = inf_ini
inf_v1[0] = inf_v1_ini
inf_v2[0] = inf_v2_ini
rem[0] = rem_ini

# Temporary recording
sus_old = sus_ini
vac1_old = vac1_ini
vac2_old = vac2_ini
exp_old = exp_ini
exp_v1_old = exp_v1_ini
exp_v2_old = exp_v2_ini
inf_old = inf_ini
inf_v1_old = inf_v1_ini
inf_v2_old = inf_v2_ini
rem_old = rem_ini

# Nature of Disease
alpha = 0.4
alpha_v1 = 0.4
alpha_v2 = 0.4
gamma = 0.0691712976652005
gamma_v1 = 0.0691712976652005
gamma_v2 = 0.0691712976652005
beta = np.zeros(days_divided)
beta_daily = np.zeros(days)

anti_loss_r = np.zeros(days_divided)
anti_loss_r_daily = np.zeros(days)
anti_loss_r_period = 150 #days
anti_loss_v1 = np.zeros(days_divided)
anti_loss_v1_daily = np.zeros(days)
anti_loss_v1_period = 150 #days
anti_loss_v2 = np.zeros(days_divided)
anti_loss_v2_daily = np.zeros(days)
anti_loss_v2_period = 150 #days

# Conditions of Vaccine
effi_1st = 0.53
effi_2nd = 0.87
sev_rate = 0.01
sev_rate_v1 = 0.001
sev_rate_v2 = 0.001
v_speed_1st_daily = Spreading_data.iloc[0:]['v_speed_1st'] 
v_speed_2nd_daily = Spreading_data.iloc[0:]['v_speed_2nd']

v_speed_1st = np.zeros(days_divided)
v_speed_2nd = np.zeros(days_divided)

# Matrix for provisional E
exp_fm_exInew = np.zeros(days+1)
exp_fm_exInew[0] = 0.0

# Counter
counter_day = 0

# Numerical solution: Runge-Kutta method
for i in range(days_divided):
    if i % int(1.0 / dt) == 0:
        exp_fm_exInew[counter_day+1] = (Spreading_data.iloc[counter_day + 1]['I_new'] + Spreading_data.iloc[counter_day + 2]['I_new'])/(2*alpha)
        dE = exp_fm_exInew[counter_day+1] - (exp_old + exp_v1_old + exp_v2_old)
        beta_daily[counter_day] = (dE + alpha * exp_old + alpha_v1 * exp_v1_old + alpha_v2 * exp_v2_old)/((inf_old + inf_v1_old + inf_v2_old) * (sus_old + vac1_old * (1 - effi_1st) + vac2_old * (1 - effi_2nd)))
        
        if i > 0 and counter_day + anti_loss_r_period < days:
            i_r_sum = 0
            for j in range(10):
                i_r_sum += i_to_r[i+j-10] + i1_to_r[i+j-10] + i2_to_r[i+j-10]
            anti_loss_r_daily[counter_day + anti_loss_r_period] = i_r_sum
        if i > 0 and counter_day + anti_loss_v1_period < days:
            v1_e1_sum = 0
            for j in range(10):
                v1_e1_sum += v1_to_e1[i+j-10]
            anti_loss_v1_daily[counter_day + anti_loss_v1_period] = v_speed_1st_daily[counter_day] * 0.003 - v1_e1_sum
        if i > 0 and counter_day + anti_loss_v2_period < days:
            v2_e2_sum = 0
            for j in range(10):
                v2_e2_sum += v2_to_e2[i+j-10]
            anti_loss_v2_daily[counter_day + anti_loss_v2_period] = v_speed_2nd_daily[counter_day] - v2_e2_sum
        
        counter_day += 1        

    beta[i] = beta_daily[counter_day-1]
    v_speed_1st[i] = v_speed_1st_daily[counter_day-1]
    v_speed_2nd[i] = v_speed_2nd_daily[counter_day-1]
    anti_loss_r[i] = anti_loss_r_daily[counter_day-1]
    anti_loss_v1[i] = anti_loss_v1_daily[counter_day-1]
    anti_loss_v2[i] = anti_loss_v2_daily[counter_day-1]
    

    k_s_e_1 = f_s_e(sus_old, inf_old, inf_v1_old, inf_v2_old, beta)
    k_v1_e1_1 = f_v1_e1(vac1_old, inf_old, inf_v1_old, inf_v2_old, beta, effi_1st)
    k_v2_e2_1 = f_v2_e2(vac2_old, inf_old, inf_v1_old, inf_v2_old, beta, effi_2nd)
    k_e_i_1 = f_e_i(exp_old, alpha)
    k_e1_i1_1 = f_e1_i1(exp_v1_old, alpha_v1)
    k_e2_i2_1 = f_e2_i2(exp_v2_old, alpha_v2)
    k_i_r_1 = f_i_r(inf_old, gamma)
    k_i1_r_1 = f_i1_r(inf_v1_old, gamma_v1)
    k_i2_r_1 = f_i2_r(inf_v2_old, gamma_v2)

    k_sus_1 = - k_s_e_1 - v_speed_1st[i] + anti_loss_r[i] + anti_loss_v1[i] + anti_loss_v2[i]
    k_vac1_1 = v_speed_1st[i] - v_speed_2nd[i] - k_v1_e1_1 - anti_loss_v1[i]
    k_vac2_1 = v_speed_2nd[i] - k_v2_e2_1 - anti_loss_v2[i]
    k_exp_1 = k_s_e_1 - k_e_i_1
    k_exp_v1_1 = k_v1_e1_1 - k_e1_i1_1
    k_exp_v2_1 = k_v2_e2_1 - k_e2_i2_1
    k_inf_1 = k_e_i_1 - k_i_r_1
    k_inf_v1_1 = k_e1_i1_1 - k_i1_r_1
    k_inf_v2_1 = k_e2_i2_1 - k_i2_r_1
    k_rem_1 = k_i_r_1 + k_i1_r_1 + k_i2_r_1 - anti_loss_r[i]


    k_s_e_2 = f_s_e(sus_old + k_sus_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, inf_v1_old + k_inf_v1_1 * 0.5 * dt, inf_v2_old + k_inf_v2_1 * 0.5 * dt, beta)
    k_v1_e1_2 = f_v1_e1(vac1_old + k_vac1_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, inf_v1_old + k_inf_v1_1 * 0.5 * dt, inf_v2_old + k_inf_v2_1 * 0.5 * dt, beta, effi_1st)
    k_v2_e2_2 = f_v2_e2(vac2_old + k_vac2_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, inf_v1_old + k_inf_v1_1 * 0.5 * dt, inf_v2_old + k_inf_v2_1 * 0.5 * dt, beta, effi_2nd)
    k_e_i_2 = f_e_i(exp_old + k_exp_1 * 0.5 * dt, alpha)
    k_e1_i1_2 = f_e1_i1(exp_v1_old + k_exp_v1_1 * 0.5 * dt, alpha_v1)
    k_e2_i2_2 = f_e2_i2(exp_v2_old + k_exp_v2_1 * 0.5 * dt, alpha_v2)
    k_i_r_2 = f_i_r(inf_old + k_inf_1 * 0.5 * dt, gamma)
    k_i1_r_2 = f_i1_r(inf_v1_old + k_inf_v1_1 * 0.5 * dt, gamma_v1)
    k_i2_r_2 = f_i2_r(inf_v2_old + k_inf_v2_1 * 0.5 * dt, gamma_v2)
    
    k_sus_2 = - k_s_e_2 - v_speed_1st[i] + anti_loss_r[i] + anti_loss_v1[i] + anti_loss_v2[i]
    k_vac1_2 = v_speed_1st[i] - v_speed_2nd[i] - k_v1_e1_2 - anti_loss_v1[i]
    k_vac2_2 = v_speed_2nd[i] - k_v2_e2_2 - anti_loss_v2[i]
    k_exp_2 = k_s_e_2 - k_e_i_2
    k_exp_v1_2 = k_v1_e1_2 - k_e1_i1_2
    k_exp_v2_2 = k_v2_e2_2 - k_e2_i2_2
    k_inf_2 = k_e_i_2 - k_i_r_2
    k_inf_v1_2 = k_e1_i1_2 - k_i1_r_2
    k_inf_v2_2 = k_e2_i2_2 - k_i2_r_2
    k_rem_2 = k_i_r_2 + k_i1_r_2 + k_i2_r_2 - anti_loss_r[i]


    k_s_e_3 = f_s_e(sus_old + k_sus_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, inf_v1_old + k_inf_v1_2 * 0.5 * dt, inf_v2_old + k_inf_v2_2 * 0.5 * dt, beta)
    k_v1_e1_3 = f_v1_e1(vac1_old + k_vac1_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, inf_v1_old + k_inf_v1_2 * 0.5 * dt, inf_v2_old + k_inf_v2_2 * 0.5 * dt, beta, effi_1st)
    k_v2_e2_3 = f_v2_e2(vac2_old + k_vac2_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, inf_v1_old + k_inf_v1_2 * 0.5 * dt, inf_v2_old + k_inf_v2_2 * 0.5 * dt, beta, effi_2nd)
    k_e_i_3 = f_e_i(exp_old + k_exp_2 * 0.5 * dt, alpha)
    k_e1_i1_3 = f_e1_i1(exp_v1_old + k_exp_v1_2 * 0.5 * dt, alpha_v1)
    k_e2_i2_3 = f_e2_i2(exp_v2_old + k_exp_v2_2 * 0.5 * dt, alpha_v2)
    k_i_r_3 = f_i_r(inf_old + k_inf_2 * 0.5 * dt, gamma)
    k_i1_r_3 = f_i1_r(inf_v1_old + k_inf_v1_2 * 0.5 * dt, gamma_v1)
    k_i2_r_3 = f_i2_r(inf_v2_old + k_inf_v2_2 * 0.5 * dt, gamma_v2)
    
    k_sus_3 = - k_s_e_3 - v_speed_1st[i] + anti_loss_r[i] + anti_loss_v1[i] + anti_loss_v2[i]
    k_vac1_3 = v_speed_1st[i] - v_speed_2nd[i] - k_v1_e1_3 - anti_loss_v1[i]
    k_vac2_3 = v_speed_2nd[i] - k_v2_e2_3 - anti_loss_v2[i]
    k_exp_3 = k_s_e_3 - k_e_i_3
    k_exp_v1_3 = k_v1_e1_3 - k_e1_i1_3
    k_exp_v2_3 = k_v2_e2_3 - k_e2_i2_3
    k_inf_3 = k_e_i_3 - k_i_r_3
    k_inf_v1_3 = k_e1_i1_3 - k_i1_r_3
    k_inf_v2_3 = k_e2_i2_3 - k_i2_r_3
    k_rem_3 = k_i_r_3 + k_i1_r_3 + k_i2_r_3 - anti_loss_r[i]


    k_s_e_4 = f_s_e(sus_old + k_sus_3 * dt, inf_old + k_inf_3 * dt, inf_v1_old + k_inf_v1_3 * dt, inf_v2_old + k_inf_v2_3 * dt, beta)
    k_v1_e1_4 = f_v1_e1(vac1_old + k_vac1_3 * dt, inf_old + k_inf_3 * dt, inf_v1_old + k_inf_v1_3 * dt, inf_v2_old + k_inf_v2_3 * dt, beta, effi_1st)
    k_v2_e2_4 = f_v2_e2(vac2_old + k_vac2_3 * dt, inf_old + k_inf_3 * dt, inf_v1_old + k_inf_v1_3 * dt, inf_v2_old + k_inf_v2_3 * dt, beta, effi_2nd)
    k_e_i_4 = f_e_i(exp_old + k_exp_3 * dt, alpha)
    k_e1_i1_4 = f_e1_i1(exp_v1_old + k_exp_v1_3 * dt, alpha_v1)
    k_e2_i2_4 = f_e2_i2(exp_v2_old + k_exp_v2_3 * dt, alpha_v2)
    k_i_r_4 = f_i_r(inf_old + k_inf_3 * dt, gamma)
    k_i1_r_4 = f_i1_r(inf_v1_old + k_inf_v1_3 * dt, gamma_v1)
    k_i2_r_4 = f_i2_r(inf_v2_old + k_inf_v2_3 * dt, gamma_v2)
    
    k_sus_4 = - k_s_e_4 - v_speed_1st[i] + anti_loss_r[i] + anti_loss_v1[i] + anti_loss_v2[i]
    k_vac1_4 = v_speed_1st[i] - v_speed_2nd[i] - k_v1_e1_4 - anti_loss_v1[i]
    k_vac2_4 = v_speed_2nd[i] - k_v2_e2_4 - anti_loss_v2[i]
    k_exp_4 = k_s_e_4 - k_e_i_4
    k_exp_v1_4 = k_v1_e1_4 - k_e1_i1_4
    k_exp_v2_4 = k_v2_e2_4 - k_e2_i2_4
    k_inf_4 = k_e_i_4 - k_i_r_4
    k_inf_v1_4 = k_e1_i1_4 - k_i1_r_4
    k_inf_v2_4 = k_e2_i2_4 - k_i2_r_4
    k_rem_4 = k_i_r_4 + k_i1_r_4 + k_i2_r_4 - anti_loss_r[i]


    sus_new = sus_old + dt / 6.0 * (k_sus_1 + 2.0 * k_sus_2 + 2.0 * k_sus_3 + k_sus_4)
    vac1_new = vac1_old + dt / 6.0 * (k_vac1_1 + 2.0 * k_vac1_2 + 2.0 * k_vac1_3 + k_vac1_4)
    vac2_new = vac2_old + dt / 6.0 * (k_vac2_1 + 2.0 * k_vac2_2 + 2.0 * k_vac2_3 + k_vac2_4)
    exp_new = exp_old + dt / 6.0 * (k_exp_1 + 2.0 * k_exp_2 + 2.0 * k_exp_3 + k_exp_4)
    exp_v1_new = exp_v1_old + dt / 6.0 * (k_exp_v1_1 + 2.0 * k_exp_v1_2 + 2.0 * k_exp_v1_3 + k_exp_v1_4)
    exp_v2_new = exp_v2_old + dt / 6.0 * (k_exp_v2_1 + 2.0 * k_exp_v2_2 + 2.0 * k_exp_v2_3 + k_exp_v2_4)
    inf_new = inf_old + dt / 6.0 * (k_inf_1 + 2.0 * k_inf_2 + 2.0 * k_inf_3 + k_inf_4)
    inf_v1_new = inf_v1_old + dt / 6.0 * (k_inf_v1_1 + 2.0 * k_inf_v1_2 + 2.0 * k_inf_v1_3 + k_inf_v1_4)
    inf_v2_new = inf_v2_old + dt / 6.0 * (k_inf_v2_1 + 2.0 * k_inf_v2_2 + 2.0 * k_inf_v2_3 + k_inf_v2_4)
    rem_new = rem_old + dt / 6.0 * (k_rem_1 + 2.0 * k_rem_2 + 2.0 * k_rem_3 + k_rem_4)
    
    
    if sus_new < 0:
        sus_new = 0.0
    if vac1_new < 0:
        vac1_new = 0.0
    if vac2_new < 0:
        vac2_new = 0.0
    if exp_new < 0:
        exp_new = 0.0
    if exp_v1_new < 0:
        exp_v1_new = 0.0
    if exp_v2_new < 0:
        exp_v2_new = 0.0
    if inf_new < 0:
        inf_new = 0.0
    if inf_v1_new < 0:
        inf_v1_new = 0.0
    if inf_v2_new < 0:
        inf_v2_new = 0.0
    if rem_new < 0:
        rem_new = 0.0
    
    sus_old = sus_new
    vac1_old = vac1_new
    vac2_old = vac2_new
    exp_old = exp_new
    exp_v1_old = exp_v1_new
    exp_v2_old = exp_v2_new
    inf_old = inf_new
    inf_v1_old = inf_v1_new
    inf_v2_old = inf_v2_new
    rem_old = rem_new
    
    sus[i+1] = sus_old
    vac1[i+1] = vac1_old
    vac2[i+1] = vac2_old
    exp[i+1] = exp_old
    exp_v1[i+1] = exp_v1_old
    exp_v2[i+1] = exp_v2_old
    inf[i+1] = inf_old
    inf_v1[i+1] = inf_v1_old
    inf_v2[i+1] = inf_v2_old
    rem[i+1] = rem_old 

    v1_to_e1[i] = dt / 6.0 * (k_v1_e1_1 + 2.0 * k_v1_e1_2 + 2.0 * k_v1_e1_3 + k_v1_e1_4)
    v2_to_e2[i] = dt / 6.0 * (k_v2_e2_1 + 2.0 * k_v2_e2_2 + 2.0 * k_v2_e2_3 + k_v2_e2_4)
    e_to_i[i] = dt / 6.0 * (k_e_i_1 + 2.0 * k_e_i_2 + 2.0 * k_e_i_3 + k_e_i_4)
    e1_to_i1[i] = dt / 6.0 * (k_e1_i1_1 + 2.0 * k_e1_i1_2 + 2.0 * k_e1_i1_3 + k_e1_i1_4)
    e2_to_i2[i] = dt / 6.0 * (k_e2_i2_1 + 2.0 * k_e2_i2_2 + 2.0 * k_e2_i2_3 + k_e2_i2_4)
    i_to_r[i] = dt / 6.0 * (k_i_r_1 + 2.0 * k_i_r_2 + 2.0 * k_i_r_3 + k_i_r_4)
    i1_to_r[i] = dt / 6.0 * (k_i1_r_1 + 2.0 * k_i1_r_2 + 2.0 * k_i1_r_3 + k_i1_r_4) 
    i2_to_r[i] = dt / 6.0 * (k_i2_r_1 + 2.0 * k_i2_r_2 + 2.0 * k_i2_r_3 + k_i2_r_4)


k = 0

for j in range(0, days_divided+1):
    if j % int(1.0 / dt) == 0:
        sus_cal[k] = sus[j] 
        vac1_cal[k] = vac1[j]
        vac2_cal[k] = vac2[j]
        exp_cal[k] = exp[j] 
        exp_v1_cal[k] = exp_v1[j]
        exp_v2_cal[k] = exp_v2[j]
        inf_cal[k] = inf[j] 
        inf_v1_cal[k] = inf_v1[j]
        inf_v2_cal[k] = inf_v2[j]
        rem_cal[k] = rem[j] 
        k += 1

k = 0

for j in range(0, days_divided):
    if j % int(1.0 / dt) == 0:
        k += 1
    inf_new_cal[k-1] += e_to_i[j]
    inf_v1_new_cal[k-1] += e1_to_i1[j]
    inf_v2_new_cal[k-1] += e2_to_i2[j]

beta_data = pd.DataFrame(beta_daily * n_total, columns= ['beta'])
date_data = pd.DataFrame(Spreading_data.iloc[0:len(Spreading_data) - 2]['date'])

R2_I = r2_score(Spreading_data.iloc[0:days + 1]['I_now'], inf_cal)
R2_I_new = r2_score(Spreading_data.iloc[0:days + 1]['I_new'], inf_new_cal)
R2_R = r2_score(Spreading_data.iloc[0:days + 1]['Rem'], rem_cal)
print(R2_I)
print(R2_I_new)
print(R2_R)

outputdata_beta = pd.concat([date_data, beta_data], axis=1)
outputdata_beta.to_csv("beta_fm_exIandR.csv", index=False)


sus_output = pd.DataFrame(sus_cal, columns = ['S'])
vac1_output = pd.DataFrame(vac1_cal, columns = ['V1'])
vac2_output = pd.DataFrame(vac2_cal, columns = ['V2'])
exp_output = pd.DataFrame(exp_cal, columns = ['E'])
exp_v1_output = pd.DataFrame(exp_v1_cal, columns = ['E1'])
exp_v2_output = pd.DataFrame(exp_v2_cal, columns = ['E1'])
inf_output = pd.DataFrame(inf_cal, columns = ['I'])
inf_v1_output = pd.DataFrame(inf_v1_cal, columns = ['I1'])
inf_v2_output = pd.DataFrame(inf_v2_cal, columns = ['I2'])
rem_output = pd.DataFrame(rem_cal, columns = ['R'])

outputdata_SEIR = pd.concat([sus_output, vac1_output,vac2_output, exp_output, exp_v1_output,  exp_v2_output, inf_output, inf_v1_output, inf_v2_output, rem_output], axis=1)
outputdata_SEIR.to_csv("SEIR_change.csv", index=False)

inf_new_output = pd.DataFrame(inf_new_cal, columns = ['I_new'])
inf_v1_new_output = pd.DataFrame(inf_v1_new_cal, columns = ['I1_new'])
inf_v2_new_output = pd.DataFrame(inf_v2_new_cal, columns = ['I2_new'])

outputdata_I_new = pd.concat([inf_new_output, inf_v1_new_output, inf_v2_new_output], axis=1)
outputdata_I_new.to_csv("I_new.csv", index=False)




