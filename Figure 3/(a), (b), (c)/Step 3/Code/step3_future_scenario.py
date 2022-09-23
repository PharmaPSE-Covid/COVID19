
'''Calling library'''
import numpy as np
import pandas as pd
import datetime

print('start_____' + str(datetime.datetime.now().hour) + ':' + str(datetime.datetime.now().minute) + ':' + str(datetime.datetime.now().second))

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
beta_fm_exIandR = pd.read_csv('beta_fm_exIandR.csv')

# Time and the number of days
days_pred = 180
days_divided_pred = days_pred * 10
days = len(Spreading_data) - 2
days_divided = days * 10
dt = 1.0 / 10.0 # [day]


'''2020/1/22～2022/1/1'''

# Initial conditions of compartment
n_total = 1.263 * 1.0e+8
    # V,E,I,R
exp_ini =  0.1786
exp_v1_ini = 0.0
exp_v2_ini = 0.0
inf_ini = Spreading_data.iloc[0]['I_now']
inf_v1_ini = 0.0
inf_v2_ini = 0.0
rem_ini = Spreading_data.iloc[0]['Rem']
vac1_ini = 0.0
vac2_ini = 0.0
    # S
sus_ini = n_total - vac1_ini - vac2_ini - exp_ini - exp_v1_ini- exp_v2_ini - inf_ini - inf_v1_ini - inf_v2_ini - rem_ini


# Nature of disease before 2022/1/1
    #alpha
alpha = 0.40
alpha_v1 = 0.40
alpha_v2 = 0.40
    #gamma
gamma = 0.0691712976652005
gamma_v1 = 0.0691712976652005
gamma_v2 = 0.0691712976652005
    #beta
beta = np.zeros(days_divided)
day_counter = -1
for i in range(days_divided):
    if i % int(1.0 / dt) == 0:
        day_counter += 1
    beta[i] = beta_fm_exIandR.iloc[day_counter]['beta'] / n_total
    #severity rate
sev_rate = 0.01 
sev_rate_v1 = 0.001 
sev_rate_v2 = 0.001 


# Conditions about vaccine and antibody before 2022/1/1
    #efficacy
effi_1st = 0.53
effi_2nd = 0.87
    #antibody loss
anti_loss_r_period = 150 #days
anti_loss_v1_period = 150 #days
anti_loss_v2_period = 150 #days
anti_loss_v3_period = 150 #days
    #speed of vaccination [shot/days]
v_speed_1st = np.zeros(days_divided)
v_speed_2nd = np.zeros(days_divided)
day_counter = -1
for i in range(days_divided):
    if i % int(1.0 / dt) == 0:
        day_counter += 1
    v_speed_1st[i] = Spreading_data.iloc[day_counter]['v_speed_1st']
    v_speed_2nd[i] = Spreading_data.iloc[day_counter]['v_speed_2nd']


# Matrix for recording
    #about compartment
sus = np.zeros(days_divided + 1)
vac1 = np.zeros(days_divided + 1)
vac2 = np.zeros(days_divided + 1)
exp = np.zeros(days_divided + 1)
exp_v1 = np.zeros(days_divided + 1)
exp_v2 = np.zeros(days_divided + 1)
inf = np.zeros(days_divided + 1)
inf_v1 = np.zeros(days_divided + 1)
inf_v2 = np.zeros(days_divided + 1)
rem = np.zeros(days_divided + 1)
    #daily compartment
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
inf_total_cal = np.zeros(days+1)
sev_cal = np.zeros(days+1)
    #about arrow
s_to_e = np.zeros(days_divided)
v1_to_e1 = np.zeros(days_divided)
v2_to_e2 = np.zeros(days_divided)
e_to_i = np.zeros(days_divided)
e1_to_i1 = np.zeros(days_divided)
e2_to_i2 = np.zeros(days_divided)
i_to_r = np.zeros(days_divided)
i1_to_r = np.zeros(days_divided)
i2_to_r = np.zeros(days_divided)
anti_loss_r = np.zeros(days_divided)
anti_loss_v1 = np.zeros(days_divided)
anti_loss_v2 = np.zeros(days_divided)
    #daily arrow
inf_new_cal = np.zeros(days)
inf_v1_new_cal = np.zeros(days)
inf_v2_new_cal = np.zeros(days)
inf_total_new_cal = np.zeros(days)
sev_new_cal = np.zeros(days)
anti_loss_r_daily = np.zeros(days)
anti_loss_v1_daily = np.zeros(days)
anti_loss_v2_daily = np.zeros(days)
    #connect
anti_loss_r_daily_connect = np.zeros(days_pred)
anti_loss_v1_daily_connect = np.zeros(days_pred)
anti_loss_v2_daily_connect = np.zeros(days_pred)


# Prepare for calculations
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


# Numerical solution: Runge-Kutta method
day_counter = 0
for i in range(days_divided):
    
    if i % int(1.0 / dt) == 0:

        if i > 0 and day_counter + anti_loss_r_period < days:
            i_r_sum = 0
            for j in range(10):
                i_r_sum += i_to_r[i+j-10] + i1_to_r[i+j-10] + i2_to_r[i+j-10]
            if i_r_sum > 0:
                anti_loss_r_daily[day_counter + anti_loss_r_period] = i_r_sum

        if i > 0 and day_counter + anti_loss_r_period >= days and day_counter + anti_loss_r_period - days < days_pred:
            i_r_sum = 0
            for j in range(10):
                i_r_sum += i_to_r[i+j-10] + i1_to_r[i+j-10] + i2_to_r[i+j-10]
            if i_r_sum > 0:
                anti_loss_r_daily_connect[day_counter + anti_loss_r_period - days] = i_r_sum

        if i > 0 and day_counter + anti_loss_v1_period < days:
            v1_e1_sum = 0
            for j in range(10):
                v1_e1_sum += v1_to_e1[i+j-10]
            if Spreading_data.iloc[day_counter-1]['v_speed_1st'] * 0.003 - v1_e1_sum > 0:
                anti_loss_v1_daily[day_counter + anti_loss_v1_period] = Spreading_data.iloc[day_counter-1]['v_speed_1st'] * 0.003 - v1_e1_sum #0.003では無い方が良いかも？

        if i > 0 and day_counter + anti_loss_v1_period >= days and day_counter + anti_loss_v1_period - days < days_pred :
            v1_e1_sum = 0
            for j in range(10):
                v1_e1_sum += v1_to_e1[i+j-10]
            if Spreading_data.iloc[day_counter-1]['v_speed_1st'] * 0.003 - v1_e1_sum > 0:
                anti_loss_v1_daily_connect[day_counter + anti_loss_v1_period - days] = Spreading_data.iloc[day_counter-1]['v_speed_1st'] * 0.003 - v1_e1_sum #0.003では無い方が良いかも？

        if i > 0 and day_counter + anti_loss_v2_period < days:
            v2_e2_sum = 0
            for j in range(10):
                v2_e2_sum += v2_to_e2[i+j-10]
            if Spreading_data.iloc[day_counter-1]['v_speed_2nd'] - v2_e2_sum > 0:
                anti_loss_v2_daily[day_counter + anti_loss_v2_period] = Spreading_data.iloc[day_counter-1]['v_speed_2nd'] - v2_e2_sum

        if i > 0 and day_counter + anti_loss_v2_period >= days and day_counter + anti_loss_v2_period - days < days_pred :
            v2_e2_sum = 0
            for j in range(10):
                v2_e2_sum += v2_to_e2[i+j-10]
            if Spreading_data.iloc[day_counter-1]['v_speed_2nd'] - v2_e2_sum > 0:
                anti_loss_v2_daily_connect[day_counter + anti_loss_v2_period - days] = Spreading_data.iloc[day_counter-1]['v_speed_2nd'] - v2_e2_sum
        
        day_counter += 1  

    anti_loss_r[i] = anti_loss_r_daily[day_counter-1]
    anti_loss_v1[i] = anti_loss_v1_daily[day_counter-1]
    anti_loss_v2[i] = anti_loss_v2_daily[day_counter-1]

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
            
    s_to_e[i] = dt / 6.0 * (k_s_e_1 + 2.0 * k_s_e_2 + 2.0 * k_s_e_3 + k_s_e_4)
    v1_to_e1[i] = dt / 6.0 * (k_v1_e1_1 + 2.0 * k_v1_e1_2 + 2.0 * k_v1_e1_3 + k_v1_e1_4)
    v2_to_e2[i] = dt / 6.0 * (k_v2_e2_1 + 2.0 * k_v2_e2_2 + 2.0 * k_v2_e2_3 + k_v2_e2_4)
    e_to_i[i] = dt / 6.0 * (k_e_i_1 + 2.0 * k_e_i_2 + 2.0 * k_e_i_3 + k_e_i_4)
    e1_to_i1[i] = dt / 6.0 * (k_e1_i1_1 + 2.0 * k_e1_i1_2 + 2.0 * k_e1_i1_3 + k_e1_i1_4)
    e2_to_i2[i] = dt / 6.0 * (k_e2_i2_1 + 2.0 * k_e2_i2_2 + 2.0 * k_e2_i2_3 + k_e2_i2_4)
    i_to_r[i] = dt / 6.0 * (k_i_r_1 + 2.0 * k_i_r_2 + 2.0 * k_i_r_3 + k_i_r_4)
    i1_to_r[i] = dt / 6.0 * (k_i1_r_1 + 2.0 * k_i1_r_2 + 2.0 * k_i1_r_3 + k_i1_r_4) 
    i2_to_r[i] = dt / 6.0 * (k_i2_r_1 + 2.0 * k_i2_r_2 + 2.0 * k_i2_r_3 + k_i2_r_4)


# daily recording
k = 0
for j in range(0, days_divided + 1):
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
        inf_total_cal[k] = inf[j] + inf_v1[j] + inf_v2[j]
        sev_cal[k] = inf[j] * sev_rate + inf_v1[j] * sev_rate_v1 + inf_v2[j] * sev_rate_v2
        k += 1

k = 0
for j in range(0, days_divided):
    if j % int(1.0 / dt) == 0:
        k += 1
    inf_new_cal[k-1] += e_to_i[j]
    inf_v1_new_cal[k-1] += e1_to_i1[j]
    inf_v2_new_cal[k-1] += e2_to_i2[j]
    inf_total_new_cal[k-1] += e_to_i[j] + e1_to_i1[j] + e2_to_i2[j]
    sev_new_cal[k-1] += e_to_i[j] * sev_rate + e1_to_i1[j] * sev_rate_v1 + e2_to_i2[j] * sev_rate_v2


'''2022/1/1～'''

print('_______________________________')
print('step1_____' + str(datetime.datetime.now().hour) + ':' + str(datetime.datetime.now().minute) + ':' + str(datetime.datetime.now().second))


'''Functions'''
def f_s_e_pred(sus_old, inf_old, inf_v1_old, inf_v2_old, inf_v3_old, beta):
    x = beta * sus_old * (inf_old + inf_v1_old + inf_v2_old + inf_v3_old)
    return x

def f_v1_e1_pred(vac1_old, inf_old, inf_v1_old, inf_v2_old, inf_v3_old, beta, effi_1st):
    x = beta * vac1_old * (inf_old + inf_v1_old + inf_v2_old + inf_v3_old) * (1 - effi_1st)
    return x

def f_v2_e2_pred(vac2_old, inf_old, inf_v1_old, inf_v2_old, inf_v3_old, beta, effi_2nd):
    x = beta * vac2_old * (inf_old + inf_v1_old + inf_v2_old + inf_v3_old) * (1 - effi_2nd)
    return x

def f_v3_e3_pred(vac3_old, inf_old, inf_v1_old, inf_v2_old, inf_v3_old, beta, effi_3rd):
    x = beta * vac3_old * (inf_old + inf_v1_old + inf_v2_old + inf_v3_old) * (1 - effi_3rd)
    return x

def f_e3_i3_pred(exp_v3_old, alpha_v3):
    x = alpha_v3 * exp_v3_old
    return x

def f_i3_r_pred(inf_v3_old, gamma_v3):
    x = gamma_v3 * inf_v3_old
    return x


'''Main'''

# Upload csv files
scenario_data = pd.read_csv('scenario_data.csv', index_col=0)

number_of_scenario = len(scenario_data.columns)
if number_of_scenario > 10:
    timecounter = number_of_scenario // 10
if number_of_scenario <= 10:
    timecounter = 1
print('number of scenario = ' + str(number_of_scenario))
print('_______________________________')


#Record Scenario
clm_I = scenario_data.columns[0:]
ind_I = np.arange(pd.to_datetime(beta_fm_exIandR.iloc[-1]['date']), pd.to_datetime(beta_fm_exIandR.iloc[-1]['date']) + np.timedelta64(days_pred,'D'), np.timedelta64(1,'D'), dtype='datetime64')
I_new_pred = pd.DataFrame(index = ind_I, columns = clm_I)
I_now_pred = pd.DataFrame(index = ind_I, columns = clm_I)
sev_now_pred = pd.DataFrame(index = ind_I, columns = clm_I)


for scenario_number in range(number_of_scenario):
    if scenario_number % timecounter == 0:
        print(str(scenario_number) + '/' + str(number_of_scenario) + ' --- ' + str(datetime.datetime.now().hour) + ':' + str(datetime.datetime.now().minute) + ':' + str(datetime.datetime.now().second))

    # Nature of disease after 2022/1/1
        #alpha
    alpha = 1 / scenario_data.iloc[0, scenario_number]
    alpha_v1 = 1 / scenario_data.iloc[1, scenario_number]
    alpha_v2 = 1 / scenario_data.iloc[2, scenario_number]
    alpha_v3 = 1 / scenario_data.iloc[3, scenario_number]
        #gamma
    gamma = 1 / scenario_data.iloc[5, scenario_number]
    gamma_v1 = 1 / scenario_data.iloc[6, scenario_number]
    gamma_v2 = 1 / scenario_data.iloc[7, scenario_number]
    gamma_v3 = 1 / scenario_data.iloc[8, scenario_number]
        #beta
    beta = scenario_data.iloc[4, scenario_number] / n_total
    beta_emergency = 0.64
    emergency_start = 1
    emergency_end = 200
        #severity rate
    sev_rate = scenario_data.iloc[12, scenario_number]
    sev_rate_v1 = scenario_data.iloc[13, scenario_number]
    sev_rate_v2 = scenario_data.iloc[14, scenario_number]
    sev_rate_v3 = scenario_data.iloc[15, scenario_number]

    # Conditions about vaccine and antibody after 2022/1/1
        #efficacy
    effi_1st = scenario_data.iloc[9, scenario_number]
    effi_2nd = scenario_data.iloc[10, scenario_number]
    effi_3rd = scenario_data.iloc[11, scenario_number]

        #speed of vaccination [shot/days]
    v_speed_1st = 0
    v_speed_2nd = 0
    v_speed_3rd_t = scenario_data.iloc[16, scenario_number]
    v_speed_3rd = 0
    v3_delay = scenario_data.iloc[17, scenario_number]
    ve_period = scenario_data.iloc[18, scenario_number]
    v3_end =  v3_delay + 95

    # Matrix for recording
        #about compartment
    sus_pred = np.zeros(days_divided_pred + 1)
    vac1_pred = np.zeros(days_divided_pred + 1)
    vac2_pred = np.zeros(days_divided_pred + 1)
    vac3_pred = np.zeros(days_divided_pred + 1)
    exp_pred = np.zeros(days_divided_pred + 1)
    exp_v1_pred = np.zeros(days_divided_pred + 1)
    exp_v2_pred = np.zeros(days_divided_pred + 1)
    exp_v3_pred = np.zeros(days_divided_pred + 1)
    inf_pred = np.zeros(days_divided_pred + 1)
    inf_v1_pred = np.zeros(days_divided_pred + 1)
    inf_v2_pred = np.zeros(days_divided_pred + 1)
    inf_v3_pred = np.zeros(days_divided_pred + 1)
    rem_pred = np.zeros(days_divided_pred + 1)
        #daily compartment
    sus_cal_pred = np.zeros(days_pred+1)
    vac1_cal_pred = np.zeros(days_pred+1)
    vac2_cal_pred = np.zeros(days_pred+1)
    vac3_cal_pred = np.zeros(days_pred+1)
    exp_cal_pred = np.zeros(days_pred+1)
    exp_v1_cal_pred = np.zeros(days_pred+1)
    exp_v2_cal_pred = np.zeros(days_pred+1)
    exp_v3_cal_pred = np.zeros(days_pred+1)
    inf_cal_pred = np.zeros(days_pred+1)
    inf_v1_cal_pred = np.zeros(days_pred+1)
    inf_v2_cal_pred = np.zeros(days_pred+1)
    inf_v3_cal_pred = np.zeros(days_pred+1)
    rem_cal_pred = np.zeros(days_pred+1)
    inf_total_cal_pred = np.zeros(days_pred+1)
    sev_cal_pred = np.zeros(days_pred+1)
        #about arrow
    s_to_e_pred = np.zeros(days_divided_pred)
    v1_to_e1_pred = np.zeros(days_divided_pred)
    v2_to_e2_pred = np.zeros(days_divided_pred)
    v3_to_e3_pred = np.zeros(days_divided_pred)
    e_to_i_pred = np.zeros(days_divided_pred)
    e1_to_i1_pred = np.zeros(days_divided_pred)
    e2_to_i2_pred = np.zeros(days_divided_pred)
    e3_to_i3_pred = np.zeros(days_divided_pred)
    i_to_r_pred = np.zeros(days_divided_pred)
    i1_to_r_pred = np.zeros(days_divided_pred)
    i2_to_r_pred = np.zeros(days_divided_pred)
    i3_to_r_pred = np.zeros(days_divided_pred)
    anti_loss_r_pred = np.zeros(days_divided_pred)
    anti_loss_v1_pred = np.zeros(days_divided_pred)
    anti_loss_v2_pred = np.zeros(days_divided_pred)
    anti_loss_v3_pred = np.zeros(days_divided_pred)
        #daily arrow
    inf_new_cal_pred = np.zeros(days_pred)
    inf_v1_new_cal_pred = np.zeros(days_pred)
    inf_v2_new_cal_pred = np.zeros(days_pred)
    inf_v3_new_cal_pred = np.zeros(days_pred)
    inf_total_new_cal_pred = np.zeros(days_pred)
    sev_new_cal_pred = np.zeros(days_pred)
    anti_loss_r_daily_pred = np.zeros(days_pred)
    anti_loss_v1_daily_pred = np.zeros(days_pred)
    anti_loss_v2_daily_pred = np.zeros(days_pred)
    anti_loss_v3_daily_pred = np.zeros(days_pred)

    
    # Prepare for calculations
    sus_old = sus[-1]
    vac1_old = vac1[-1]
    vac2_old = vac2[-1]
    vac3_old = 0
    exp_old = exp[-1]
    exp_v1_old = exp_v1[-1]
    exp_v2_old = exp_v2[-1]
    exp_v3_old = 0
    inf_old = inf[-1]
    inf_v1_old = inf_v1[-1]
    inf_v2_old = inf_v2[-1]
    inf_v3_old = 0
    rem_old = rem[-1]
    
    sus_pred[0] = sus[-1]
    vac1_pred[0] = vac1[-1]
    vac2_pred[0] = vac2[-1]
    vac3_pred[0] = 0
    exp_pred[0] = exp[-1]
    exp_v1_pred[0] = exp_v1[-1]
    exp_v2_pred[0] = exp_v2[-1]
    exp_v3_pred[0] = 0
    inf_pred[0] = inf[-1]
    inf_v1_pred[0] = inf_v1[-1]
    inf_v2_pred[0] = inf_v2[-1]
    inf_v2_pred[0] = 0
    rem_pred[0] = rem[-1]


    # Numerical solution: Runge-Kutta method
    day_counter = 0
    for i in range(days_divided_pred):
        
        if i % int(1.0 / dt) == 0:
            
            if day_counter == emergency_start - 1:
                beta *= beta_emergency
            if day_counter == emergency_end - 1:
                beta /= beta_emergency 
            
            if day_counter == v3_delay:
                v_speed_3rd = v_speed_3rd_t
            if day_counter == v3_end:
                v_speed_3rd = 0
            
            if i > 0 and day_counter + anti_loss_r_period < days_pred:
                i_r_sum = 0
                for j in range(10):
                    i_r_sum += i_to_r_pred[i+j-10] + i1_to_r_pred[i+j-10] + i2_to_r_pred[i+j-10] + i3_to_r_pred[i+j-10]
                anti_loss_r_daily_pred[day_counter + anti_loss_r_period] = i_r_sum
            
            if i > 0 and day_counter + anti_loss_v1_period < days_pred:
                v1_e1_sum = 0
                for j in range(10):
                    v1_e1_sum += v1_to_e1_pred[i+j-10]
                anti_loss_v1_daily_pred[day_counter + anti_loss_v1_period] = 0
            
            if i > 0 and day_counter + anti_loss_v2_period < days_pred:
                v2_e2_sum = 0
                for j in range(10):
                    v2_e2_sum += v2_to_e2_pred[i+j-10]
                anti_loss_v2_daily_pred[day_counter + anti_loss_v2_period] = 0
                
            if i > 0 and day_counter + anti_loss_v3_period < days_pred:
                v3_e3_sum = 0
                for j in range(10):
                    v3_e3_sum += v3_to_e3_pred[i+j-10]
                anti_loss_v3_daily_pred[day_counter + anti_loss_v3_period] = 0
            
            day_counter += 1  
        
        if anti_loss_r_daily_pred[day_counter-1] == 0:
            anti_loss_r_daily_pred[day_counter-1] = anti_loss_r_daily_connect[day_counter-1]
        if anti_loss_v1_daily_pred[day_counter-1] == 0:
            anti_loss_v1_daily_pred[day_counter-1] = anti_loss_v1_daily_connect[day_counter-1]
        if anti_loss_v2_daily_pred[day_counter-1] == 0:
            anti_loss_v2_daily_pred[day_counter-1] = anti_loss_v2_daily_connect[day_counter-1]
        
        anti_loss_r_pred[i] = anti_loss_r_daily_pred[day_counter-1]
        anti_loss_v1_pred[i] = anti_loss_v1_daily_pred[day_counter-1]
        anti_loss_v2_pred[i] = anti_loss_v2_daily_pred[day_counter-1]
        anti_loss_v3_pred[i] = anti_loss_v3_daily_pred[day_counter-1]
    
        k_s_e_1 = f_s_e_pred(sus_old, inf_old, inf_v1_old, inf_v2_old, inf_v3_old, beta)
        k_v1_e1_1 = f_v1_e1_pred(vac1_old, inf_old, inf_v1_old, inf_v2_old, inf_v3_old, beta, effi_1st)
        k_v2_e2_1 = f_v2_e2_pred(vac2_old, inf_old, inf_v1_old, inf_v2_old, inf_v3_old, beta, effi_2nd)
        k_v3_e3_1 = f_v3_e3_pred(vac3_old, inf_old, inf_v1_old, inf_v2_old, inf_v3_old, beta, effi_3rd)
        k_e_i_1 = f_e_i(exp_old, alpha)
        k_e1_i1_1 = f_e1_i1(exp_v1_old, alpha_v1)
        k_e2_i2_1 = f_e2_i2(exp_v2_old, alpha_v2)
        k_e3_i3_1 = f_e3_i3_pred(exp_v3_old, alpha_v3)
        k_i_r_1 = f_i_r(inf_old, gamma)
        k_i1_r_1 = f_i1_r(inf_v1_old, gamma_v1)
        k_i2_r_1 = f_i2_r(inf_v2_old, gamma_v2)
        k_i3_r_1 = f_i3_r_pred(inf_v3_old, gamma_v3)
    
        k_sus_1 = - k_s_e_1 - v_speed_1st - v_speed_3rd + anti_loss_r_pred[i] + anti_loss_v1_pred[i] + anti_loss_v2_pred[i] + anti_loss_v3_pred[i]
        k_vac1_1 = v_speed_1st - v_speed_2nd - k_v1_e1_1 - anti_loss_v1_pred[i]
        k_vac2_1 = v_speed_2nd - k_v2_e2_1 - anti_loss_v2_pred[i]
        k_vac3_1 = v_speed_3rd - k_v3_e3_1 - anti_loss_v3_pred[i]
        k_exp_1 = k_s_e_1 - k_e_i_1
        k_exp_v1_1 = k_v1_e1_1 - k_e1_i1_1
        k_exp_v2_1 = k_v2_e2_1 - k_e2_i2_1
        k_exp_v3_1 = k_v3_e3_1 - k_e3_i3_1
        k_inf_1 = k_e_i_1 - k_i_r_1
        k_inf_v1_1 = k_e1_i1_1 - k_i1_r_1
        k_inf_v2_1 = k_e2_i2_1 - k_i2_r_1
        k_inf_v3_1 = k_e3_i3_1 - k_i3_r_1
        k_rem_1 = k_i_r_1 + k_i1_r_1 + k_i2_r_1 + k_i3_r_1 - anti_loss_r_pred[i]
    
        k_s_e_2 = f_s_e_pred(sus_old + k_sus_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, inf_v1_old + k_inf_v1_1 * 0.5 * dt, inf_v2_old + k_inf_v2_1 * 0.5 * dt, inf_v3_old + k_inf_v3_1 * 0.5 * dt, beta)
        k_v1_e1_2 = f_v1_e1_pred(vac1_old + k_vac1_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, inf_v1_old + k_inf_v1_1 * 0.5 * dt, inf_v2_old + k_inf_v2_1 * 0.5 * dt, inf_v3_old + k_inf_v3_1 * 0.5 * dt, beta, effi_1st)
        k_v2_e2_2 = f_v2_e2_pred(vac2_old + k_vac2_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, inf_v1_old + k_inf_v1_1 * 0.5 * dt, inf_v2_old + k_inf_v2_1 * 0.5 * dt, inf_v3_old + k_inf_v3_1 * 0.5 * dt, beta, effi_2nd)
        k_v3_e3_2 = f_v3_e3_pred(vac3_old + k_vac3_1 * 0.5 * dt, inf_old + k_inf_1 * 0.5 * dt, inf_v1_old + k_inf_v1_1 * 0.5 * dt, inf_v2_old + k_inf_v2_1 * 0.5 * dt, inf_v3_old + k_inf_v3_1 * 0.5 * dt, beta, effi_3rd)
        k_e_i_2 = f_e_i(exp_old + k_exp_1 * 0.5 * dt, alpha)
        k_e1_i1_2 = f_e1_i1(exp_v1_old + k_exp_v1_1 * 0.5 * dt, alpha_v1)
        k_e2_i2_2 = f_e2_i2(exp_v2_old + k_exp_v2_1 * 0.5 * dt, alpha_v2)
        k_e3_i3_2 = f_e3_i3_pred(exp_v3_old + k_exp_v3_1 * 0.5 * dt, alpha_v3)
        k_i_r_2 = f_i_r(inf_old + k_inf_1 * 0.5 * dt, gamma)
        k_i1_r_2 = f_i1_r(inf_v1_old + k_inf_v1_1 * 0.5 * dt, gamma_v1)
        k_i2_r_2 = f_i2_r(inf_v2_old + k_inf_v2_1 * 0.5 * dt, gamma_v2)
        k_i3_r_2 = f_i3_r_pred(inf_v3_old + k_inf_v3_1 * 0.5 * dt, gamma_v3)
        
        k_sus_2 = - k_s_e_2 - v_speed_1st - v_speed_3rd + anti_loss_r_pred[i] + anti_loss_v1_pred[i] + anti_loss_v2_pred[i] + anti_loss_v3_pred[i]
        k_vac1_2 = v_speed_1st - v_speed_2nd - k_v1_e1_2 - anti_loss_v1_pred[i]
        k_vac2_2 = v_speed_2nd - k_v2_e2_2 - anti_loss_v2_pred[i]
        k_vac3_2 = v_speed_3rd - k_v3_e3_2 - anti_loss_v3_pred[i]
        k_exp_2 = k_s_e_2 - k_e_i_2
        k_exp_v1_2 = k_v1_e1_2 - k_e1_i1_2
        k_exp_v2_2 = k_v2_e2_2 - k_e2_i2_2
        k_exp_v3_2 = k_v3_e3_2 - k_e3_i3_2
        k_inf_2 = k_e_i_2 - k_i_r_2
        k_inf_v1_2 = k_e1_i1_2 - k_i1_r_2
        k_inf_v2_2 = k_e2_i2_2 - k_i2_r_2
        k_inf_v3_2 = k_e3_i3_2 - k_i3_r_2
        k_rem_2 = k_i_r_2 + k_i1_r_2 + k_i2_r_2 + k_i3_r_2 - anti_loss_r_pred[i]
    
        k_s_e_3 = f_s_e_pred(sus_old + k_sus_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, inf_v1_old + k_inf_v1_2 * 0.5 * dt, inf_v2_old + k_inf_v2_2 * 0.5 * dt, inf_v3_old + k_inf_v3_2 * 0.5 * dt, beta)
        k_v1_e1_3 = f_v1_e1_pred(vac1_old + k_vac1_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, inf_v1_old + k_inf_v1_2 * 0.5 * dt, inf_v2_old + k_inf_v2_2 * 0.5 * dt, inf_v3_old + k_inf_v3_2 * 0.5 * dt, beta, effi_1st)
        k_v2_e2_3 = f_v2_e2_pred(vac2_old + k_vac2_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, inf_v1_old + k_inf_v1_2 * 0.5 * dt, inf_v2_old + k_inf_v2_2 * 0.5 * dt, inf_v3_old + k_inf_v3_2 * 0.5 * dt, beta, effi_2nd)
        k_v3_e3_3 = f_v3_e3_pred(vac3_old + k_vac3_2 * 0.5 * dt, inf_old + k_inf_2 * 0.5 * dt, inf_v1_old + k_inf_v1_2 * 0.5 * dt, inf_v2_old + k_inf_v2_2 * 0.5 * dt, inf_v3_old + k_inf_v3_2 * 0.5 * dt, beta, effi_3rd)
        k_e_i_3 = f_e_i(exp_old + k_exp_2 * 0.5 * dt, alpha)
        k_e1_i1_3 = f_e1_i1(exp_v1_old + k_exp_v1_2 * 0.5 * dt, alpha_v1)
        k_e2_i2_3 = f_e2_i2(exp_v2_old + k_exp_v2_2 * 0.5 * dt, alpha_v2)
        k_e3_i3_3 = f_e3_i3_pred(exp_v3_old + k_exp_v3_2 * 0.5 * dt, alpha_v3)
        k_i_r_3 = f_i_r(inf_old + k_inf_2 * 0.5 * dt, gamma)
        k_i1_r_3 = f_i1_r(inf_v1_old + k_inf_v1_2 * 0.5 * dt, gamma_v1)
        k_i2_r_3 = f_i2_r(inf_v2_old + k_inf_v2_2 * 0.5 * dt, gamma_v2)
        k_i3_r_3 = f_i3_r_pred(inf_v3_old + k_inf_v3_2 * 0.5 * dt, gamma_v3)
        
        k_sus_3 = - k_s_e_3 - v_speed_1st - v_speed_3rd + anti_loss_r_pred[i] + anti_loss_v1_pred[i] + anti_loss_v2_pred[i] + anti_loss_v3_pred[i]
        k_vac1_3 = v_speed_1st - v_speed_2nd - k_v1_e1_3 - anti_loss_v1_pred[i]
        k_vac2_3 = v_speed_2nd - k_v2_e2_3 - anti_loss_v2_pred[i]
        k_vac3_3 = v_speed_3rd - k_v3_e3_3 - anti_loss_v3_pred[i]
        k_exp_3 = k_s_e_3 - k_e_i_3
        k_exp_v1_3 = k_v1_e1_3 - k_e1_i1_3
        k_exp_v2_3 = k_v2_e2_3 - k_e2_i2_3
        k_exp_v3_3 = k_v3_e3_3 - k_e3_i3_3
        k_inf_3 = k_e_i_3 - k_i_r_3
        k_inf_v1_3 = k_e1_i1_3 - k_i1_r_3
        k_inf_v2_3 = k_e2_i2_3 - k_i2_r_3
        k_inf_v3_3 = k_e3_i3_3 - k_i3_r_3
        k_rem_3 = k_i_r_3 + k_i1_r_3 + k_i2_r_3 + k_i3_r_3 - anti_loss_r_pred[i]
    
        k_s_e_4 = f_s_e_pred(sus_old + k_sus_3 * dt, inf_old + k_inf_3 * dt, inf_v1_old + k_inf_v1_3 * dt, inf_v2_old + k_inf_v2_3 * dt, inf_v3_old + k_inf_v3_3 * dt, beta)
        k_v1_e1_4 = f_v1_e1_pred(vac1_old + k_vac1_3 * dt, inf_old + k_inf_3 * dt, inf_v1_old + k_inf_v1_3 * dt, inf_v2_old + k_inf_v2_3 * dt, inf_v3_old + k_inf_v3_3 * dt, beta, effi_1st)
        k_v2_e2_4 = f_v2_e2_pred(vac2_old + k_vac2_3 * dt, inf_old + k_inf_3 * dt, inf_v1_old + k_inf_v1_3 * dt, inf_v2_old + k_inf_v2_3 * dt, inf_v3_old + k_inf_v3_3 * dt, beta, effi_2nd)
        k_v3_e3_4 = f_v3_e3_pred(vac3_old + k_vac3_3 * dt, inf_old + k_inf_3 * dt, inf_v1_old + k_inf_v1_3 * dt, inf_v2_old + k_inf_v2_3 * dt, inf_v3_old + k_inf_v3_3 * dt, beta, effi_3rd)
        k_e_i_4 = f_e_i(exp_old + k_exp_3 * dt, alpha)
        k_e1_i1_4 = f_e1_i1(exp_v1_old + k_exp_v1_3 * dt, alpha_v1)
        k_e2_i2_4 = f_e2_i2(exp_v2_old + k_exp_v2_3 * dt, alpha_v2)
        k_e3_i3_4 = f_e3_i3_pred(exp_v3_old + k_exp_v3_3 * dt, alpha_v3)
        k_i_r_4 = f_i_r(inf_old + k_inf_3 * dt, gamma)
        k_i1_r_4 = f_i1_r(inf_v1_old + k_inf_v1_3 * dt, gamma_v1)
        k_i2_r_4 = f_i2_r(inf_v2_old + k_inf_v2_3 * dt, gamma_v2)
        k_i3_r_4 = f_i3_r_pred(inf_v3_old + k_inf_v3_3 * dt, gamma_v3)
        
        k_sus_4 = - k_s_e_4 - v_speed_1st - v_speed_3rd + anti_loss_r_pred[i] + anti_loss_v1_pred[i] + anti_loss_v2_pred[i] + anti_loss_v3_pred[i]
        k_vac1_4 = v_speed_1st - v_speed_2nd - k_v1_e1_4 - anti_loss_v1_pred[i]
        k_vac2_4 = v_speed_2nd - k_v2_e2_4 - anti_loss_v2_pred[i]
        k_vac3_4 = v_speed_3rd - k_v3_e3_4 - anti_loss_v3_pred[i]
        k_exp_4 = k_s_e_4 - k_e_i_4
        k_exp_v1_4 = k_v1_e1_4 - k_e1_i1_4
        k_exp_v2_4 = k_v2_e2_4 - k_e2_i2_4
        k_exp_v3_4 = k_v3_e3_4 - k_e3_i3_4
        k_inf_4 = k_e_i_4 - k_i_r_4
        k_inf_v1_4 = k_e1_i1_4 - k_i1_r_4
        k_inf_v2_4 = k_e2_i2_4 - k_i2_r_4
        k_inf_v3_4 = k_e3_i3_4 - k_i3_r_4
        k_rem_4 = k_i_r_4 + k_i1_r_4 + k_i2_r_4 + k_i3_r_4 - anti_loss_r_pred[i]
    
        sus_new = sus_old + dt / 6.0 * (k_sus_1 + 2.0 * k_sus_2 + 2.0 * k_sus_3 + k_sus_4)
        vac1_new = vac1_old + dt / 6.0 * (k_vac1_1 + 2.0 * k_vac1_2 + 2.0 * k_vac1_3 + k_vac1_4)
        vac2_new = vac2_old + dt / 6.0 * (k_vac2_1 + 2.0 * k_vac2_2 + 2.0 * k_vac2_3 + k_vac2_4)
        vac3_new = vac3_old + dt / 6.0 * (k_vac3_1 + 2.0 * k_vac3_2 + 2.0 * k_vac3_3 + k_vac3_4)
        exp_new = exp_old + dt / 6.0 * (k_exp_1 + 2.0 * k_exp_2 + 2.0 * k_exp_3 + k_exp_4)
        exp_v1_new = exp_v1_old + dt / 6.0 * (k_exp_v1_1 + 2.0 * k_exp_v1_2 + 2.0 * k_exp_v1_3 + k_exp_v1_4)
        exp_v2_new = exp_v2_old + dt / 6.0 * (k_exp_v2_1 + 2.0 * k_exp_v2_2 + 2.0 * k_exp_v2_3 + k_exp_v2_4)
        exp_v3_new = exp_v3_old + dt / 6.0 * (k_exp_v3_1 + 2.0 * k_exp_v3_2 + 2.0 * k_exp_v3_3 + k_exp_v3_4)
        inf_new = inf_old + dt / 6.0 * (k_inf_1 + 2.0 * k_inf_2 + 2.0 * k_inf_3 + k_inf_4)
        inf_v1_new = inf_v1_old + dt / 6.0 * (k_inf_v1_1 + 2.0 * k_inf_v1_2 + 2.0 * k_inf_v1_3 + k_inf_v1_4)
        inf_v2_new = inf_v2_old + dt / 6.0 * (k_inf_v2_1 + 2.0 * k_inf_v2_2 + 2.0 * k_inf_v2_3 + k_inf_v2_4)
        inf_v3_new = inf_v3_old + dt / 6.0 * (k_inf_v3_1 + 2.0 * k_inf_v3_2 + 2.0 * k_inf_v3_3 + k_inf_v3_4)
        rem_new = rem_old + dt / 6.0 * (k_rem_1 + 2.0 * k_rem_2 + 2.0 * k_rem_3 + k_rem_4)
        
        if sus_new < 0:
            sus_new = 0.0
        if vac1_new < 0:
            vac1_new = 0.0
        if vac2_new < 0:
            vac2_new = 0.0
        if vac3_new < 0:
            vac3_new = 0.0 
        if exp_new < 0:
            exp_new = 0.0
        if exp_v1_new < 0:
            exp_v1_new = 0.0
        if exp_v2_new < 0:
            exp_v2_new = 0.0
        if exp_v3_new < 0:
            exp_v3_new = 0.0
        if inf_new < 0:
            inf_new = 0.0
        if inf_v1_new < 0:
            inf_v1_new = 0.0
        if inf_v2_new < 0:
            inf_v2_new = 0.0
        if inf_v3_new < 0:
            inf_v3_new = 0.0
        if rem_new < 0:
            rem_new = 0.0
        
        sus_old = sus_new
        vac1_old = vac1_new
        vac2_old = vac2_new
        vac3_old = vac3_new
        exp_old = exp_new
        exp_v1_old = exp_v1_new
        exp_v2_old = exp_v2_new
        exp_v3_old = exp_v3_new
        inf_old = inf_new
        inf_v1_old = inf_v1_new
        inf_v2_old = inf_v2_new
        inf_v3_old = inf_v3_new
        rem_old = rem_new
    
        sus_pred[i+1] = sus_old
        vac1_pred[i+1] = vac1_old
        vac2_pred[i+1] = vac2_old
        vac3_pred[i+1] = vac3_old
        exp_pred[i+1] = exp_old
        exp_v1_pred[i+1] = exp_v1_old
        exp_v2_pred[i+1] = exp_v2_old
        exp_v3_pred[i+1] = exp_v3_old
        inf_pred[i+1] = inf_old
        inf_v1_pred[i+1] = inf_v1_old
        inf_v2_pred[i+1] = inf_v2_old
        inf_v3_pred[i+1] = inf_v3_old
        rem_pred[i+1] = rem_old 
                
        s_to_e_pred[i] = dt / 6.0 * (k_s_e_1 + 2.0 * k_s_e_2 + 2.0 * k_s_e_3 + k_s_e_4)
        v1_to_e1_pred[i] = dt / 6.0 * (k_v1_e1_1 + 2.0 * k_v1_e1_2 + 2.0 * k_v1_e1_3 + k_v1_e1_4)
        v2_to_e2_pred[i] = dt / 6.0 * (k_v2_e2_1 + 2.0 * k_v2_e2_2 + 2.0 * k_v2_e2_3 + k_v2_e2_4)
        v3_to_e3_pred[i] = dt / 6.0 * (k_v3_e3_1 + 2.0 * k_v3_e3_2 + 2.0 * k_v3_e3_3 + k_v3_e3_4)
        e_to_i_pred[i] = dt / 6.0 * (k_e_i_1 + 2.0 * k_e_i_2 + 2.0 * k_e_i_3 + k_e_i_4)
        e1_to_i1_pred[i] = dt / 6.0 * (k_e1_i1_1 + 2.0 * k_e1_i1_2 + 2.0 * k_e1_i1_3 + k_e1_i1_4)
        e2_to_i2_pred[i] = dt / 6.0 * (k_e2_i2_1 + 2.0 * k_e2_i2_2 + 2.0 * k_e2_i2_3 + k_e2_i2_4)
        e3_to_i3_pred[i] = dt / 6.0 * (k_e3_i3_1 + 2.0 * k_e3_i3_2 + 2.0 * k_e3_i3_3 + k_e3_i3_4)
        i_to_r_pred[i] = dt / 6.0 * (k_i_r_1 + 2.0 * k_i_r_2 + 2.0 * k_i_r_3 + k_i_r_4)
        i1_to_r_pred[i] = dt / 6.0 * (k_i1_r_1 + 2.0 * k_i1_r_2 + 2.0 * k_i1_r_3 + k_i1_r_4) 
        i2_to_r_pred[i] = dt / 6.0 * (k_i2_r_1 + 2.0 * k_i2_r_2 + 2.0 * k_i2_r_3 + k_i2_r_4)
        i3_to_r_pred[i] = dt / 6.0 * (k_i3_r_1 + 2.0 * k_i3_r_2 + 2.0 * k_i3_r_3 + k_i3_r_4)


    # daily recording
    k = 0
    for j in range(0, days_divided_pred + 1):
        if j % int(1.0 / dt) == 0:
            sus_cal_pred[k] = sus_pred[j] 
            vac1_cal_pred[k] = vac1_pred[j]
            vac2_cal_pred[k] = vac2_pred[j]
            vac3_cal_pred[k] = vac3_pred[j]
            exp_cal_pred[k] = exp_pred[j] 
            exp_v1_cal_pred[k] = exp_v1_pred[j]
            exp_v2_cal_pred[k] = exp_v2_pred[j]
            exp_v3_cal_pred[k] = exp_v3_pred[j]
            inf_cal_pred[k] = inf_pred[j] 
            inf_v1_cal_pred[k] = inf_v1_pred[j]
            inf_v2_cal_pred[k] = inf_v2_pred[j]
            inf_v3_cal_pred[k] = inf_v3_pred[j]
            rem_cal_pred[k] = rem_pred[j]
            inf_total_cal_pred[k] = inf_pred[j] + inf_v1_pred[j] + inf_v2_pred[j] + inf_v3_pred[j]
            sev_cal_pred[k] = inf_pred[j] * sev_rate + inf_v1_pred[j] * sev_rate_v1 + inf_v2_pred[j] * sev_rate_v2 + inf_v3_pred[j] * sev_rate_v3
            k += 1
    
    k = 0
    for j in range(0, days_divided_pred):
        if j % int(1.0 / dt) == 0:
            k += 1
        inf_new_cal_pred[k-1] += e_to_i_pred[j]
        inf_v1_new_cal_pred[k-1] += e1_to_i1_pred[j]
        inf_v2_new_cal_pred[k-1] += e2_to_i2_pred[j]
        inf_v3_new_cal_pred[k-1] += e3_to_i3_pred[j]
        inf_total_new_cal_pred[k-1] += e_to_i_pred[j] + e1_to_i1_pred[j] + e2_to_i2_pred[j] + e3_to_i3_pred[j]
        sev_new_cal_pred[k-1] += e_to_i_pred[j] * sev_rate + e1_to_i1_pred[j] * sev_rate_v1 + e2_to_i2_pred[j] * sev_rate_v2 + e3_to_i3_pred[j] * sev_rate_v3

    for i in range(days_pred):
        I_new_pred.iloc[i, scenario_number] = inf_total_new_cal_pred[i]
        I_now_pred.iloc[i, scenario_number] = inf_total_cal_pred[i]
        sev_now_pred.iloc[i, scenario_number] = sev_cal_pred[i]


I_new_file = 'I_new_' + str(number_of_scenario) + 'scenario.csv'
I_now_file = 'I_now_' + str(number_of_scenario) + 'scenario.csv'
sev_now_file = 'sev_now_' + str(number_of_scenario) + 'scenario.csv'
I_new_pred.to_csv(I_new_file)
I_now_pred.to_csv(I_now_file)
sev_now_pred.to_csv(sev_now_file)

print('_______________________________')
print('end_____' + str(datetime.datetime.now().hour) + ':' + str(datetime.datetime.now().minute) + ':' + str(datetime.datetime.now().second))



