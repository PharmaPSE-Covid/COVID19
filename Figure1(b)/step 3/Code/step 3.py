
'''Calling library'''
import numpy as np
import pandas as pd
import datetime

print('start  ---  ' + str(datetime.datetime.now().hour) + ':' + str(datetime.datetime.now().minute) + ':' + str(datetime.datetime.now().second))

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
beta_fm_exIandR = pd.read_csv('beta_fm_exIandR.csv')
beta_pred = pd.read_csv('beta_pred.csv')
v_scenario_1st_daily = pd.read_csv('vaccine_scenario_1st.csv')
v_scenario_2nd_daily = pd.read_csv('vaccine_scenario_2nd.csv')

# Time and the number of days

n_of_preddays = len(beta_pred)
diff_startday = int((pd.to_datetime(beta_pred.iloc[0][0]) - pd.to_datetime(beta_fm_exIandR.iloc[0]['date'])) / pd.Timedelta(days = 1))
till_inputend = n_of_preddays + diff_startday

dt = 1.0 / 10
days_divided = int((till_inputend - 1) / dt)

# Initial conditions
n_of_scenario = len(v_scenario_1st_daily.columns) - 1
timecounter = n_of_scenario // 20
print('number of scenario = ' + str(n_of_scenario))
print('_______________________________')

alpha = 0.40
gamma = 0.0691712976652005

effi_1st = 0.53
effi_2nd = 0.87

n_total = 1.263 * 1.0e+8
exp_ini =  0.1786
inf_ini = Spreading_data.iloc[0]['I_now']
rem_ini = Spreading_data.iloc[0]['R']

clm_I = v_scenario_1st_daily.columns[1:]
ind_I = np.arange(pd.to_datetime(beta_fm_exIandR.iloc[0]['date']), pd.to_datetime(beta_pred.iloc[-1][0]) + np.timedelta64(1,'D'), np.timedelta64(1,'D'), dtype='datetime64')
I_new_pred = pd.DataFrame(index = ind_I, columns = clm_I)
I_now_pred = pd.DataFrame(index = ind_I, columns = clm_I)

beta = np.zeros(days_divided)
day_counter = 0
for i in range(days_divided):
    if i % int(1.0 / dt) == 0:
        day_counter += 1
    if i < int((diff_startday) / dt):
        beta[i] = beta_fm_exIandR.iloc[day_counter - 1]['beta'] / n_total
    else:
        beta[i] = beta_pred.iloc[day_counter - diff_startday]['beta'] / n_total



for scenario_number in range(n_of_scenario):
    
    vac1_ini = v_scenario_1st_daily.iloc[0, scenario_number + 1]
    vac2_ini = v_scenario_2nd_daily.iloc[0, scenario_number + 1]
    sus_ini = n_total - exp_ini - inf_ini -rem_ini - vac1_ini - vac2_ini

    sus_old = sus_ini
    vac1_old = vac1_ini
    vac2_old = vac2_ini
    exp_old = exp_ini
    inf_old = inf_ini
    rem_old = rem_ini


    v_speed_1st = np.zeros(days_divided)
    v_speed_2nd = np.zeros(days_divided)

    day_counter = 0
    for i in range(days_divided):
        if i % int(1.0 / dt) == 0:
            day_counter += 1
        v_speed_1st[i] = v_scenario_1st_daily.iloc[day_counter - 1, scenario_number + 1]
        v_speed_2nd[i] = v_scenario_2nd_daily.iloc[day_counter - 1, scenario_number + 1]

    # Matrix for output
    sus = np.zeros(days_divided + 1)
    vac1 = np.zeros(days_divided + 1)
    vac2 = np.zeros(days_divided + 1)
    exp = np.zeros(days_divided + 1)
    inf = np.zeros(days_divided + 1)
    rem = np.zeros(days_divided + 1)
    inf_tot_cal = np.zeros(int(float(days_divided) * dt) + 1)
    inf_new_cal = np.zeros(int(float(days_divided) * dt) + 1)
    inf_now_cal = np.zeros(int(float(days_divided) * dt) + 1)
    sus[0] = sus_ini
    vac1[0] = vac1_ini
    vac2[0] = vac2_ini
    exp[0] = exp_ini
    inf[0] = inf_ini
    rem[0] = rem_ini

    # Numerical solution: Runge-Kutta method
    for i in range(days_divided):
        
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
            
        sus[i + 1] = sus_old
        vac1[i + 1] = vac1_old
        vac2[i + 1] = vac2_old
        exp[i + 1] = exp_old
        inf[i + 1] = inf_old
        rem[i + 1] = rem_old 
        
        i += 1
    
    k = 0
    for j in range(0, days_divided + 1):
        if j % int(1.0 / dt) == 0:
            inf_now_cal[k] = inf[j]
            inf_tot_cal[k] = (inf[j] + rem[j]) 
            k += 1                            
    

    I_new_pred.iloc[0, scenario_number] = inf_new_cal[0]
    I_now_pred.iloc[0, scenario_number] = inf_now_cal[0]

    for m in range(1, int((float(days_divided) * dt) + 1.0)):
        inf_new_cal[m] = inf_tot_cal[m] - inf_tot_cal[m - 1] 
        I_new_pred.iloc[m, scenario_number] = inf_new_cal[m]
        I_now_pred.iloc[m, scenario_number] = inf_now_cal[m]


I_new_file = 'I_new.csv'
I_now_file = 'I_now.csv'
I_new_pred.to_csv(I_new_file)
I_now_pred.to_csv(I_now_file)

print('_______________________________')
print('end_____' + str(datetime.datetime.now().hour) + ':' + str(datetime.datetime.now().minute) + ':' + str(datetime.datetime.now().second))





