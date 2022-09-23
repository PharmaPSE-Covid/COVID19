
'''Calling library'''
import numpy as np
import pandas as pd
import datetime
import math

'''Main'''

# Upload csv files
beta_fm_exIandR = pd.read_csv('beta_fm_exIandR.csv')
beta_pred = pd.read_csv('beta_pred.csv')
sumple_scenario = pd.read_csv('param_values.csv')

print(sumple_scenario)
print(len(sumple_scenario))

# Time and the number of days
n_of_preddays = len(beta_pred)
diff_startday = int((pd.to_datetime(beta_pred.iloc[0][0]) - pd.to_datetime(beta_fm_exIandR.iloc[0]['date'])) / pd.Timedelta(days = 1))
till_inputend = n_of_preddays + diff_startday

n_total = 1.263 * 1.0e+8

ind_sce = np.arange(pd.to_datetime(beta_fm_exIandR.iloc[0]['date']), pd.to_datetime(beta_pred.iloc[-1][0]) + np.timedelta64(1,'D'), np.timedelta64(1,'D'), dtype='datetime64')

clm_sce = []
for i in range(len(sumple_scenario)):
    clm_sce.append('scenario_' + str(i))

v_speed_1st_daily = pd.DataFrame(index = ind_sce, columns = clm_sce)
v_speed_2nd_daily = pd.DataFrame(index = ind_sce, columns = clm_sce)

vaccine_1st_startday_base = 473 #2020/1/22 --> 2021/5/9
vaccine_1st_speed_base = 632000
vaccine_1st_maxrate_base = 80

for i in range(len(sumple_scenario)):
    if i % 100 == 0:
        print(str(i))
    
    startday = sumple_scenario.iloc[i]['startday']
    vaccine_1st_startday = vaccine_1st_startday_base + math.floor(startday)
    vaccine_2nd_startday = vaccine_1st_startday + 23
    
    cumulative_1st = 0
    cumulative_2nd = 0
    
    speed_1st = sumple_scenario.iloc[i]['speed']
    speed_2nd = speed_1st * 0.997

    maxrate = sumple_scenario.iloc[i]['maxrate']    
    max_1st = maxrate * n_total / 100
    max_2nd = max_1st * 0.975
    
    for day in range(till_inputend):
        if day < vaccine_1st_startday:
            v_speed_1st_daily.iloc[day, i] = 0
            v_speed_2nd_daily.iloc[day, i] = 0
        if day == vaccine_1st_startday:
            v_speed_1st_daily.iloc[day, i] = speed_1st * (startday - math.floor(startday))
            v_speed_2nd_daily.iloc[day, i] = 0
            cumulative_1st = cumulative_1st + speed_1st * (startday - math.floor(startday))
        if day > vaccine_1st_startday and day < vaccine_2nd_startday:
            v_speed_1st_daily.iloc[day, i] = speed_1st
            v_speed_2nd_daily.iloc[day, i] = 0
            cumulative_1st = cumulative_1st + speed_1st
        if day == vaccine_2nd_startday:
            v_speed_1st_daily.iloc[day, i] = speed_1st
            v_speed_2nd_daily.iloc[day, i] = speed_2nd * (startday - math.floor(startday))
            cumulative_1st = cumulative_1st + speed_1st
            cumulative_2nd = cumulative_2nd + speed_2nd * (startday - math.floor(startday))
        if day > vaccine_2nd_startday and cumulative_1st < max_1st:
            v_speed_1st_daily.iloc[day, i] = speed_1st
            v_speed_2nd_daily.iloc[day, i] = speed_2nd
            cumulative_1st = cumulative_1st + speed_1st
            cumulative_2nd = cumulative_2nd + speed_2nd
        if cumulative_1st >= max_1st and cumulative_2nd < max_2nd:
            v_speed_1st_daily.iloc[day, i] = 0
            v_speed_2nd_daily.iloc[day, i] = speed_2nd
            cumulative_2nd = cumulative_2nd + speed_2nd
        if cumulative_2nd >= max_2nd:
            v_speed_1st_daily.iloc[day, i] = 0
            v_speed_2nd_daily.iloc[day, i] = 0

scenario1_file = 'vaccine_scenario_1st.csv'
v_speed_1st_daily.to_csv(scenario1_file)

scenario2_file = 'vaccine_scenario_2nd.csv'
v_speed_2nd_daily.to_csv(scenario2_file)

