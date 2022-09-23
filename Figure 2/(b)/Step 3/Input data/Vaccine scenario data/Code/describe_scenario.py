
'''Calling library'''
import numpy as np
import pandas as pd
import datetime

'''Main'''

# Upload csv files
beta_fm_exIandR = pd.read_csv('beta_fm_exIandR.csv')
beta_pred = pd.read_csv('beta_pred.csv')

# Time and the number of days
n_of_preddays = len(beta_pred)
diff_startday = int((pd.to_datetime(beta_pred.iloc[0][0]) - pd.to_datetime(beta_fm_exIandR.iloc[0]['date'])) / pd.Timedelta(days = 1))
till_inputend = n_of_preddays + diff_startday

n_total = 1.263 * 1.0e+8

ind_sce = np.arange(pd.to_datetime(beta_fm_exIandR.iloc[0]['date']), pd.to_datetime(beta_pred.iloc[-1][0]) + np.timedelta64(1,'D'), np.timedelta64(1,'D'), dtype='datetime64')
v_sce_startday_pattern = 17
v_sce_startday_interval = 7 #days, -56 ~ 5/9 ~ +56
v_sce_speed_pattern = 17
v_sce_speed_interval = 50000 #shot/day, 232000 ~ 1032000
v_sce_maxrate_pattern = 11
v_sce_maxrate_interval = 5 #%, 40 ~ 100

clm_sce = []
for i in range(v_sce_startday_pattern):
    for j in range(v_sce_speed_pattern):
        for k in range(v_sce_maxrate_pattern):
            startday = (i - 8) * v_sce_startday_interval
            speed = 2 + j * 0.5
            maxrate = 40 + k * v_sce_maxrate_interval
            clm_sce.append(str(startday) + '_' + str(speed) + '_' + str(maxrate))

v_speed_1st_daily = pd.DataFrame(index = ind_sce, columns = clm_sce)
v_speed_2nd_daily = pd.DataFrame(index = ind_sce, columns = clm_sce)

vaccine_1st_startday_base = 473 #2020/1/22 --> 2021/5/9
vaccine_1st_speed_base = 632000
vaccine_1st_maxrate_base = 80

for i in range(v_sce_startday_pattern):
    vaccine_1st_startday = vaccine_1st_startday_base + (i -8) * v_sce_startday_interval
    vaccine_2nd_startday = vaccine_1st_startday + 23
    print(str(i) + '_' + str(datetime.datetime.now().hour) + ':' + str(datetime.datetime.now().minute) + ':' + str(datetime.datetime.now().second))
    for j in range(v_sce_speed_pattern):
        shot_per_day_1st = vaccine_1st_speed_base + (j - 8) * v_sce_speed_interval
        shot_per_day_2nd = shot_per_day_1st * 0.997
        for k in range(v_sce_maxrate_pattern):
            maxrate_1st = vaccine_1st_maxrate_base + (k - 8) * v_sce_maxrate_interval
            maxrate_2nd = maxrate_1st * 0.975
            
            period_1st = (n_total * maxrate_1st / 100) // shot_per_day_1st
            period_2nd = (n_total * maxrate_2nd / 100) // shot_per_day_2nd
            vaccine_1st_endday = vaccine_1st_startday + period_1st
            vaccine_2nd_endday = vaccine_2nd_startday + period_2nd
            scenario_number = k + v_sce_maxrate_pattern * j + v_sce_speed_pattern * v_sce_maxrate_pattern * i

            for day in range(till_inputend):
                if day < vaccine_1st_startday:
                    v_speed_1st_daily.iloc[day, scenario_number] = 0
                    v_speed_2nd_daily.iloc[day, scenario_number] = 0
                if day >= vaccine_1st_startday and day < vaccine_2nd_startday:
                    v_speed_1st_daily.iloc[day, scenario_number] = shot_per_day_1st
                    v_speed_2nd_daily.iloc[day, scenario_number] = 0
                if day >= vaccine_2nd_startday and day <= vaccine_1st_endday:
                    v_speed_1st_daily.iloc[day, scenario_number] = shot_per_day_1st
                    v_speed_2nd_daily.iloc[day, scenario_number] = shot_per_day_2nd
                if day > vaccine_1st_endday and day <= vaccine_2nd_endday:
                    v_speed_1st_daily.iloc[day, scenario_number] = 0
                    v_speed_2nd_daily.iloc[day, scenario_number] = shot_per_day_2nd
                if day > vaccine_2nd_endday:
                    v_speed_1st_daily.iloc[day, scenario_number] = 0
                    v_speed_2nd_daily.iloc[day, scenario_number] = 0


scenario1_file = 'vaccine_scenario_1st.csv'
v_speed_1st_daily.to_csv(scenario1_file)

scenario2_file = 'vaccine_scenario_2nd.csv'
v_speed_2nd_daily.to_csv(scenario2_file)


