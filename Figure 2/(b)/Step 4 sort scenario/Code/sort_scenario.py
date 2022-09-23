
'''Calling library'''
import numpy as np
import pandas as pd
import datetime

'''Main'''

# Upload csv files
I_new_pred = pd.read_csv('I_new.csv')
I_now_pred = pd.read_csv('I_now.csv')

results = pd.DataFrame(index = (['sum_I_new', 'max_I_now']), columns = (I_new_pred.columns[1:]))

n_of_scenario = len(I_new_pred.columns) - 1

sum_of_I_new = I_new_pred[417:].sum()[1:]
max_of_I_now = I_now_pred[417:].max()[1:]

for scenario_number in range(n_of_scenario):
    results.iloc[0, scenario_number] = sum_of_I_new[scenario_number]
    results.iloc[1, scenario_number] = max_of_I_now[scenario_number]

sourcedata = 'sourcedata.csv'
results.to_csv(sourcedata)

v_sce_startday_pattern = 17
v_sce_startday_interval = 7 #days, -56 ~ 5/9 ~ +56
v_sce_speed_pattern = 17
v_sce_speed_interval = 50000 #shot/day, 232000 ~ 1032000
v_sce_maxrate_pattern = 11
v_sce_maxrate_interval = 5 #%, 40 ~ 100

for maxrate in range(v_sce_maxrate_pattern):
    table_sum = pd.DataFrame(index = range(v_sce_speed_pattern), columns = range(v_sce_startday_pattern))
    table_max = pd.DataFrame(index = range(v_sce_speed_pattern), columns = range(v_sce_startday_pattern))
    for y_axis in range(v_sce_speed_pattern):
        for x_axis in range(v_sce_startday_pattern):
            x_refer = x_axis * v_sce_speed_pattern * v_sce_maxrate_pattern + (v_sce_speed_pattern - 1 - y_axis) * v_sce_maxrate_pattern + maxrate
            table_sum.iloc[y_axis, x_axis] = results.iloc[0, x_refer]
            table_max.iloc[y_axis, x_axis] = results.iloc[1, x_refer]
    rate = 40 + maxrate * 5
    sum_file = 'sum_I_new_maxrate_' + str(rate) + '.csv'
    max_file = 'max_I_now_maxrate_' + str(rate) + '.csv'
    table_sum.to_csv(sum_file, index = False, header = False)
    table_max.to_csv(max_file, index = False, header = False)

