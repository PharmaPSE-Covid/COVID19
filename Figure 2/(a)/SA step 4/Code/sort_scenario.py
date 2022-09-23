
'''Calling library'''
import numpy as np
import pandas as pd
import datetime

'''Main'''

# Upload csv files
I_new_pred = pd.read_csv('I_new_10000scenario.csv')
I_now_pred = pd.read_csv('I_now_10000scenario.csv')

results = pd.DataFrame(index = (I_new_pred.columns[1:]), columns = (['sum_I_new', 'max_I_now']))

n_of_scenario = len(I_new_pred.columns) - 1

sum_of_I_new = I_new_pred[417:].sum()[1:]
max_of_I_now = I_now_pred[417:].max()[1:]

for scenario_number in range(n_of_scenario):
    results.iloc[scenario_number, 0] = sum_of_I_new[scenario_number]
    results.iloc[scenario_number, 1] = max_of_I_now[scenario_number]

sourcedata = 'run_model_gsa.csv'
results.to_csv(sourcedata)
