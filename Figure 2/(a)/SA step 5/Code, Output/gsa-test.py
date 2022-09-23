
from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
import numpy as np
import pandas as pd

# Define the model inputs
problem = {
    'num_vars': 4,
    'names': ['startday', 'speed', 'maxrate', 'efficacy'],
    'bounds': [[-28, 28],
               [432000, 832000],
               [30, 60],
               [0.57, 0.97]]
}

# Generate samples
param_values = saltelli.sample(problem, 1000)

# Run model

run_model_gsa = pd.read_csv('run_model_gsa.csv')

Y1 = run_model_gsa['sum_I_new'].values
Y2 = run_model_gsa['max_I_now'].values

# Perform analysis
print('------------------------------------------')
print('About sum of I_new')
Si = sobol.analyze(problem, Y1, print_to_console=True)
print('------------------------------------------')
print('About max of I_now')
Si = sobol.analyze(problem, Y2, print_to_console=True)

