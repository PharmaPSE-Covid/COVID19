
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
param_values_dataframe = pd.DataFrame(param_values, columns = ['startday', 'speed', 'maxrate', 'efficacy'])
param_values_dataframe.to_csv('param_values.csv')
