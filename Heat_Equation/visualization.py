import math
import pandas as pd
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt


def sum(t, x, n):
    res = x*t*t
    for i in range(0, n):
        l = math.pi * (1 / 2 + i)
        res += 2 / (l**6) * (math.exp(-l*l*t) + l*l*t - 1) * math.cos(l*x)
    return res

# analytical_sol_50 = []
# for t in np.arange(0, 10.02, 0.01):
#     for x in np.linspace(0, 1, 50):
#         analytical_sol_50.append(sum(t, x, 100))

def read_blocks(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    blocks = content.strip().split('\n\n')
    all_data = []

    for block in blocks:
        lines = block.strip().split('\n')[1:]
        for line in lines:
            parts = line.split()
            if len(parts) > 1:
                all_data.append(float(parts[1]))
    return all_data

# Type file name you want to plot
data_50 = read_blocks("data50.txt")

def max_abs_difference(array1, array2):
    if len(array1) != len(array2):
        raise ValueError("Arrays must be equal length")
    differences = [abs(a - b) for a, b in zip(array1, array2)]
    return max(differences)

def get_max_error(nx, file):
    analytical_sol = []
    for t in np.arange(0, 10.02, 0.01):
        for x in np.linspace(0, 1, nx):
            analytical_sol.append(sum(t, x, 100))
    data = read_blocks(file)
    return max_abs_difference(analytical_sol, data)

max_errors = np.array([0.003863264830977986, 0.0013532650364793675, 0.0005965983013078713, 0.00046741730154131744, 0.0003506826468751001, 0.00024071412145509896])
h = [1 / 30, 1 / 50, 1/75, 1 / 85, 1 / 100, 1/125]
k = np.polyfit(np.log(h), np.log(max_errors), 1)

plt.scatter(np.log(h), np.log(max_errors))
plt.grid()
plt.xlabel("Log(step)")
plt.ylabel("Log(error)")
plt.show()

print(k)
