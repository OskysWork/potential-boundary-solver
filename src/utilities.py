import numpy as np
import os

# Setting working directory to script location
os.chdir(os.path.dirname(os.path.abspath(__file__)))

print("Now working in:", os.getcwd())

def read_aerofoil(file_name):
    """
    Read coords and retun Nx2 array
    """
    data = np.loadtxt(file_name, skiprows=1)
    return data[:,0], data[:, 1]