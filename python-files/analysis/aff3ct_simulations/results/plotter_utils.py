from distutils.command.config import config
import numpy as np
import pandas as pd
from config import cfg

def trim_data(lines):
    data = []
    #Only append data lines
    for line in lines:
        if line.startswith(' '):
            line = line.strip()
            line = line.split("|")
            data.append(line)
    
    # Further cleanup
    data = np.array(data, dtype=str)
    data = np.char.strip(data,' ')
    data = data[:,:-1]
    bool_data = np.invert(np.isin(data,['']))[0]
    data = data[:,bool_data]
    data = data.astype(np.float64)

    return data

def create_dataframe(data):
    if data.shape[1] == 8:
        columns = ["Es/N0", "Eb/N0", "FRA", "BE", "FE", "BER", "FER", "SIM_THR"]
    elif data.shape[1] == 7:
        columns = ["RBER", "FRA", "BE", "FE", "BER", "FER", "SIM_THR"]
    df = pd.DataFrame(data,columns=columns)
    return df

# Returns a dict with file name as key and corresponding data (data frame) as value 
def read_data(input_files):
    data_dict = {}
    for fname in input_files:
        with open(fname,'r') as f:
            lines = f.readlines()
        
        data = trim_data(lines)
        df = create_dataframe(data)
        data_dict[fname] = df

    return data_dict

        
