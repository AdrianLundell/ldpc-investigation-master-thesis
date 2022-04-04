#%%
import numpy as np
ints = range(8)
gray_codes = {i^(i>>1) : i for i in ints}
voltage_levels = [0,2,4,6,8,10,12,14]

for i in ints:
    
    #Generate binary representation of i
    s = ''
    for j in range(3):
        s = str((i>>j) & 1) + s

    print(f"#{s} : Voltage level {gray_codes[i]}")
    print(voltage_levels[gray_codes[i]])
# %%
