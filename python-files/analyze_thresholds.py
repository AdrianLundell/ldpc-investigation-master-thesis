#%%
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("Thresholds.csv", sep=", ", header=0)
print(data)
#%%Histograms to see range of data
llrs = ["llr1", "llr2", "llr3", "llr4"]
thresholds = ["middle", "neg_offset", "pos_offset"]
data.hist(llrs)
data.hist(thresholds)

#%% Threshold-noise correlation'
data.plot.scatter("noise_db", "middle")
data.plot.scatter("noise_db", "neg_offset")
data.plot.scatter("noise_db", "pos_offset")

data.plot.scatter("sigma_ratio", "middle")
data.plot.scatter("sigma_ratio", "neg_offset")
data.plot.scatter("sigma_ratio", "pos_offset")

data.plot.scatter(["sigma_ratio", "sigma_ratio", "sigma_ratio", "sigma_ratio"], llrs)
data.plot.scatter(["noise_db", "noise_db", "noise_db", "noise_db"], llrs)


#%%MI-noise correlation
data.plot.scatter("noise_db", "max_mi")
data.plot.scatter("sigma_ratio", "max_mi")

# %%
