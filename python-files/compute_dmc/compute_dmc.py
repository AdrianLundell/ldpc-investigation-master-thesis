#%%
import yaml
import numpy as np 
import dmc_utils

#Load settings
try:
    with open("config.yml", "r") as ymlfile:
        settings = yaml.safe_load(ymlfile)
except:
    print("Config not found, using default settings.")
    settings = {
                "n_thresholds" : 3,
                "skew" : 0.5,
                "rber_range" : [0.001, 0.5],
                "n_rber" : 10,
                "output" : "4_05_AWGN_test.csv",
                "symmetric_thresholds" : True
                }

try:
    rber_arr = np.linspace(settings["rber_range"][0], settings["rber_range"][1], settings["n_rber"], False)
    n_columns = 5 + (1 + 2*settings["n_thresholds"])
    result = np.zeros((rber_arr.size, n_columns))

    print("Calculating thresholds...")

    for i, rber in enumerate(rber_arr):
        sigma = dmc_utils.rber_to_sigma(rber, settings["skew"])
        sigma1 = (1-settings["skew"])*sigma
        sigma2 = settings["skew"]*sigma

        if settings["n_thresholds"] == 3:
            t, capacity = dmc_utils.optimize_thresholds(sigma1, sigma2, symmetric=settings["symmetric_thresholds"])
            llr = dmc_utils.llrs(t, sigma1, sigma2)
        elif settings["n_thresholds"] == 1:
            t = dmc_utils.mid_point(sigma1, sigma2)
            capacity = 0
            llr = dmc_utils.llrs([t], sigma1, sigma2)


        row = np.hstack((rber, sigma, sigma1, sigma2, capacity, t, llr))
        result[i, :] = row

        status = f"{i+1}/{settings['n_rber']} rbers computed.  "
        print(status, end = "\r", flush=True)
except:
    print("Threshold calculation failed, writing completed values to file.")
finally:
    header = "RBER, SIGMA, SIGMA1, SIGMA2, CAPACITY, "
    for i in range(settings["n_thresholds"]):
        header += f"T{i}, "
    for i in range(settings["n_thresholds"]+1):
        header += f"LLR{i}, "
    header = header[:-2]

    np.savetxt(settings["output"], result, header = header)

    print(f"Computation done, file saved to {settings['output']}.")
# %%
