#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 

data = pd.read_csv("python-files/Thresholds.csv", sep=", ", header=0)
data["x1"] = data["sigma1"]**2/(data["sigma2"]**2) #Ratio
data["x2"] = 10*np.log10(2**2/(data["sigma1"]**2 + data["sigma2"]**2)) #Total SnR

#Plot 2D
fig1 = plt.figure()
ax = fig1.add_subplot()
ax.scatter(data.x1, data.middle, label="middle")
ax.scatter(data.x1, data.middle - data.neg_offset, label="negative")
ax.scatter(data.x1, data.middle + data.pos_offset, label = "positive")
ax.set_xlabel('Ratio')
ax.set_ylabel('Values')
ax.set_title("Ratio")
ax.legend()

fig2 = plt.figure()
ax = fig2.add_subplot()
ax.scatter(data.x2, data.middle)
ax.scatter(data.x2, data.middle - data.neg_offset)
ax.scatter(data.x2, data.middle + data.pos_offset)
ax.set_xlabel('Total noise')
ax.set_ylabel('Middle point')
ax.set_title('Total noise')

#Plot 2D
fig3 = plt.figure()
ax = fig3.add_subplot()
ax.scatter(data.x1, data.llr1)
ax.scatter(data.x1, data.llr2)
ax.scatter(data.x1, data.llr3)
ax.scatter(data.x1, data.llr4)
ax.set_xlabel('Ratio')
ax.set_ylabel('Values')
ax.set_title("LLRs")
ax.legend()

fig4 = plt.figure()
ax = fig4.add_subplot()
ax.scatter(data.x2, data.llr1)
ax.scatter(data.x2, data.llr2)
ax.scatter(data.x2, data.llr3)
ax.scatter(data.x2, data.llr4)
ax.set_xlabel('Total noise')
ax.set_ylabel('Values')
ax.set_title("LLRs")
ax.legend()

plt.show()

#Plot 3D
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(data.x1, data.x2, data.middle)
ax.set_xlabel('Ratio')
ax.set_ylabel('Total noise')
ax.set_zlabel('Middle point')
ax.set_title("Middle point")

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(data.x1, data.x2, data.neg_offset)
ax.set_xlabel('Ratio')
ax.set_ylabel('Total noise')
ax.set_zlabel('Negative offset')
ax.set_title("Negative offset")

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(data.x1, data.x2, data.pos_offset)
ax.set_xlabel('Ratio')
ax.set_ylabel('Total noise')
ax.set_zlabel('Positive offset')
ax.set_title("Positive offset")


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(data.x1, data.x2, data.llr1)
ax.scatter(data.x1, data.x2, data.llr2)
ax.scatter(data.x1, data.x2, data.llr3)
ax.scatter(data.x1, data.x2, data.llr4)

ax.set_xlabel('Ratio')
ax.set_ylabel('Total noise')
ax.set_zlabel('LLr')
ax.set_title("LLr")

plt.show()
