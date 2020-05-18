import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




# Read in the angle integrated cross-sections calculated via azure calculate w/o data
colNames = ['E','Counts']
# 260,409
df = pd.read_table('E_cal_spectras/Spectra_run_260_det-0',sep='\s+',names=colNames)
dff = pd.read_table('E_cal_spectras/BGSubSpectra_run_260_det-0',sep='\s+',names=colNames)

plt.scatter(df['E'],df['Counts'],c='k',s=1)
plt.scatter(dff['E'],dff['Counts'],c='r',s=1)
plt.xlim(0,1800)
plt.ylim(1,1e6)
# plt.yscale('log')
plt.yscale('log',basey=10)
plt.show()
