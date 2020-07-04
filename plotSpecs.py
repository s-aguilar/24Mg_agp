import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




# Read in the angle integrated cross-sections calculated via azure calculate w/o data
colNames = ['E','Counts']

# df = pd.read_table('E_cal_spectras/summedSpectrasALL_det-0',sep='\s+',names=colNames)
# dff = pd.read_table('E_cal_spectras/summedSpectrasBGsubALL_det-0',sep='\s+',names=colNames)


# 260,409 #238
# df = pd.read_table('E_cal_spectras/Spectra_run_260_det-0',sep='\s+',names=colNames)
# dff = pd.read_table('E_cal_spectras/BGSubSpectra_run_260_det-0',sep='\s+',names=colNames)

# df = pd.read_table('E_cal_spectras/summedSpectrasALL_det-0',sep='\s+',names=colNames)
# dff = pd.read_table('E_cal_spectras/summedSpectrasBGsubALL_det-0',sep='\s+',names=colNames)

df = pd.read_table('E_cal_spectras/summedSpectrasSingleResonance_det-0',sep='\s+',names=colNames)
dff = pd.read_table('E_cal_spectras/summedSpectrasBGsubSingleResonance_det-0',sep='\s+',names=colNames)

plt.scatter(df['E'],df['Counts'],c='k',s=1,label='Raw')
plt.scatter(dff['E'],dff['Counts'],c='r',s=1,label='BG Sub')
plt.xlim(0,2400)
plt.ylim(1,1e7)
# plt.yscale('log')
plt.yscale('log',basey=10)
plt.legend()
plt.show()
