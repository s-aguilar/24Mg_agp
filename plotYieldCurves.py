import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Read in the data into dataframe
colNames = ['Energy','Angle','Yield','Error']
p1 = pd.read_table('rMatrix/24Mg_rMatrix_p1_angInt.dat',names=colNames)
p2 = pd.read_table('rMatrix/24Mg_rMatrix_p2_angInt.dat',names=colNames)
a1 = pd.read_table('rMatrix/24Mg_rMatrix_a1_angInt.dat',names=colNames)
o17_1 = pd.read_table('rMatrix/24Mg_rMatrix_o17_1_angInt.dat',names=colNames)
o17_ng = pd.read_table('rMatrix/24Mg_rMatrix_o17_ng_angInt.dat',names=colNames)


# Get the yields
p1_nrg = p1['Energy'].to_numpy()
p1_y = p1['Yield'].to_numpy()
p1_yerr = p1['Error'].to_numpy()

p2_nrg = p2['Energy'].to_numpy()
p2_y = p2['Yield'].to_numpy()
p2_yerr = p2['Error'].to_numpy()

a1_nrg = a1['Energy'].to_numpy()
a1_y = a1['Yield'].to_numpy()
a1_yerr = a1['Error'].to_numpy()

o17_1_nrg = o17_1['Energy'].to_numpy()
o17_1_y = o17_1['Yield'].to_numpy()
o17_1_yerr = o17_1['Error'].to_numpy()

o17_ng_nrg = o17_ng['Energy'].to_numpy()
o17_ng_y = o17_ng['Yield'].to_numpy()
o17_ng_yerr = o17_ng['Error'].to_numpy()


# Plot yield curves
# plt.errorbar(x = p1_nrg, y = p1_y, yerr = p1_yerr, fmt = 'k.', markersize = '4', label = 'p1')
# plt.errorbar(x = p2_nrg, y = p2_y, yerr = p2_yerr, fmt = 'b.', markersize = '4', label = 'p2')
plt.errorbar(x = a1_nrg, y = a1_y, yerr = a1_yerr, fmt = 'm.', markersize = '4', label = 'a1')
plt.errorbar(x = o17_1_nrg, y = o17_1_y, yerr = o17_1_yerr, fmt = 'r.', markersize = '4', label = 'o17_1')
plt.errorbar(x = o17_ng_nrg, y = o17_ng_y, yerr = o17_ng_yerr, fmt = 'c.', markersize = '4', label = 'o17_ng')
plt.xlabel('Energy lab (MeV)')
plt.ylabel('Yield (arb. units)')
plt.yscale('log')
plt.title('Yield Curve')
plt.legend()
plt.show()
