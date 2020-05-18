import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from AtomicMassTable import GetElement


# Read in the data into dataframe
colNames = ['E','Ang','XS','Err']

# GetElement returns (index,Z,A,Mass,'element')
_m1H = GetElement(1,1)[3]
_m4He = GetElement(2,4)[3]
_m24Mg = GetElement(12,24)[3]
_m27Al = GetElement(13,27)[3]

# Separation energies (MeV), used to convert to excitation energy
s_a = 9.98414
s_p = 11.5849


# In Lab energy
ap1 = pd.read_table('rMatrix/24Mg_rMatrix_p1.dat',names=colNames)
ap2 = pd.read_table('rMatrix/24Mg_rMatrix_p2.dat',names=colNames)
aa1 = pd.read_table('rMatrix/24Mg_rMatrix_a1.dat',names=colNames)

nel_pp1 = pd.read_table('rMatrix/Nelson_pp1_EXFOR.dat',names=colNames)
nel_pp2 = pd.read_table('rMatrix/Nelson_pp2_EXFOR.dat',names=colNames)
nel_pa1 = pd.read_table('rMatrix/Nelson_pa1_EXFOR.dat',names=colNames)


# Conversion factor: lab to CM energy
labToCM_me = (_m24Mg/(_m24Mg+_m4He))
labToCM_nel = (_m27Al/(_m27Al+_m1H))

# Lab to excitation energy
ap1['E'] = ap1['E']*labToCM_me+s_a
ap2['E'] = ap2['E']*labToCM_me+s_a
aa1['E'] = aa1['E']*labToCM_me+s_a

nel_pp1['E'] = nel_pp1['E']*labToCM_nel+s_p
nel_pp2['E'] = nel_pp2['E']*labToCM_nel+s_p
nel_pa1['E'] = nel_pa1['E']*labToCM_nel+s_p

# Select angle
angle = 60
ap1 = ap1.query('Ang == @angle')
ap2 = ap2.query('Ang == @angle')
aa1 = aa1.query('Ang == @angle')

# angle = 135
# nel_pp1 = nel_pp1.query('Ang == @angle')
# nel_pp2 = nel_pp2.query('Ang == @angle')
# nel_pa1 = nel_pa1.query('Ang == @angle')





# Plot Excitation Curves
plt.errorbar(x = ap1['E'], y = ap1['XS'], yerr = ap1['Err'], fmt = 'r.', markersize = '4', label = 'ap1')
# plt.errorbar(x = ap2['E'], y = ap2['XS'], yerr = ap2['Err'], fmt = 'g.', markersize = '4', label = 'ap2')
# plt.errorbar(x = aa1['E'], y = aa1['XS'], yerr = aa1['Err'], fmt = 'b.', markersize = '4', label = 'aa1')

plt.errorbar(x = nel_pp1['E'], y = nel_pp1['XS'], yerr = nel_pp1['Err'], fmt = 'k.', markersize = '4', label = 'nel_pp1')
# plt.errorbar(x = nel_pp2['E'], y = nel_pp2['XS'], yerr = nel_pp2['Err'], fmt = 'm.', markersize = '4', label = 'nel_pp2')
# plt.errorbar(x = nel_pa1['E'], y = nel_pa1['XS'], yerr = nel_pa1['Err'], fmt = 'c.', markersize = '4', label = 'nel_pa1')

plt.xlabel('Energy(MeV)')
plt.ylabel('XS (barns)')
plt.yscale('log')
plt.title('Excitation Curve')
plt.legend()
plt.show()
