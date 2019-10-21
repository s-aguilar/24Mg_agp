import numpy as np

# rutherford theoretical
ang = 145*np.pi/180
unit = 1.6E-19
eps = 8.85418782E-12
KE_a = 5.13E6*unit
term1 = 2*79*unit*unit/(4*np.pi*eps)
term2 = 1/(4*KE_a)

diff_cross = term1*term1*term2*term2/(np.sin(ang/2))**4
cross = diff_cross
# print(cross)



n_I = 38422E-10/(2.*unit)

avo = 6.022E23
mass = 197/avo
thicc = 103E-6/0.0001
area = np.pi*(.005)**2

n_T = (thicc/mass)

_yield = 1538
eff = _yield/(cross*n_I*n_T)
# print(eff)



# Experimental
n_I = 969778E-10/(2.*unit)
_yield = 2359

# steradian = np.pi*.0025**2/.6**2
# print(steradian)
cross = 0.05E-28#*4.5E-5
# print(cross)

mass1 = 24/avo
calc_thick = (_yield/(eff*cross*n_I))*mass1
print(calc_thick*1E6/1E4)
