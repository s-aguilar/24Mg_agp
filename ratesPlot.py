import numpy as np
import pandas as pd
import matplotlib.pyplot as plt





# """
col = ["T9","Rate"]

myp0rates = pd.read_excel('legendre_out/DATA/analytically/p0/a0/p0_rates.xlsx',sheet_name='Sheet1',header=0)
myp1rates = pd.read_excel('legendre_out/DATA/analytically/p1/a0/p1_rates.xlsx',sheet_name='Sheet1',header=0)
myp2rates = pd.read_excel('legendre_out/DATA/analytically/p2/a0/p2_rates.xlsx',sheet_name='Sheet1',header=0)

azp0rates = pd.read_table('p0rates.out',sep='\s+')
azp1rates = pd.read_table('p1rates.out',sep='\s+')
azp2rates = pd.read_table('p2rates.out',sep='\s+')

# print(myp1rates.columns)
# print(azp1rates.columns)

temp0 = myp0rates['T9'].values
temp1 = myp1rates['T9'].values
temp2 = myp2rates['T9'].values
aztemp0 = azp0rates['T9'].values
aztemp1 = azp1rates['T9'].values
aztemp2 = azp2rates['T9'].values

rate0 = myp0rates['Rate'].values
rate1 = myp1rates['Rate'].values
rate2 = myp2rates['Rate'].values
azrate0 = azp0rates['Rate'].values
azrate1 = azp1rates['Rate'].values
azrate2 = azp2rates['Rate'].values

plt.scatter(temp0,rate0,c='k')
plt.plot(temp0,rate0,color='k',label='p$_0$')

plt.scatter(temp1,rate1,c='g')
plt.plot(temp1,rate1,color='g',label='p$_1$')

plt.scatter(temp2,rate2,c='b')
plt.plot(temp2,rate2,color='b',label='p$_2$')

# plt.scatter(aztemp0,azrate0,c='g')
# plt.plot(aztemp0,azrate0,color='g',label='azure rate p0')
#
# plt.scatter(aztemp1,azrate1,c='fuchsia')
# plt.plot(aztemp1,azrate1,color='fuchsia',label='azure rate p1')
#
# plt.scatter(aztemp2,azrate2,c='b')
# plt.plot(aztemp2,azrate2,color='b',label='azure rate p2')

# line = [1000]*len(aztemp1)
#
# plt.plot(aztemp1,line,color='y')
# plt.grid(True)

plt.ylim(1e-30,1e6)
plt.yscale('log')
plt.xlim(0,10)

plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)')
plt.xlabel('Temperature (T9)')
plt.legend()
plt.tight_layout()
plt.savefig('AzureFitCompareReactionRates.png',dpi=300)
plt.show()
