import numpy as np
import pandas as pd
import matplotlib.pyplot as plt





# """
col = ["T9","Rate"]

# myp0rates = pd.read_excel('legendre_out/DATA/p0/a0/p0_rates.xlsx',sheet_name='Sheet1',header=0)
myp0rates = pd.read_excel('24MgREACLIB.xlsx',sheet_name='il10',header=0)
myp1rates = pd.read_excel('legendre_out/DATA/p1/a0/p1_rates.xlsx',sheet_name='Sheet1',header=0)
# myp1rates = myp1rates.query('T9 != 4.5')
myp2rates = pd.read_excel('legendre_out/DATA/p2/a0/p2_rates.xlsx',sheet_name='Sheet1',header=0)
# myp2rates = myp2rates.query('T9 != 4.5')

# azp0rates = pd.read_table('p0rates.out',sep='\s+')
# azp1rates = pd.read_table('p1rates.out',sep='\s+')
# azp2rates = pd.read_table('p2rates.out',sep='\s+')

# print(myp1rates.columns)
# print(azp1rates.columns)

temp0 = myp0rates['T9'].to_numpy()
temp1 = myp1rates['T9'].to_numpy()
temp2 = myp2rates['T9'].to_numpy()

# aztemp0 = azp0rates['T9'].to_numpy()
# aztemp1 = azp1rates['T9'].to_numpy()
# aztemp2 = azp2rates['T9'].to_numpy()

rate0 = myp0rates['Rate'].to_numpy()
rate1 = myp1rates['Rate'].to_numpy()
rate2 = myp2rates['Rate'].to_numpy()
# azrate0 = azp0rates['Rate'].to_numpy()
# azrate1 = azp1rates['Rate'].to_numpy()
# azrate2 = azp2rates['Rate'].to_numpy()

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

plt.ylim(1e-30,1e8)
plt.yscale('log')
plt.xlim(0,10)

plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)')
plt.xlabel('Temperature (T9)')
plt.legend()
plt.tight_layout()
plt.savefig('compareReactionRates.png',dpi=300)
plt.show()

# myp0rates = myp0rates.query('T9 in @temp1')
# # print(myp0rates.head(24))
# temp0 = myp0rates['T9'].to_numpy()
# rate0 = myp0rates['Rate'].to_numpy()

print(temp0,'\n',temp1)
###
plt.clf()
plt.scatter(temp1,rate1/rate0,c='g')
plt.plot(temp1,rate1/rate0,color='g',label='p$_1$ ratio')

plt.scatter(temp2,rate2/rate0,c='b')
plt.plot(temp2,rate2/rate0,color='b',label='p$_2$ ratio')

plt.scatter(temp2,(rate1+rate2)/rate0,c='k')
plt.plot(temp2,(rate1+rate2)/rate0,color='k',label='p$_1$ + p$_2$ ratio')
plt.ylim(1e-9,1e0)
plt.yscale('log')
plt.xlim(0,10)

plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)')
plt.xlabel('Temperature (T9)')
plt.legend()
plt.tight_layout()
plt.savefig('rxnRateRatios.png',dpi=300)
plt.show()


print('p1 / p0 ratio')
for ind in range(len(rate0)):
    if(temp0[ind]>=1 and temp0[ind]<=10):
        print('T = %.2f GK\t\t%.2E'%(temp0[ind],rate1[ind]/rate0[ind]))

print('\np2 / p0 ratio')
for ind in range(len(rate0)):
    if(temp0[ind]>=1 and temp0[ind]<=10):
        print('T = %.2f GK\t\t%.2E'%(temp0[ind],rate2[ind]/rate0[ind]))

print('\np1+p2 / p0 ratio')
for ind in range(len(rate0)):
    if(temp0[ind]>=1 and temp0[ind]<=10):
        print('T = %.2f GK\t\t%.2E'%(temp0[ind],(rate1[ind]+rate2[ind])/rate0[ind]))
