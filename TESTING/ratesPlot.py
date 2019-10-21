import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

def rateFcn(a0,a1,a2,a3,a4,a5,a6,T):
    """ Take JINA REACLIB coefficients and calculate a reaction rate """
    return np.exp(a0+a1/T+a2/T**(1/3)+a3*T**(1/3)+a4*T+a5*T**(5/3)+a6*np.log(T))



# """
col = ["T9","Rate"]

sheetNames = ['nacr','il01','rath','laur','ths8','il10']
colors = ['b','g','r','c','m','y']

for ind,names in enumerate(sheetNames):
    df = pd.read_excel('/Users/sebastian/Desktop/24Mg_agp/24MgREACLIB.xlsx',sheet_name='%s'%names,header=0)
    _temp = df['T9'].values
    _rate = df['Rate'].values
    col = colors[ind]

    plt.plot(_temp,_rate,color=col,label=names)
    plt.ylim(1e-18,1e8)
    plt.yscale('log')
    plt.xlim(0,10)

    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    # plt.legend()
    plt.tight_layout()

# plt.show()
# exit()
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

# plt.scatter(temp0,rate0,c='k')
plt.plot(temp0,rate0,color='k',label='me')

# plt.scatter(temp1,rate1,c='r')
# plt.plot(temp1,rate1,color='r',label='$$p$_{1}$')

# plt.scatter(temp2,rate2,c='b')
# plt.plot(temp2,rate2,color='b',label='$\\gamma$p$_2$')

# plt.scatter(aztemp0,azrate0,c='g')
# plt.plot(aztemp0,azrate0,color='k',label='azure rate p0')
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

# plt.ylim(1e-18,1e4)
# plt.yscale('log')
# plt.xlim(0,10)
#
# plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
# plt.xlabel('Temperature (T9)',fontsize=14)
# plt.legend()
# plt.tight_layout()
plt.title('p$_{0}$ Reaction Rates',fontsize=20)
plt.legend()
plt.xscale('log') #
plt.xlim(1e-1,1e1)
plt.savefig('p0CompareReactionRates.png',dpi=300)
# plt.show()
exit()

plt.clf()
# Plot the rate ratios
y = (rate1+rate2)/rate0
plt.plot(temp0,y)

plt.yscale('log')
plt.xlim(0,10)
plt.ylim(1e-5,1e0)
plt.ylabel('Ratio of Reaction Rate',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.tight_layout()
plt.savefig('RatioReactionRates.png',dpi=300)
plt.show()
