import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

# Use to switch directory paths to use (for the main or TESTING)
import os
currentDir = os.path.realpath('')
parentDir = os.path.realpath('..')



"""
ALWAYS SET THIS PATH TO THE ONE YOU WANT TO USE!
"""
desiredDir = currentDir


"""
Script plots:
- p0 reaction rates of REACLIB and my own
- p0 reaction rates of my own and AZURE2's
- my reaction rates (p0,p1,p2) vs AZURE2's reaction rates (p0,p1,p2)
"""




# Currently not being used but can be implemented, only need to get the coefficients
def rateFcn(a0,a1,a2,a3,a4,a5,a6,T):
    """ Take JINA REACLIB coefficients and calculate a reaction rate """
    return np.exp(a0+a1/T+a2/T**(1/3)+a3*T**(1/3)+a4*T+a5*T**(5/3)+a6*np.log(T))


# Plot the reaclib p0 rates
col = ["T9","Rate"]
sheetNames = ['nacr','il01','rath','laur','ths8','il10']
colors = ['b','g','r','c','m','y']
for ind,names in enumerate(sheetNames):
    pathToREACLIB = os.path.join(desiredDir,'24MgREACLIB.xlsx')  # This changes
    df = pd.read_excel(pathToREACLIB,sheet_name='%s'%names,header=0)
    _temp = df['T9'].values
    _rate = df['Rate'].values
    plt.plot(_temp,_rate,color=colors[ind],label=names)

p0Path = os.path.join(desiredDir,'legendre_out/DATA/p0/a0/p0_ratesEXTRAP.xlsx')
myp0rates = pd.read_excel(p0Path,sheet_name='Sheet1',header=0)

p1Path = os.path.join(desiredDir,'legendre_out/DATA/p1/a0/p1_ratesEXTRAP.xlsx')
myp1rates = pd.read_excel(p1Path,sheet_name='Sheet1',header=0)

p2Path = os.path.join(desiredDir,'legendre_out/DATA/p2/a0/p2_ratesEXTRAP.xlsx')
myp2rates = pd.read_excel(p2Path,sheet_name='Sheet1',header=0)

azp0Path = os.path.join(desiredDir,'p0rates.out')
azp0rates = pd.read_table(azp0Path,sep='\s+')

azp1Path = os.path.join(desiredDir,'p1rates.out')
azp1rates = pd.read_table(azp1Path,sep='\s+')

azp2Path = os.path.join(desiredDir,'p2rates.out')
azp2rates = pd.read_table(azp2Path,sep='\s+')



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
# plt.plot(aztemp0,azrate0,color='g',label='azure rate p0')
#
# plt.scatter(aztemp1,azrate1,c='fuchsia')
# plt.plot(aztemp1,azrate1,color='fuchsia',label='azure rate p1')
#
# plt.scatter(aztemp2,azrate2,c='b')
# plt.plot(aztemp2,azrate2,color='b',label='azure rate p2')


plt.yscale('log')
plt.ylim(1e-20,1e8) # 1e-18
# plt.xscale('log')
# plt.xlim(1e-1,1e1)
plt.xlim(0,1e1)
plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.title('p$_{0}$ Reaction Rates',fontsize=20)
plt.legend()
plt.grid(b=True, which='both', axis='both')
savePath = os.path.join(desiredDir,'p0CompareReactionRates.png')
plt.savefig(savePath,dpi=300)
# plt.show()


# Plot my p0 rates and Azures
plt.clf()

plt.plot(temp0,rate0,color='k',label='me')
plt.plot(aztemp0,azrate0,color='b',label='azure rate p0')

plt.yscale('log')
plt.ylim(1e-20,1e8) # 1e-18
# plt.xscale('log')
# plt.xlim(1e-1,1e1)
plt.xlim(0,1e1)
plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.title('p$_{0}$ Reaction Rates',fontsize=20)
plt.legend()
# plt.grid(b=True, which='both', axis='both')
savePath = os.path.join(desiredDir,'p0CompareReactionRatesAzure.png')
plt.savefig(savePath,dpi=300)


# Plot my p1 rates and Azures
plt.clf()

plt.plot(temp1,rate1,color='k',label='me')
plt.plot(aztemp1,azrate1,color='b',label='azure rate p1')

plt.yscale('log')
plt.ylim(1e-20,1e8) # 1e-18
# plt.xscale('log')
# plt.xlim(1e-1,1e1)
plt.xlim(0,1e1)
plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.title('p$_{1}$ Reaction Rates',fontsize=20)
plt.legend()
# plt.grid(b=True, which='both', axis='both')
savePath = os.path.join(desiredDir,'p1CompareReactionRatesAzure.png')
plt.savefig(savePath,dpi=300)


# Plot my p1 rates and Azures
plt.clf()

plt.plot(temp2,rate2,color='k',label='me')
plt.plot(aztemp2,azrate2,color='b',label='azure rate p2')

plt.yscale('log')
plt.ylim(1e-20,1e8) # 1e-18
# plt.xscale('log')
# plt.xlim(1e-1,1e1)
plt.xlim(0,1e1)
plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.title('p$_{2}$ Reaction Rates',fontsize=20)
plt.legend()
# plt.grid(b=True, which='both', axis='both')
savePath = os.path.join(desiredDir,'p2CompareReactionRatesAzure.png')
plt.savefig(savePath,dpi=300)



# Plot the rate ratios of azure wrt my p0
plt.clf()
y = azrate0/rate0
plt.plot(temp0,y,color='k')
# plt.yscale('log')
plt.xlim(0,10)
# plt.ylim(1e-2,1e1)
plt.ylim(0,3)
plt.ylabel('Ratio of Reaction Rate',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.title('p$_{0}$ Reaction Rates Ratio',fontsize=20)
eq = r'$\frac{Rate~p_{{0}_{AZURE2}}}{Rate~p_{0}}$'
plt.text(9, 1, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
# plt.grid(b=True, which='both', axis='both')

savePath = os.path.join(desiredDir,'RatioReactionRatesAzure.png')
plt.savefig(savePath,dpi=300)



# Plot the rate ratios wrt my p0
plt.clf()
y = (rate1+rate2)/rate0
plt.plot(temp0,y,color='k')
plt.yscale('log')
plt.xlim(0,10)
plt.ylim(1e-5,1e0)
plt.ylabel('Ratio of Reaction Rate',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.title('p$_{0}$ Reaction Rates Ratio',fontsize=20)
eq = r'$\frac{Rate~p_{1}~+~Rate~p_{2}}{Rate~p_{0}}$'
plt.text(9, 1e-3, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
# plt.grid(b=True, which='both', axis='both')

savePath = os.path.join(desiredDir,'RatioReactionRates.png')
plt.savefig(savePath,dpi=300)

# plt.show()
