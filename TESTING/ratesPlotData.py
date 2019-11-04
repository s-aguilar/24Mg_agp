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
# for ind,names in enumerate(sheetNames):
#     pathToREACLIB = os.path.join(desiredDir,'24MgREACLIB.xlsx')  # This changes
#     df = pd.read_excel(pathToREACLIB,sheet_name='%s'%names,header=0)
#     _temp = df['T9'].values
#     _rate = df['Rate'].values
#     plt.plot(_temp,_rate,color=colors[ind],label=names)

p1Path = os.path.join(desiredDir,'legendre_out/DATA/p1/a0/p1_rates.xlsx')
myp1rates = pd.read_excel(p1Path,sheet_name='Sheet1',header=0)

p2Path = os.path.join(desiredDir,'legendre_out/DATA/p2/a0/p2_rates.xlsx')
myp2rates = pd.read_excel(p2Path,sheet_name='Sheet1',header=0)

azp0Path = os.path.join(desiredDir,'p0rates.out')
azp0rates = pd.read_table(azp0Path,sep='\s+')

azp1Path = os.path.join(desiredDir,'p1rates.out')
azp1rates = pd.read_table(azp1Path,sep='\s+')

azp2Path = os.path.join(desiredDir,'p2rates.out')
azp2rates = pd.read_table(azp2Path,sep='\s+')



temp1 = myp1rates['T9'].values
temp2 = myp2rates['T9'].values
aztemp0 = azp0rates['T9'].values
aztemp1 = azp1rates['T9'].values
aztemp2 = azp2rates['T9'].values


rate1 = myp1rates['Rate'].values
rate2 = myp2rates['Rate'].values
azrate0 = azp0rates['Rate'].values
azrate1 = azp1rates['Rate'].values
azrate2 = azp2rates['Rate'].values

plt.scatter(temp1,rate1,c='r')
plt.plot(temp1,rate1,color='r',label='$\\gamma$p$_{1}$')

plt.scatter(temp2,rate2,c='g')
plt.plot(temp2,rate2,color='g',label='$\\gamma$p$_2$')

plt.scatter(aztemp0,azrate0,c='b')
plt.plot(aztemp0,azrate0,color='b',label='azure rate p0')

plt.scatter(aztemp1,azrate1,c='fuchsia')
plt.plot(aztemp1,azrate1,color='fuchsia',label='azure rate p1')

plt.scatter(aztemp2,azrate2,c='k')
plt.plot(aztemp2,azrate2,color='k',label='azure rate p2')


plt.yscale('log')
plt.ylim(1e-20,1e8) # 1e-18
# plt.xscale('log')
# plt.xlim(1e-1,1e1)
plt.xlim(0,1e1)
plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.title('Reaction Rates',fontsize=20)
plt.legend()
# plt.grid(b=True, which='both', axis='both')
savePath = os.path.join(desiredDir,'DataCompareReactionRates.png')
plt.savefig(savePath,dpi=300)
# plt.show()



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
savePath = os.path.join(desiredDir,'p1DataCompareReactionRatesAzure.png')
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
savePath = os.path.join(desiredDir,'p2DataCompareReactionRatesAzure.png')
plt.savefig(savePath,dpi=300)



# Plot the rate ratios of my p1 and p2 wrt iliadis p0
plt.clf()

pathToREACLIB = os.path.join(desiredDir,'24MgREACLIB.xlsx')
df = pd.read_excel(pathToREACLIB,sheet_name='il10',header=0)
_tempIliadis = df['T9'].values
_rateIliadis = df['Rate'].values

df = pd.read_excel(pathToREACLIB,sheet_name='il10_coef',header=0)
_a0 = df['a0'].values
_a1 = df['a1'].values
_a2 = df['a2'].values
_a3 = df['a3'].values
_a4 = df['a4'].values
_a5 = df['a5'].values
_a6 = df['a6'].values

fname = 'rate_Temps.dat'
temperature = np.loadtxt( fname )

colors = ['r','g','b']
for ind in range(len(_a0)):
    _rate_ = rateFcn(_a0[ind],_a1[ind],_a2[ind],_a3[ind],_a4[ind],_a5[ind],_a6[ind],temperature)
    plt.plot(temperature,_rate_,c=colors[ind])

# plt.scatter(_temp,_rate,c='k')
plt.plot(_tempIliadis,_rateIliadis,c='k',label='Iliadis Points')
# plt.ylim(1e-20,1e8) # 1e-18
plt.yscale('log')
plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.legend()
savePath = os.path.join(desiredDir,'DataRatioIliadis.png')
plt.savefig(savePath,dpi=300)
# plt.show()

plt.clf()

for ind in range(len(_a0)):
    _rate_ = rateFcn(_a0[ind],_a1[ind],_a2[ind],_a3[ind],_a4[ind],_a5[ind],_a6[ind],temperature)
    plt.plot(temperature,_rate_,c=colors[ind])

# plt.scatter(_temp,_rate,c='k')
plt.plot(_tempIliadis,_rateIliadis,c='k',label='Iliadis Points')
plt.ylim(1e-20,1e8) # 1e-18
plt.yscale('log')
plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.legend()
savePath = os.path.join(desiredDir,'DataRatioIliadis1.png')
plt.savefig(savePath,dpi=300)
# plt.show()

plt.clf()
# exit()


ratiop1 = (rate1)/_rateIliadis
plt.plot(_tempIliadis,ratiop1,color='k')
plt.yscale('log')
plt.xlim(0,10)
plt.ylim(1e-5,1e0)
plt.ylabel('Ratio of Reaction Rate',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
# plt.title('p$_{0}$ Reaction Rates Ratio',fontsize=20)
eq = r'$\frac{ Rate~p_{1_{Data}} }{ Rate~p_{0_{Iliadis}} }$'
plt.text(9, 1e-3, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
# plt.grid(b=True, which='both', axis='both')

savePath = os.path.join(desiredDir,'DataRatiop1ReactionRates.png')
plt.savefig(savePath,dpi=300)

df = pd.DataFrame(data=_tempIliadis,columns=['T9'])
df = df.assign(Ratio=pd.Series(ratiop1,index=df.index).values)
csvPath = os.path.join(desiredDir,'iliadisRatiop1.csv')
df.to_csv(csvPath)
excelPath = os.path.join(desiredDir,'iliadisRatiop1.xlsx')
df.to_excel(excelPath)

plt.clf()


ratiop2 = (rate2)/_rateIliadis
plt.plot(_tempIliadis,ratiop2,color='k')
plt.yscale('log')
plt.xlim(0,10)
plt.ylim(1e-5,1e0)
plt.ylabel('Ratio of Reaction Rate',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
# plt.title('p$_{0}$ Reaction Rates Ratio',fontsize=20)
eq = r'$\frac{ Rate~p_{2_{Data}} }{ Rate~p_{0_{Iliadis}} }$'
plt.text(9, 1e-3, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
# plt.grid(b=True, which='both', axis='both')

savePath = os.path.join(desiredDir,'DataRatiop2ReactionRates.png')
plt.savefig(savePath,dpi=300)

df = pd.DataFrame(data=_tempIliadis,columns=['T9'])
df = df.assign(Ratio=pd.Series(ratiop2,index=df.index).values)
csvPath = os.path.join(desiredDir,'iliadisRatiop2.csv')
df.to_csv(csvPath)
excelPath = os.path.join(desiredDir,'iliadisRatiop2.xlsx')
df.to_excel(excelPath)

plt.clf()


ratiosum = (rate1+rate2)/_rateIliadis
plt.plot(_tempIliadis,ratiosum,color='k')
plt.yscale('log')
plt.xlim(0,10)
plt.ylim(1e-5,1e0)
plt.ylabel('Ratio of Reaction Rate',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
# plt.title('p$_{0}$ Reaction Rates Ratio',fontsize=20)
eq = r'$\frac{ Rate~p_{1_{Data}}~+~Rate~p_{2_{Data}} }{ Rate~p_{0_{Iliadis}} }$'
plt.text(9, 1e-3, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
# plt.grid(b=True, which='both', axis='both')

savePath = os.path.join(desiredDir,'DataRatioSumReactionRates.png')
plt.savefig(savePath,dpi=300)

df = pd.DataFrame(data=_tempIliadis,columns=['T9'])
df = df.assign(Ratio=pd.Series(ratiosum,index=df.index).values)
csvPath = os.path.join(desiredDir,'iliadisRatioSum.csv')
df.to_csv(csvPath)
excelPath = os.path.join(desiredDir,'iliadisRatioSum.xlsx')
df.to_excel(excelPath)

plt.clf()


plt.plot(_tempIliadis,ratiosum,color='k',label='Sum')
plt.plot(_tempIliadis,ratiop1,color='r',label='p1')
plt.plot(_tempIliadis,ratiop2,color='g',label='p2')
plt.yscale('log')
plt.xlim(0,10)
plt.ylim(1e-5,1e0)
plt.ylabel('Ratios of Reaction Rate',fontsize=14)
plt.xlabel('Temperature (T9)',fontsize=14)
plt.legend()
savePath = os.path.join(desiredDir,'DataRatioSumReactionRates.png')
plt.savefig(savePath,dpi=300)
