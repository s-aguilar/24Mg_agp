import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def newActivity(act,t,half):
    lam = np.log(2)/half
    N_0 = act/lam
    N_t = N_0*np.exp(-lam*t)
    return lam*N_t


"""############################ MAIN #######################################"""
# Old source activity and date measured
activity_60Co = 34.51e3     # kBq
date_60Co = "08/01/2016"

activity_137Cs = 3.752e3
date_137Cs = "08/01/2007"


# Source half-lives
half_60Co = 5.2747          # years
half_60Co *= 365*24*60*60   # seconds

half_137Cs = 30.08
half_137Cs *= 365*24*60*60


# Calculate source activity  on date of experiment
date_exp = "09/15/2018"

span_60Co = len( pd.date_range(date_60Co,date_exp) )    # of days
span_60Co *= 24*60*60                                   # of seconds

span_137Cs = len( pd.date_range(date_137Cs,date_exp) )
span_137Cs *= 24*60*60

currentActivity_60Co = newActivity(activity_60Co,span_60Co,half_60Co)
currentActivity_137Cs = newActivity(activity_137Cs,span_137Cs,half_137Cs)


# Read in calibration information in pandas DataFrame
df_60Co_1173 = pd.read_csv('60Co_1173cal.csv')
df_60Co_1332 = pd.read_csv('60Co_1332cal.csv')
df_137Cs_661 = pd.read_csv('137Cs_661cal.csv')

# print (df_60Co_1173.head())


# Collect Runtime info, 60Co peaks calculated from same run
runtime_60Co = df_60Co_1173['Runtime']
runtime_137Cs = df_137Cs_661['Runtime']


# Calculate total number of decay events for the run
totalDecayEvents_60Co = np.array(runtime_60Co*currentActivity_60Co)
totalDecayEvents_137Cs = np.array(runtime_137Cs*currentActivity_137Cs)



# Gather peak areas and their errors
area_1173 = np.array(df_60Co_1173['Area'])
areaErr_1173 = np.array(df_60Co_1173['Area err'])

area_1332 = np.array(df_60Co_1332['Area'])
areaErr_1332 = np.array(df_60Co_1332['Area err'])

area_661 = np.array(df_137Cs_661['Area'])
areaErr_661 = np.array(df_137Cs_661['Area err'])


# Calculate Efficiencies and Efficiency Err
eff_1173 = area_1173/totalDecayEvents_60Co
effErr_1173 = areaErr_1173/totalDecayEvents_60Co

eff_1332 = area_1332/totalDecayEvents_60Co
effErr_1332 = areaErr_1332/totalDecayEvents_60Co

eff_661 = area_661/totalDecayEvents_137Cs
effErr_661 = areaErr_661/totalDecayEvents_137Cs

avgEff = (eff_661+eff_1332+eff_1173)/3
avgEffErr = (effErr_661+effErr_1332+effErr_1173)/3


# Angle (deg) for each detector, negative is beam left, positive is beam right
#         00  01 02 03 04 05 06 07  10  11  12   13   14
angle = np.array([120,105,90,45,30,15,0,-15,-30,-45,-90,-105,-120])

"""
# Plot efficiencies as function of angle for each calibration run
plt.errorbar(angle,eff_1173,yerr=effErr_1173,fmt='.')
plt.xlabel('Angle (deg)')
plt.xlim(-125,125)
plt.ylabel('Efficiency')
plt.ylim(0,.0006)
plt.title('$^{60}$Co - E$_{\\gamma}$ = 1173 keV - Efficiencies')
# plt.tight_layout()
plt.savefig('effPlots/1173_angular_eff.png',dpi=1200)
# plt.show()

plt.clf()
plt.errorbar(angle,eff_1332,yerr=effErr_1332,fmt='.')
plt.xlabel('Angle (deg)')
plt.xlim(-125,125)
plt.ylabel('Efficiency')
plt.ylim(0,.0006)
plt.title('$^{60}$Co - E$_{\\gamma}$ = 1332 keV - Efficiencies')
# plt.tight_layout()
plt.savefig('effPlots/1332_angular_eff.png',dpi=1200)
# plt.show()

plt.clf()
plt.errorbar(angle,eff_661,yerr=effErr_661,fmt='.')
plt.xlabel('Angle (deg)')
plt.xlim(-125,125)
plt.ylabel('Efficiency')
plt.ylim(0,.0006)
plt.title('$^{137}$Cs - E$_{\\gamma}$ = 661 keV - Efficiencies')
# plt.tight_layout()
plt.savefig('effPlots/661_angular_eff.png',dpi=1200)
# plt.show()
#"""

# plot the efficiency curve for each detector
det_eff = np.vstack((eff_661,eff_1173,eff_1332))
det_effErr = np.vstack((effErr_661,effErr_1173,effErr_1332))

split_eff = np.hsplit(det_eff,13)
split_effErr = np.hsplit(det_effErr,13)


detName = np.array(df_60Co_1173['Detector'])
e_gam = np.array([661.657,1173.228,1332.492])
labels = ["$10^{-4}$","$10^{-3}$","$10^{-2}$"]
labelsx = ["$10^{2}$","$10^{3}$","$10^{4}$"]
fit_x = np.array([100,e_gam[0],e_gam[1],e_gam[2],2000])
# """
for xx in range(13):
    plt.clf()
    plt.yscale('log')
    plt.xscale('log')
    z = interp1d(e_gam,split_eff[xx].flatten(),fill_value="extrapolate")
    fit_y = z(fit_x)
    plt.plot(fit_x,fit_y,color='r',alpha=.75)
    plt.errorbar(e_gam,split_eff[xx],yerr=split_effErr[xx],fmt='.')
    plt.xlim(600,2000)
    # plt.xticks((100,1000,10000),labelsx)
    plt.xlabel('Energy (KeV)')
    plt.yticks((.0001,.001,.01),labels)
    plt.ylim(.0001,.001)
    plt.ylabel('Efficiency')
    plt.grid(b=True,which='both',axis='y',alpha=.5)
    # plt.tight_layout()
    plt.title('%s Efficiency Curve'%detName[xx])
    plt.savefig('effPlots/detEffCurve/%s.png'%detName[xx],dpi=1200)
    # plt.show()
    plt.clf()


# """

###############################################################################
# Estimate the p1, p2 and a1 channel location per detector via linear interpolation
###############################################################################

# Get the centroid location info
cent_1173 = np.array(df_60Co_1173['Centroid'])
cent_1332 = np.array(df_60Co_1332['Centroid'])
cent_661 = np.array(df_137Cs_661['Centroid'])

centroids = np.vstack((cent_661,cent_1173,cent_1332))
split_cent = np.hsplit(centroids,13)

e_p1 = 843.76   # keV
e_p2 = 1014.56
e_a1 = 1368.67

p1 = []
p2 = []
a1 = []

# linear interpolation to relate energy to channel
for yy in range(13):
    det_interp = interp1d(e_gam,split_cent[yy].T,fill_value="extrapolate")
    p1.append(det_interp(e_p1)[0])
    p2.append(det_interp(e_p2)[0])
    a1.append(det_interp(e_a1)[0])


d = {'Det':detName,'Angle':angle,'p1':p1,'p2':p2,'a1':a1,'eff':avgEff,'eff err':avgEffErr}
df = pd.DataFrame(data=d)
# print(df.head(13))


# Save new DataFrame to a csv file
df.to_csv('detectorEfficiencies.csv',sep=',',index=False)
