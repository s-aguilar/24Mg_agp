import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

print('Attempting to create required directories: ')
try:
    os.mkdir('effPlots')
    os.mkdir('effPlots/comparison')
    os.mkdir('effPlots/detEffCurve')
    os.mkdir('effPlots/detEffCurve/log')
    os.mkdir('effPlots/detEffCurve/normal')
    os.mkdir('effPlots/detEffCurve/corresponding')
    os.mkdir('effPlots/detEffCurve/corresponding/log')
    os.mkdir('effPlots/detEffCurve/corresponding/normal')
except OSError:
    print ("Directories already exist!")
else:
    print('DONE!')



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
df_60Co_1173 = pd.read_csv('csv/60Co_1173cal.csv')
df_60Co_1332 = pd.read_csv('csv/60Co_1332cal.csv')
df_137Cs_661 = pd.read_csv('csv/137Cs_661cal.csv')


# Collect Runtime info, 60Co peaks calculated from same run
runtime_60Co = df_60Co_1173['Runtime']
runtime_137Cs = df_137Cs_661['Runtime']


# Calculate total number of decay events for the run
totalDecayEvents_60Co = np.array(runtime_60Co*currentActivity_60Co)
totalDecayEvents_137Cs = np.array(runtime_137Cs*currentActivity_137Cs)



# """
# Plot change in peak areas using different peak area methods
print('\nPlotting change in efficiency using different peak area methods')
start_time = time.time()

angle = np.array([120,105,90,45,30,15,0,-15,-30,-45,-90,-105,-120])
percentChange1173 = np.array(df_60Co_1173['percentChange'])
percentChange1332 = np.array(df_60Co_1332['percentChange'])
percentChange661 = np.array(df_137Cs_661['percentChange'])

plt.plot(angle,percentChange1173)
plt.xlabel('Angle (deg)')
plt.ylabel('Percent Change')
plt.title("60Co - 1173 keV - Percent Change")
plt.tight_layout()
plt.savefig('effPlots/comparison/percentchange1173.png')

plt.clf()
plt.plot(angle,percentChange1332)
plt.xlabel('Angle (deg)')
plt.ylabel('Percent Change')
plt.title("60Co - 1332 keV - Percent Change")
plt.tight_layout()
plt.savefig('effPlots/comparison/percentchange1332.png')

plt.clf()
plt.plot(angle,percentChange661)
plt.xlabel('Angle (deg)')
plt.ylabel('Percent Change')
plt.title("137Cs - 661 keV - Percent Change")
plt.tight_layout()
plt.savefig('effPlots/comparison/percentchange661.png')

end_time = time.time()
print('DONE!\t Process required: %f seconds'%(end_time - start_time))
# """



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

# """
# Plot efficiencies as function of angle for each calibration run
print('\nPlotting efficiencies as function of angle for each calibration run')
start_time = time.time()

plt.clf()
plt.errorbar(angle,eff_1173,yerr=effErr_1173,fmt='.')
plt.xlabel('Angle (deg)')
plt.xlim(-125,125)
plt.ylabel('Efficiency')
plt.ylim(0,.0006)
plt.title('$^{60}$Co - E$_{\\gamma}$ = 1173 keV - Efficiencies')
# plt.tight_layout()
plt.savefig('effPlots/1173_angular_eff.png',dpi=1200)

plt.clf()
plt.errorbar(angle,eff_1332,yerr=effErr_1332,fmt='.')
plt.xlabel('Angle (deg)')
plt.xlim(-125,125)
plt.ylabel('Efficiency')
plt.ylim(0,.0006)
plt.title('$^{60}$Co - E$_{\\gamma}$ = 1332 keV - Efficiencies')
# plt.tight_layout()
plt.savefig('effPlots/1332_angular_eff.png',dpi=1200)

plt.clf()
plt.errorbar(angle,eff_661,yerr=effErr_661,fmt='.')
plt.xlabel('Angle (deg)')
plt.xlim(-125,125)
plt.ylabel('Efficiency')
plt.ylim(0,.0006)
plt.title('$^{137}$Cs - E$_{\\gamma}$ = 661 keV - Efficiencies')
# plt.tight_layout()
plt.savefig('effPlots/661_angular_eff.png',dpi=1200)

end_time = time.time()
print('DONE!\t Process required: %f seconds'%(end_time - start_time))
# """


# Plot the efficiency curve for each detector
print('\nPlotting the efficiency curve for each detector')
start_time = time.time()

det_eff = np.vstack((eff_661,eff_1173,eff_1332))
det_effErr = np.vstack((effErr_661,effErr_1173,effErr_1332))

split_eff = np.hsplit(det_eff,13)
split_effErr = np.hsplit(det_effErr,13)

detName = np.array(df_60Co_1173['Detector'])
e_gam = np.array([661.657,1173.228,1332.492])
labels = ["$10^{-4}$","$10^{-3}$","$10^{-2}$"]
labelsx = ["$10^{2}$","$10^{3}$","$10^{4}$"]
fit_x = np.array([e_gam[0],e_gam[1],e_gam[2]])

# print(e_gam,split_eff[0])
# """
for xx in range(13):

    # Normal plot
    plt.clf()
    z = interp1d(e_gam,split_eff[xx].flatten(),fill_value="extrapolate")
    fit_y = z(fit_x)
    plt.plot(fit_x,fit_y,color='r',alpha=.75)
    plt.errorbar(e_gam,split_eff[xx],yerr=split_effErr[xx],fmt='.')
    plt.xlim(600,2000)
    plt.xlabel('Energy (KeV)')
    plt.ylabel('Efficiency')
    plt.grid(b=True,which='both',axis='y',alpha=.5)
    # plt.tight_layout()
    plt.title('%s Efficiency Curve'%detName[xx])
    plt.savefig('effPlots/detEffCurve/normal/%s.png'%detName[xx],dpi=1200)

    # Log plot
    plt.clf()
    plt.yscale('log')
    # plt.xscale('log')
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
    plt.savefig('effPlots/detEffCurve/log/log_%s.png'%detName[xx],dpi=1200)

end_time = time.time()
print('DONE!\t Process required: %f seconds'%(end_time - start_time))
# """


# """
# Plot corresponding detector efficiencies
print('\nPlotting matching detector efficiencies')
start_time = time.time()

ii = int(0)
while ii <= 5:

    # Normal plot
    plt.clf()
    plt.errorbar(e_gam,split_eff[ii],yerr=split_effErr[ii],fmt='.',color='b',label='Beam Right')
    plt.errorbar(e_gam,split_eff[12-ii],yerr=split_effErr[12-ii],fmt='.',color='r',label='Beam Left')
    plt.xlabel('Energy (KeV)')
    plt.ylabel('Efficiency')
    plt.grid(b=True,which='both',axis='y',alpha=.5)
    plt.title('Efficiency Curve - %ideg'%angle[ii])
    plt.legend()
    # plt.tight_layout()
    plt.savefig('effPlots/detEffCurve/corresponding/normal/%ieff.png'%angle[ii],dpi=1200)

    # Log plot
    plt.clf()
    plt.yscale('log')
    # plt.xscale('log')
    plt.errorbar(e_gam,split_eff[ii],yerr=split_effErr[ii],fmt='.',color='b',label='Beam Right')
    plt.errorbar(e_gam,split_eff[12-ii],yerr=split_effErr[12-ii],fmt='.',color='r',label='Beam Left')
    plt.xlabel('Energy (KeV)')
    plt.yticks((.0001,.001,.01),labels)
    plt.ylim(.0001,.001)
    plt.ylabel('Efficiency')
    plt.grid(b=True,which='both',axis='y',alpha=.5)
    plt.title('Efficiency Curve - %ideg'%angle[ii])
    plt.legend()
    # plt.tight_layout()
    plt.savefig('effPlots/detEffCurve/corresponding/log/log_%ieff.png'%angle[ii],dpi=1200)
    ii+=1

end_time = time.time()
print('DONE!\t Process required: %f seconds'%(end_time - start_time))
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
e_p2 = 1014.52
e_a1 = 1368.626

p1 = []
p2 = []
a1 = []

# linear interpolation to relate energy to channel
for yy in range(13):
    det_interp = interp1d(e_gam,split_cent[yy].T,fill_value="extrapolate")
    p1.append(det_interp(e_p1)[0])
    p2.append(det_interp(e_p2)[0])
    a1.append(det_interp(e_a1)[0])


d = {'Det':detName,'Angle':angle,'p1':p1,'p2':p2,'a1':a1}
df = pd.DataFrame(data=d)


# Save new DataFrame to a csv file
df.to_csv('csv/detectorCalibration.csv',sep=',',index=False)


########################################################################

# eff_peak = []
fit_x = np.array([e_p1,e_p2,e_a1])
p1eff = []
p2eff = []
a1eff = []
for zz in range(13):
    z = interp1d(e_gam,split_eff[zz].flatten(),fill_value="extrapolate")
    fit_y = z(fit_x)
    p1eff.append(fit_y[0])
    p2eff.append(fit_y[1])
    a1eff.append(fit_y[2])



# for x in range(13):
#     p1eff.append(eff_peak[x][0])
#     p2eff.append(eff_peak[x][1])
#     a1eff.append(eff_peak[x][2])

detName = list(detName)*229
angle = list(angle)*229
dd = {'Det':detName,'Angle':angle,'p1':p1eff*229,'p2':p2eff*229,'a1':a1eff*229}
dff = pd.DataFrame(data=dd)

# Save new DataFrame to a csv file
dff.to_csv('csv/detectorEfficiencies.csv',sep=',',index=False)
