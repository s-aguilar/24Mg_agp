import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# print('Attempting to create required directories: ')
# try:
#     os.mkdir('rMatrix')
#     os.mkdir('crossSection')
#     os.mkdir('crossSection/P1')
#     os.mkdir('crossSection/P2')
#     os.mkdir('crossSection/A1')
# except OSError:
#     print ("Directories already exist!")
# else:
#     print('DONE!')

detectors = ['det_h0-0','det_h0-1','det_h0-2','det_h0-3','det_h0-4','det_h0-5',
            'det_h0-6','det_h0-7','det_h1-0','det_h1-1','det_h1-2','det_h1-3',
            'det_h1-4']

# Read in the data into dataframe
dfeff = pd.read_csv('calibration/csv/detectorEfficiencies.csv')
df1 = pd.read_csv('Yields/P1/p1Yields.csv')
df2 = pd.read_csv('Yields/P2/p2Yields.csv')
df3 = pd.read_csv('Yields/A1/a1Yields.csv')

effp1 = dfeff['p1'].values
effp2 = dfeff['p2'].values
effa1 = dfeff['a1'].values
Angle = dfeff['Angle'].values

# Angle (deg) for each detector, negative is beam left, positive is beam right
#         00  01 02 03 04 05 06 07  10  11  12   13   14
angle = np.array([120,105,90,45,30,15,0,-15,-30,-45,-90,-105,-120])


thickness = 25/(1e6)                         # 25ug/cm^2
numOfTarget = thickness*(1/24)*6.022e23     # thickness * (mol/24 g) * N_a

q_e = 1.6e-19
scale = 1e-8    # 10^-8 C/pulse
q_corr = scale/(2*q_e)
barn_conv = 1/(1e-24)
solidAngle = 4*np.pi

# """
# Extract the columns of the DataFrame as numpy arrays


p1Run = df1['Run'].values
p1Det = df1['Detector'].values

p1Yield = df1['Yield'].values / q_corr
p1Yield_err = df1['Yield err'].values / q_corr

p1Yield_effcor = p1Yield / effp1
p1Yield_err_effcor = p1Yield_err / effp1

p1Cross = p1Yield_effcor / numOfTarget * barn_conv / solidAngle
p1Cross_err = p1Yield_err_effcor / numOfTarget* barn_conv / solidAngle

p1Fit = df1['Fit Status'].values
p1Ealpha = df1['Ea'].values/1000    # Convert keV to MeV


# Fit Status == 0 -> Good Fit
# Fit Status == 1 -> Bad Fit
#
# Mask for which the fit was bad

# """
# mask1Fit = (df1['Fit Status'] == 0) # Not currently in use
# mask1Fit = (df1['Area'] > 1000)
# for det in range(len(detectors)):
#
#     maskDet = ((df1['Detector']==detectors[det]) & mask1Fit)
#
#     plt.clf()
#     plt.errorbar(p1Ealpha[maskDet],p1Cross[maskDet],yerr=p1Cross_err[maskDet],fmt='b.',markersize='2')
#     plt.yscale('log')
#     plt.ylim(1e-6,1)
#     plt.xlim(4,5.6)
#     plt.xlabel('$E_{\\alpha}$ (MeV)')
#     plt.ylabel('Cross-Section (barns)')
#     plt.title('p1 %s    %d$^{\circ}$'%(detectors[det],angle[det]))
#     plt.savefig('yieldPlots/P1/p1_%s.png'%detectors[det],dpi=300)

# with open("rMatrix_p1.dat","w") as f:
#     for loop in range(len(p1Cross)):
#         printOut= '%f \t %d \t %.8f \t %.8f \n' %(p1Ealpha[loop],Angle[loop],p1Cross[loop],p1Cross_err[loop])
#         f.write(printOut)

AnglesList=['0','15','30','45','90','105','120']

# test2 = []
f = open("rMatrix/rMatrix_p1.dat","w")
f.close()
for ang in AnglesList:
    _p1Ealpha = []
    _Angle = []
    _p1Cross = []
    _p1Cross_err = []

    # Average out over same angle
    for x in range(229):    # total of 229 runs
        _p1Ealpha.append(p1Ealpha[int(13*x)])
        if ang == '0':
            # print(p1Cross[x*13+6])
            _Angle.append(Angle[x*13+6])
            _p1Cross.append(p1Cross[x*13+6])
            # test2.append(p1Cross[x*13+6])
            err = p1Cross_err[x*13+6]
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5  )   # inflating the error bars for 5% systematic uncertainty
        elif ang == '15':
            _Angle.append( (abs(Angle[x*13+5])+abs(Angle[x*13+7]))/2 )
            _p1Cross.append( (p1Cross[x*13+5]+p1Cross[x*13+7])/2 )
            # _p1Cross_err.append( (p1Cross_err[x*13+5]+p1Cross_err[x*13+7])/2 )
            err = (p1Cross_err[x*13+5]**2+p1Cross_err[x*13+7]**2)**.5
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '30':
            _Angle.append( int( (abs(Angle[x*13+4])+abs(Angle[x*13+8]))/2 ) )
            _p1Cross.append( (p1Cross[x*13+4]+p1Cross[x*13+8])/2 )
            err = (p1Cross_err[x*13+4]**2+p1Cross_err[x*13+8]**2)**.5
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '45':
            _Angle.append( int( (abs(Angle[x*13+3])+abs(Angle[x*13+9]))/2 ) )
            _p1Cross.append( (p1Cross[x*13+3]+p1Cross[x*13+9])/2 )
            err = (p1Cross_err[x*13+3]**2+p1Cross_err[x*13+9]**2)**.5
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '90':
            _Angle.append( int( (abs(Angle[x*13+2])+abs(Angle[x*13+10]))/2 ) )
            _p1Cross.append( (p1Cross[x*13+2]+p1Cross[x*13+10])/2 )
            err = (p1Cross_err[x*13+2]**2+p1Cross_err[x*13+10]**2)**.5
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '105':
            # _Angle.append( int( (abs(Angle[x*13+1])+abs(Angle[x*13+11]))/2 ) )
            _Angle.append( int( 75 ) )
            _p1Cross.append( (p1Cross[x*13+1]+p1Cross[x*13+11])/2 )
            err = (p1Cross_err[x*13+1]**2+p1Cross_err[x*13+11]**2)**.5
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '120':
            # _Angle.append( int( (abs(Angle[x*13])+abs(Angle[x*13+12]))/2 ) )
            _Angle.append( int( 60 ) )
            _p1Cross.append( (p1Cross[x*13]+p1Cross[x*13+12])/2 )
            err = (p1Cross_err[x*13]**2+p1Cross_err[x*13+12]**2)**.5
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )


    _p1Ealpha = np.array(_p1Ealpha)
    _Angle = np.array(_Angle)
    _p1Cross = np.array(_p1Cross)
    _p1Cross_err = np.array(_p1Cross_err)

    # Sort by energy, keeping others consistent!
    ind = _p1Ealpha.argsort()
    _p1Ealpha = _p1Ealpha[ind]
    _Angle = _Angle[ind]
    _p1Cross = _p1Cross[ind]
    _p1Cross_err = _p1Cross_err[ind]


    # Make the Cross-Section plot
    plt.clf()
    plt.plot(_p1Ealpha,_p1Cross)
    plt.errorbar(_p1Ealpha,_p1Cross,yerr=_p1Cross_err,fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(4,5.6)
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Differential Cross-Section (barns/sr)')
    plt.title('$^{24}$Mg($\\mathrm{\\alpha,\\gamma p_{1}}$)$^{27}$Al\t%s$^{\circ}$'%ang)
    plt.savefig('crossSection/P1/p1_%s.png'%ang,dpi=300)
    plt.clf()


    with open("rMatrix/rMatrix_p1.dat","a") as f:
        for loop in range(229):
            printOut= '%f \t %d \t %.8f \t %.8f \n' %(_p1Ealpha[loop],_Angle[loop],_p1Cross[loop],_p1Cross_err[loop])
            f.write(printOut)

# test1 = set(p1Cross[(df1['Detector']=='det_h0-6')][:16])
# print(test1)
# print('\n',test2[:16])
# test2 = set(test2[:16])
# if (test1==test2): print("its the same!")
# else: print("fuck its different!")

# """
# Extract the columns of the DataFrame as numpy arrays
p2Run = df2['Run'].values
p2Det = df2['Detector'].values

p2Yield = df2['Yield'].values / q_corr
p2Yield_err = df2['Yield err'].values / q_corr

p2Yield_effcor = p2Yield / effp2
p2Yield_err_effcor = p2Yield_err / effp2

p2Cross = p2Yield_effcor / numOfTarget * barn_conv / solidAngle
p2Cross_err = p2Yield_err_effcor / numOfTarget* barn_conv / solidAngle

p2Fit = df2['Fit Status'].values
p2Ealpha = df2['Ea'].values/1000    # Convert keV to MeV

# Fit Status == 0 -> Good Fit
# Fit Status == 1 -> Bad Fit
#
# Mask for which the fit was bad
# mask2Fit = (df2['Fit Status'] == 0)
# mask2Fit = (df2['Area'] > 1000)
#
# for det in range(len(detectors)):
#
#     maskDet = ((df2['Detector']==detectors[det]) & mask1Fit)
#
#     plt.clf()
#     plt.errorbar(p2Ealpha[maskDet],p2Cross[maskDet],yerr=p2Cross_err[maskDet],fmt='b.',markersize='2')
#     plt.yscale('log')
#     plt.ylim(1e-6,1)
#     plt.xlim(4,5.6)
#     plt.xlabel('$E_{\\alpha}$ (MeV)')
#     plt.ylabel('Cross-Section (barns/sr)')
#     plt.title('p2 %s    %d$^{\circ}$'%(detectors[det],angle[det]))
#     plt.savefig('yieldPlots/P2/p2_%s.png'%detectors[det],dpi=300)


f = open("rMatrix/rMatrix_p2.dat","w")
f.close()
for ang in AnglesList:
    _p2Ealpha = []
    _Angle = []
    _p2Cross = []
    _p2Cross_err = []

    # Average out over same angle
    for x in range(229):    # total of 229 runs
        _p2Ealpha.append(p2Ealpha[int(13*x)])
        if ang == '0':
            # print(p2Cross[x*13+6])
            _Angle.append(Angle[x*13+6])
            _p2Cross.append(p2Cross[x*13+6])
            err = p2Cross_err[x*13+6]
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '15':
            _Angle.append( (abs(Angle[x*13+5])+abs(Angle[x*13+7]))/2 )
            _p2Cross.append( (p2Cross[x*13+5]+p2Cross[x*13+7])/2 )
            err = (p2Cross_err[x*13+5]**2+p2Cross_err[x*13+7]**2)**.5
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '30':
            _Angle.append( int( (abs(Angle[x*13+4])+abs(Angle[x*13+8]))/2 ) )
            _p2Cross.append( (p2Cross[x*13+4]+p2Cross[x*13+8])/2 )
            err = (p2Cross_err[x*13+4]**2+p2Cross_err[x*13+8]**2)**.5
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '45':
            _Angle.append( int( (abs(Angle[x*13+3])+abs(Angle[x*13+9]))/2 ) )
            _p2Cross.append( (p2Cross[x*13+3]+p2Cross[x*13+9])/2 )
            err = (p2Cross_err[x*13+3]**2+p2Cross_err[x*13+9]**2)**.5
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '90':
            _Angle.append( int( (abs(Angle[x*13+2])+abs(Angle[x*13+10]))/2 ) )
            _p2Cross.append( (p2Cross[x*13+2]+p2Cross[x*13+10])/2 )
            err = (p2Cross_err[x*13+2]**2+p2Cross_err[x*13+10]**2)**.5
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '105':
            # _Angle.append( int( (abs(Angle[x*13+1])+abs(Angle[x*13+11]))/2 ) )
            _Angle.append( int( 75 ) )
            _p2Cross.append( (p2Cross[x*13+1]+p2Cross[x*13+11])/2 )
            err = (p2Cross_err[x*13+1]**2+p2Cross_err[x*13+11]**2)**.5
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '120':
            # _Angle.append( int( (abs(Angle[x*13])+abs(Angle[x*13+12]))/2 ) )
            _Angle.append( int( 60 ) )
            _p2Cross.append( (p2Cross[x*13]+p2Cross[x*13+12])/2 )
            err = (p2Cross_err[x*13]**2+p2Cross_err[x*13+12]**2)**.5
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )


    _p2Ealpha = np.array(_p2Ealpha)
    _Angle = np.array(_Angle)
    _p2Cross = np.array(_p2Cross)
    _p2Cross_err = np.array(_p2Cross_err)

    # Sort by energy, keeping others consistent!
    ind = _p2Ealpha.argsort()
    _p2Ealpha = _p2Ealpha[ind]
    _Angle = _Angle[ind]
    _p2Cross = _p2Cross[ind]
    _p2Cross_err = _p2Cross_err[ind]


    # Make the Cross-Section plot
    plt.clf()
    plt.plot(_p2Ealpha,_p2Cross)
    plt.errorbar(_p2Ealpha,_p2Cross,yerr=_p2Cross_err,fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(4,5.6)
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Differential Cross-Section (barns/sr)')
    plt.title('$^{24}$Mg($\\mathrm{\\alpha,\\gamma p_{2}}$)$^{27}$Al\t%s$^{\circ}$'%ang)
    plt.savefig('crossSection/P2/p2_%s.png'%ang,dpi=300)
    plt.clf()

    with open("rMatrix/rMatrix_p2.dat","a") as f:
        for loop in range(229):
            printOut= '%f \t %d \t %.8f \t %.8f \n' %(_p2Ealpha[loop],_Angle[loop],_p2Cross[loop],_p2Cross_err[loop])
            f.write(printOut)
# """


# """
# Extract the columns of the DataFrame as numpy arrays
a1Run = df3['Run'].values
a1Det = df3['Detector'].values

a1Yield = df3['Yield'].values/ q_corr
a1Yield_err = df3['Yield err'].values/ q_corr

a1Yield_effcor = a1Yield / effa1
a1Yield_err_effcor = a1Yield_err / effa1

a1Cross = a1Yield_effcor / numOfTarget * barn_conv / solidAngle
a1Cross_err = a1Yield_err_effcor / numOfTarget* barn_conv / solidAngle

a1Fit = df3['Fit Status'].values
a1Ealpha = df3['Ea'].values/1000    # Convert keV to MeV

# Fit Status == 0 -> Good Fit
# Fit Status == 1 -> Bad Fit
#
# Mask for which the fit was bad
# mask3Fit = (df3['Fit Status'] == 0)
# mask3Fit = (df3['Area'] > 1000)

# for det in range(len(detectors)):
#
#     maskDet = ((df3['Detector']==detectors[det]) & mask3Fit)
#
#     plt.clf()
#     plt.errorbar(a1Ealpha[maskDet],a1Cross[maskDet],yerr=a1Cross_err[maskDet],fmt='b.',markersize='2')
#     plt.yscale('log')
#     plt.ylim(1e-6,1)
#     plt.xlim(4,5.6)
#     plt.xlabel('$E_{\\alpha}$ (MeV)')
#     plt.ylabel('Cross-Section (barns/sr)')
#     plt.title('a1 %s    %d$^{\circ}$'%(detectors[det],angle[det]))
#     plt.savefig('yieldPlots/A1/a1_%s.png'%detectors[det],dpi=300)

f = open("rMatrix/rMatrix_a1.dat","w")
f.close()

for ang in AnglesList:
    _a1Ealpha = []
    _Angle = []
    _a1Cross = []
    _a1Cross_err = []

    # Average out over same angle
    for x in range(229):    # total of 229 runs
        _a1Ealpha.append(a1Ealpha[int(13*x)])
        if ang == '0':
            # print(a1Cross[x*13+6])
            _Angle.append(Angle[x*13+6])
            _a1Cross.append(a1Cross[x*13+6])
            err = a1Cross_err[x*13+6]
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 ) # inflate errorbar
        elif ang == '15':
            _Angle.append( (abs(Angle[x*13+5])+abs(Angle[x*13+7]))/2 )
            _a1Cross.append( (a1Cross[x*13+5]+a1Cross[x*13+7])/2 )
            err = (a1Cross_err[x*13+5]**2+a1Cross_err[x*13+7]**2)**.5
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '30':
            _Angle.append( int( (abs(Angle[x*13+4])+abs(Angle[x*13+8]))/2 ) )
            _a1Cross.append( (a1Cross[x*13+4]+a1Cross[x*13+8])/2 )
            err = (a1Cross_err[x*13+4]**2+a1Cross_err[x*13+8]**2)**.5
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '45':
            _Angle.append( int( (abs(Angle[x*13+3])+abs(Angle[x*13+9]))/2 ) )
            _a1Cross.append( (a1Cross[x*13+3]+a1Cross[x*13+9])/2 )
            err = (a1Cross_err[x*13+3]**2+a1Cross_err[x*13+9]**2)**.5
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '90':
            _Angle.append( int( (abs(Angle[x*13+2])+abs(Angle[x*13+10])) )/2 )
            _a1Cross.append( (a1Cross[x*13+2]+a1Cross[x*13+10])/2 )
            err = (a1Cross_err[x*13+2]**2+a1Cross_err[x*13+10]**2)**.5
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '105':
            # _Angle.append( int( (abs(Angle[x*13+1])+abs(Angle[x*13+11]))/2 ) )
            _Angle.append( int( 75 ) )
            _a1Cross.append( (a1Cross[x*13+1]+a1Cross[x*13+11])/2 )
            err = (a1Cross_err[x*13+1]**2+a1Cross_err[x*13+11]**2)**.5
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '120':
            # _Angle.append( int( (abs(Angle[x*13])+abs(Angle[x*13+12]))/2 ) )
            _Angle.append( int( 60 ) )
            _a1Cross.append( (a1Cross[x*13]+a1Cross[x*13+12])/2 )
            err = (a1Cross_err[x*13]**2+a1Cross_err[x*13+12]**2)**.5
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )


    _a1Ealpha = np.array(_a1Ealpha)
    _Angle = np.array(_Angle)
    _a1Cross = np.array(_a1Cross)
    _a1Cross_err = np.array(_a1Cross_err)

    # Sort by energy, keeping others consistent!
    ind = _a1Ealpha.argsort()
    _a1Ealpha = _a1Ealpha[ind]
    _Angle = _Angle[ind]
    _a1Cross = _a1Cross[ind]
    _a1Cross_err = _a1Cross_err[ind]

    # Make the Cross-Section plot
    plt.clf()
    plt.plot(_a1Ealpha,_a1Cross)
    plt.errorbar(_a1Ealpha,_a1Cross,yerr=_a1Cross_err,fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(4,5.6)
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Differential Cross-Section (barns/sr)')
    plt.title('$^{24}$Mg($\\mathrm{\\alpha,\\gamma\\alpha_{1}}$)$^{24}$Mg\t%s$^{\circ}$'%ang)
    plt.savefig('crossSection/A1/a1_%s.png'%ang,dpi=300)
    plt.clf()

    with open("rMatrix/rMatrix_a1.dat","a") as f:
        for loop in range(229):
            printOut= '%f \t %d \t %.8f \t %.8f \n' %(_a1Ealpha[loop],_Angle[loop],_a1Cross[loop],_a1Cross_err[loop])
            f.write(printOut)
# """

print('DONE!')
