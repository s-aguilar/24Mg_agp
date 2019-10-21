import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# print('Attempting to create required directories: ')
# try:
#     os.mkdir('yieldPlots')
#     os.mkdir('yieldPlots/P1')
#     os.mkdir('yieldPlots/P2')
#     os.mkdir('yieldPlots/A1')
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
df1 = df1.sort_values(by=['Ea'])
df2 = pd.read_csv('Yields/P2/p2Yields.csv')
df2 = df2.sort_values(by=['Ea'])
df3 = pd.read_csv('Yields/A1/a1Yields.csv')
df3 = df3.sort_values(by=['Ea'])



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

# """
# Extract the columns of the DataFrame as numpy arrays


p1Run = df1['Run'].values
p1Det = df1['Detector'].values

p1Yield = df1['Yield'].values / q_corr
p1Yield_err = df1['Yield err'].values / q_corr

p1Yield_effcor = p1Yield / effp1
p1Yield_err_effcor = p1Yield_err / effp1

p1Cross = p1Yield_effcor/numOfTarget * barn_conv
p1Cross_err = p1Yield_err_effcor/numOfTarget* barn_conv

p1Fit = df1['Fit Status'].values
p1Ealpha = df1['Ea'].values/1000    # Convert keV to MeV


# Fit Status == 0 -> Good Fit
# Fit Status == 1 -> Bad Fit
#
# Mask for which the fit was bad

# mask1Fit = (df1['Fit Status'] == 0)
#
# for det in range(len(detectors)):
#
#     maskDet = ((df1['Detector']==detectors[det])) # & mask1Fit & )
#
#     plt.clf()
#     print (detectors[det])  #,len(p1Yield[maskDet])/len(p1Yield[(df1['Detector']==detectors[det])])*100)
#     # plt.scatter(p1Ealpha[maskDet],p1Yield[maskDet],c='b',marker='.')
#     plt.errorbar(p1Ealpha[maskDet],p1Cross[maskDet],yerr=p1Cross_err[maskDet],fmt='b.',markersize='1')
#     plt.yscale('log')
#     plt.ylim(1e-6,1)
#     plt.xlim(4,5.6)
#     plt.xlabel('$E_{\\alpha}$ (MeV)')
#     plt.ylabel('Cross-Section (barns)')
#     plt.title('p1 %s    %d$^{\circ}$'%(detectors[det],angle[det]))
#     plt.savefig('yieldPlots/P1/p1_%s.png'%detectors[det],dpi=300)
#     # plt.show()
# """


AnglesList=['0','15','30','45','90','105','120']


for ang in AnglesList:
    _p1Ealpha = []
    _Angle = []
    _p1Cross = []
    _p1Cross_err = []

    # Average out over same angle
    for x in range(229):    # total of 229 runs
        _p1Ealpha.append(p1Ealpha[int(13*x)])
        if ang == '0':
            _Angle.append(Angle[int((x*13)+6)])
            _p1Cross.append(p1Cross[int((x*13)+6)])
            _p1Cross_err.append(p1Cross_err[int((x*13)+6)])
        elif ang == '15':
            _Angle.append( int( (abs(Angle[int((x*13)+5)])+abs(Angle[int((x*13)+7)]))/2 ) )
            _p1Cross.append( (p1Cross[int((x*13)+5)]+p1Cross[int((x*13)+7)])/2 )
            _p1Cross_err.append( (p1Cross_err[int((x*13)+5)]+p1Cross_err[int((x*13)+7)])/2 )
        elif ang == '30':
            _Angle.append( int( (abs(Angle[int((x*13)+4)])+abs(Angle[int((x*13)+8)]))/2 ) )
            _p1Cross.append( (p1Cross[int((x*13)+4)]+p1Cross[int((x*13)+8)])/2 )
            _p1Cross_err.append( (p1Cross_err[int((x*13)+4)]+p1Cross_err[int((x*13)+8)])/2 )
        elif ang == '45':
            _Angle.append( int( (abs(Angle[int((x*13)+3)])+abs(Angle[int((x*13)+9)]))/2 ) )
            _p1Cross.append( (p1Cross[int((x*13)+3)]+p1Cross[int((x*13)+9)])/2 )
            _p1Cross_err.append( (p1Cross_err[int((x*13)+3)]+p1Cross_err[int((x*13)+9)])/2 )
        elif ang == '90':
            _Angle.append( int( (abs(Angle[int((x*13)+2)])+abs(Angle[int((x*13)+10)]))/2 ) )
            _p1Cross.append( (p1Cross[int((x*13)+2)]+p1Cross[int((x*13)+10)])/2 )
            _p1Cross_err.append( (p1Cross_err[int((x*13)+2)]+p1Cross_err[int((x*13)+10)])/2 )
        elif ang == '105':
            _Angle.append( int( (abs(Angle[int((x*13)+1)])+abs(Angle[int((x*13)+11)]))/2 ) )
            _p1Cross.append( (p1Cross[int((x*13)+1)]+p1Cross[int((x*13)+11)])/2 )
            _p1Cross_err.append( (p1Cross_err[int((x*13)+1)]+p1Cross_err[int((x*13)+11)])/2 )
        elif ang == '120':
            _Angle.append( int( (abs(Angle[int(x*13)])+abs(Angle[int(x*13+12)]))/2 ) )
            _p1Cross.append( (p1Cross[int((x*13))]+p1Cross[int((x*13)+7)])/2 )
            _p1Cross_err.append( (p1Cross_err[int((x*13))]+p1Cross_err[int((x*13)+12)])/2 )

    plt.errorbar(_p1Ealpha,_p1Cross,yerr=_p1Cross_err,fmt='b.',markersize='1')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(4,5.6)
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Cross-Section (barns)')
    plt.title('p1 %s$^{\circ}$'%ang)
    plt.savefig('yieldPlots/P1/p1_%s.png'%ang,dpi=300)
    plt.clf()

    with open("rMatrix/rMatrix_p1_%s.dat"%ang,"w") as f:
        for loop in range(229):
            printOut= '%f \t %d \t %.8f \t %.8f \n' %(_p1Ealpha[loop],_Angle[loop],_p1Cross[loop],_p1Cross_err[loop])
            f.write(printOut)


# with open("rMatrix_p1.dat","w") as f:
#     for loop in range(len(p1Cross)):
#         printOut= '%f \t %d \t %.8f \t %.8f \n' %(p1Ealpha[loop],Angle[loop],p1Cross[loop],p1Cross_err[loop])
#         f.write(printOut)



"""
# Extract the columns of the DataFrame as numpy arrays
p2Run = df2['Run'].values
p2Det = df2['Detector'].values
p2Yield = df2['Yield'].values
p2Yield_err = df2['Yield err'].values
p2Fit = df2['Fit Status'].values
p2Ealpha = df2['Ea'].values/1000    # Convert keV to MeV

# Fit Status == 0 -> Good Fit
# Fit Status == 1 -> Bad Fit
#
# Mask for which the fit was bad
mask2Fit = (df2['Fit Status'] == 0)

for det in detectors:

    maskDet = ( (df2['Detector']==det) & mask2Fit )

    plt.clf()
    print (det,len(p2Yield[maskDet])/len(p2Yield[(df2['Detector']==det)])*100)
    plt.plot(p2Ealpha[maskDet],p2Yield[maskDet],c='b',marker='.')
    plt.errorbar(p2Ealpha[maskDet],p2Yield[maskDet],yerr=p2Yield_err[maskDet],fmt='b.')
    # plt.yscale('log')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Yield')
    plt.title('p2 %s'%det)
    plt.savefig('yieldPlots/P2/p2_%s.png'%det,dpi=900)
    # plt.show()
# """

"""
# Extract the columns of the DataFrame as numpy arrays
a1Run = df3['Run'].values
a1Det = df3['Detector'].values
a1Yield = df3['Yield'].values
a1Yield_err = df3['Yield err'].values
a1Fit = df3['Fit Status'].values
a1Ealpha = df3['Ea'].values/1000    # Convert keV to MeV

# Fit Status == 0 -> Good Fit
# Fit Status == 1 -> Bad Fit
#
# Mask for which the fit was bad
mask3Fit = (df3['Fit Status'] == 0)

for det in detectors:

    maskDet = ( (df3['Detector']==det) & mask3Fit )

    plt.clf()
    print (det,len(a1Yield[maskDet])/len(a1Yield[(df3['Detector']==det)])*100)
    plt.plot(a1Ealpha[maskDet],a1Yield[maskDet],c='b',marker='.')
    plt.errorbar(a1Ealpha[maskDet],a1Yield[maskDet],yerr=a1Yield_err[maskDet],fmt='b.')
    # plt.yscale('log')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Yield')
    plt.title('a1 %s'%det)
    plt.savefig('yieldPlots/A1/a1_%s.png'%det,dpi=900)
    # plt.show()
# """


print('DONE!')
