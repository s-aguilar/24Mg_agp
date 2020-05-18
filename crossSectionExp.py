import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
I normalized the cross-sections at certain angles with some factor James
calculated from an isotropic resonance. This is to correct for issues we had
with the efficiencies. In addition, I manually suppressed some large
error bars such that it does not disable log plotting in Azure. For this I
artificially assign an error bar. This feature can easily be disabled
# """

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

# This function checks to see if the errror is larger than the measurement,
# if it is it just artificially assigns it an error
def check(x,x_err):
    if x_err >= x:
        return x*.9
    #if x_err >= x: return x_err
    else: return x_err



    print('Attempting to create required parent directories: ')
    try:
        os.mkdir('rMatrix')
        os.mkdir('crossSection')
    except OSError:
        print ("Directories already exist!")
    else:
        print('DONE!')




# Angle (deg) for each detector, negative is beam left, positive is beam right
#                  00  01 02 03 04 05 06 07  10  11  12   13   14
angle = np.array([120,105,90,45,30,15,0,-15,-30,-45,-90,-105,-120])
detectors = ['det_h0-0','det_h0-1','det_h0-2','det_h0-3','det_h0-4','det_h0-5',
            'det_h0-6','det_h0-7','det_h1-0','det_h1-1','det_h1-2','det_h1-3',
            'det_h1-4']

mirrAng = np.array([120,-120,105,-105,90,-90,45,-45,30,-30,15,-15])
mirrDet = ['det_h0-0','det_h1-4','det_h0-1','det_h1-3','det_h0-2','det_h1-2',
           'det_h0-3','det_h1-1','det_h0-4','det_h1-0','det_h0-5','det_h0-7']
norm0 = 1
norm15 = 1
norm30 = 1
norm45 = 1
norm90 = 1
norm105 = .87
norm120 = .81


thickness = 34.87/(1e6)                        # 37ug/cm^2 converted to g/cm^2
numOfTarget = thickness*(1/23.985)*6.022e23     # thickness * (mol/23.985 g) * N_a
# numOfTarget=1 ###############################################################################################################################################################
q_e = 1.6022e-19
scale = 1e-8    # 10^-8 C/pulse
q_corr = scale/(2*q_e)
barn_conv = 1/(1e-24)
solidAngle = 4*np.pi


channels = ['P1','P2','A1']
# channels = ['P3']
# channels = ['O17_1','O17_ng']
for chan in channels:

    if (chan == 'O17_1' or chan == 'O17_ng'): numOfTarget = 1

    print('Attempting to create required directories for channels: ')
    try:
        os.mkdir('crossSection/%s'%chan)
    except OSError:
        print ("Directories already exist!")
    else:
        print('DONE!')

    print('Working on %s channel'%chan)

    # Read in the data into dataframe
    dfeff = pd.read_csv('calibration/csv/detectorEfficiencies.csv')
    df = pd.read_csv('Yields/%s/%sYields.csv'%(chan,chan.lower()))

    Angle = dfeff['Angle'].to_numpy()

    # Let the O17_1 (P1 buddy peak) be same efficiency as P1
    if (chan != 'O17_1' and chan != 'O17_ng' and chan != 'P3'):
        eff = dfeff['%s'%chan.lower()].to_numpy()
    elif (chan == 'O17_1'):
        eff = dfeff['p1'].to_numpy()
    elif (chan == 'O17_ng'):
        # eff = np.ones_like(dfeff['p1'].to_numpy())
        eff = dfeff['a1'].to_numpy()
    elif (chan == 'P3'):
        # eff = np.ones_like(dfeff['p1'].to_numpy())
        eff = dfeff['a1'].to_numpy()


    # Extract the columns of the DataFrame as numpy arrays
    chanRun = df['Run'].to_numpy()
    chanDet = df['Detector'].to_numpy()

    chanYield = df['Yield'].to_numpy() / q_corr
    chanYield_err = df['Yield err'].to_numpy() / q_corr

    chanYield_effcor = chanYield / eff / solidAngle
    chanYield_err_effcor = chanYield_err / eff / solidAngle

    chanCross = chanYield_effcor / numOfTarget * barn_conv
    chanCross_err = chanYield_err_effcor / numOfTarget * barn_conv
    chanFit = df['Fit Status'].to_numpy()
    chanEalpha = df['Ea'].to_numpy()/1000    # Convert keV to MeV


    # Fit Status == 0 -> Good Fit
    # Fit Status == 1 -> Bad Fit

    # Mask for which the fit was bad
    # maskFit = (df['Fit Status'] == 0) # Not currently in use
    # maskFit = (df['Area'] > 1000)
    # for det in range(len(detectors)):
    #
    #     maskDet = ((df['Detector']==detectors[det]) & maskFit)
    #
    #     plt.clf()
    #     plt.errorbar(chanEalpha[maskDet],chanCross[maskDet],yerr=chanCross_err[maskDet],fmt='b.',markersize='2')
    #     plt.yscale('log')
    #     plt.ylim(1e-6,1)
    #     plt.xlim(4,5.6)
    #     plt.xlabel('$E_{\\alpha}$ (MeV)')
    #     plt.ylabel('Cross-Section (barns)')
    #     plt.title('%s %s    %d$^{\circ}$'%(chan,detectors[det],angle[det]))
    #     plt.savefig('yieldPlots/%s/%s_%s.png'%(chan,chan.lower(),detectors[det]),dpi=300)


    # Plot mirror detectors against each other
    # for det in range(0,len(mirrDet),2):
    #     plt.clf()
    #     mask1 = ((df['Detector']==mirrDet[det]) & maskFit)
    #     mask2 = ((df['Detector']==mirrDet[det+1]) & maskFit)
    #     plt.errorbar(chanEalpha[mask1],chanCross[mask1],yerr=chanCross_err[mask1],fmt='b.',markersize='2')
    #     plt.errorbar(chanEalpha[mask2],chanCross[mask2],yerr=chanCross_err[mask2],fmt='k.',markersize='2')
    #     plt.yscale('log')
    #     plt.ylim(5e-5,1e-2)
    #     plt.xlim(4,5.6)
    #     plt.xlabel('$E_{\\alpha}$ (MeV)')
    #     plt.ylabel('Cross-Section (barns)')
    #     plt.title('%s %s    %d$^{\circ}$'%(chan,detectors[det],mirrAng[det]))
    #     plt.savefig('yieldPlots/%s/mirr_%s_%s.png'%(chan,chan.lower(),mirrAng[det]),dpi=300)
    #     plt.show()
    # exit()

    with open("rMatrix/24Mg_rMatrix_%s_allAngles.dat"%chan.lower(),"w") as f:
        for loop in range(len(chanCross)):
            printOut= '%f \t %d \t %.8E \t %.8E \n' %(chanEalpha[loop],Angle[loop],chanCross[loop],chanCross_err[loop])
            f.write(printOut)
    # continue

    AnglesList=['0','15','30','45','90','105','120']

    # test2 = []
    f = open("rMatrix/24Mg_rMatrix_%s.dat"%chan.lower(),"w")
    f.close()
    for ang in AnglesList:
        _chanEalpha = []
        _Angle = []
        _chanCross = []
        _chanCross_err = []

        # Average out over same angle
        for x in range(222):    # total of 222 runs
            _chanEalpha.append(chanEalpha[int(13*x)])
            if ang == '0':
                # print(chanCross[x*13+6])
                _Angle.append(Angle[x*13+6])
                cross = chanCross[x*13+6]*norm0
                _chanCross.append(cross)
                # test2.append(chanCross[x*13+6])
                err = chanCross_err[x*13+6]*norm0
                err = check(cross,err)
                _chanCross_err.append( (err**2+(.05*_chanCross[-1])**2)**.5  )   # inflating the error bars for 5% systematic uncertainty
            elif ang == '15':
                _Angle.append( (abs(Angle[x*13+5])+abs(Angle[x*13+7]))/2 )
                cross = (chanCross[x*13+5]+chanCross[x*13+7])/2*norm15
                _chanCross.append(cross)
                # _chanCross_err.append( (chanCross_err[x*13+5]+chanCross_err[x*13+7])/2 )
                err = (chanCross_err[x*13+5]**2+chanCross_err[x*13+7]**2)**.5*norm15
                err = check(cross,err)
                _chanCross_err.append( (err**2+(.05*_chanCross[-1])**2)**.5 )
            elif ang == '30':
                _Angle.append( int( (abs(Angle[x*13+4])+abs(Angle[x*13+8]))/2 ) )
                cross = (chanCross[x*13+4]+chanCross[x*13+8])/2*norm30
                _chanCross.append(cross)
                err = (chanCross_err[x*13+4]**2+chanCross_err[x*13+8]**2)**.5*norm30
                err = check(cross,err)
                _chanCross_err.append( (err**2+(.05*_chanCross[-1])**2)**.5 )
            elif ang == '45':
                _Angle.append( int( (abs(Angle[x*13+3])+abs(Angle[x*13+9]))/2 ) )
                cross = (chanCross[x*13+3]+chanCross[x*13+9])/2*norm45
                _chanCross.append(cross)
                err = (chanCross_err[x*13+3]**2+chanCross_err[x*13+9]**2)**.5*norm45
                err = check(cross,err)
                _chanCross_err.append( (err**2+(.05*_chanCross[-1])**2)**.5 )
            elif ang == '120':
                # _Angle.append( int( (abs(Angle[x*13])+abs(Angle[x*13+12]))/2 ) )
                _Angle.append( int( 60 ) )
                cross = (chanCross[x*13]+chanCross[x*13+12])/2*norm120
                _chanCross.append(cross)
                err = (chanCross_err[x*13]**2+chanCross_err[x*13+12]**2)**.5*norm120
                err = check(cross,err)
                _chanCross_err.append( (err**2+(.05*_chanCross[-1])**2)**.5 )
            elif ang == '105':
                # _Angle.append( int( (abs(Angle[x*13+1])+abs(Angle[x*13+11]))/2 ) )
                _Angle.append( int( 75 ) )
                cross = (chanCross[x*13+1]+chanCross[x*13+11])/2*norm105
                _chanCross.append(cross)
                err = (chanCross_err[x*13+1]**2+chanCross_err[x*13+11]**2)**.5*norm105
                err = check(cross,err)
                _chanCross_err.append( (err**2+(.05*_chanCross[-1])**2)**.5 )
            elif ang == '90':
                _Angle.append( int( (abs(Angle[x*13+2])+abs(Angle[x*13+10]))/2 ) )
                cross = (chanCross[x*13+2]+chanCross[x*13+10])/2*norm90
                _chanCross.append(cross)
                err = (chanCross_err[x*13+2]**2+chanCross_err[x*13+10]**2)**.5*norm90
                err = check(cross,err)
                _chanCross_err.append( (err**2+(.05*_chanCross[-1])**2)**.5 )




        _chanEalpha = np.array(_chanEalpha)
        _Angle = np.array(_Angle)
        _chanCross = np.array(_chanCross)
        _chanCross_err = np.array(_chanCross_err)

        # Sort by energy, keeping others consistent!
        ind = _chanEalpha.argsort()
        _chanEalpha = _chanEalpha[ind]
        _Angle = _Angle[ind]
        _chanCross = _chanCross[ind]
        _chanCross_err = _chanCross_err[ind]


        # Make the Cross-Section plot
        plt.clf()
        # plt.plot(_chanEalpha,_chanCross)
        plt.errorbar(_chanEalpha,_chanCross,yerr=_chanCross_err,fmt='b.',markersize='2')
        plt.yscale('log')
        # plt.ylim(1e-6,1e-1)
        plt.xlim(4,5.6)
        plt.xlabel('$E_{\\alpha}$ (MeV)', fontsize=14)
        plt.ylabel('Differential Cross-Section (barns/sr)', fontsize=14)
        plt.title('$^{24}$Mg($\\alpha$,p$_{1}\\gamma$) - %s$^{\circ}$'%ang, fontsize=20)
        plt.savefig('crossSection/%s/%s_%s.png'%(chan,chan.lower(),ang),dpi=600)
        plt.clf()


        with open("rMatrix/24Mg_rMatrix_%s.dat"%chan.lower(),"a") as f:
            for loop in range(222):
                printOut= '%f \t %d \t %.8E \t %.8E \n' %(_chanEalpha[loop],_Angle[loop],_chanCross[loop],_chanCross_err[loop])
                f.write(printOut)

    print('DONE!')
