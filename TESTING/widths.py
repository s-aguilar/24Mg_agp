import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

detectors = ['det_h0-0','det_h0-1','det_h0-2','det_h0-3','det_h0-4','det_h0-5',
            'det_h0-6','det_h0-7','det_h1-0','det_h1-1','det_h1-2','det_h1-3',
            'det_h1-4']

# Read in the data into dataframe
df1 = pd.read_csv('Yields/P1/p1Yields.csv')
df1 = df1.sort_values(by=['Ea'])
df2 = pd.read_csv('Yields/P2/p2Yields.csv')
df2 = df2.sort_values(by=['Ea'])
df3 = pd.read_csv('Yields/A1/a1Yields.csv')
df3 = df3.sort_values(by=['Ea'])

# Check widths at same angle run by run
p1_widths = df1['sig1'].values
p1Ealpha = df1['Ea'].values/1000    # Convert keV to MeV

p2_widths = df2['sig1'].values
p2Ealpha = df2['Ea'].values/1000    # Convert keV to MeV

a1_widths = df3['sig1'].values
a1Ealpha = df3['Ea'].values/1000    # Convert keV to MeV

# Angle (deg) for each detector, negative is beam left, positive is beam right
#         00  01 02 03 04 05 06 07  10  11  12   13   14
angle = np.array([120,105,90,45,30,15,0,-15,-30,-45,-90,-105,-120])

for det in range(len(detectors)):

    maskDet = ((df1['Detector']==detectors[det])) # & mask1Fit & )
    plt.scatter(p1Ealpha[maskDet],p1_widths[maskDet],c='b',marker='.')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Width (keV)')
    plt.title('p1 widths %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('widths/P1/p1_%s.png'%detectors[det],dpi=300)
    plt.clf()

for det in range(len(detectors)):

    maskDet = ((df2['Detector']==detectors[det])) # & mask1Fit & )
    plt.scatter(p2Ealpha[maskDet],p2_widths[maskDet],c='b',marker='.')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Width (keV)')
    plt.title('p2 widths %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('widths/P2/p2_%s.png'%detectors[det],dpi=300)
    plt.clf()

for det in range(len(detectors)):

    maskDet = ((df3['Detector']==detectors[det])) # & mask1Fit & )
    plt.scatter(a1Ealpha[maskDet],a1_widths[maskDet],c='b',marker='.')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Width (keV)')
    plt.title('a1 widths %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('widths/A1/a1_%s.png'%detectors[det],dpi=300)
    plt.clf()




# Resolution ()%)
p1_res = 2.35*p1_widths/843.76*100
p2_res = 2.35*p2_widths/1014.52*100
a1_res = 2.35*a1_widths/1368.63*100

for det in range(len(detectors)):

    maskDet = ((df1['Detector']==detectors[det])) # & mask1Fit & )
    plt.scatter(p1Ealpha[maskDet],p1_res[maskDet],c='b',marker='.')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Resolution (%)')
    plt.title('p1 Resolution %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('resolution/P1/p1_%s.png'%detectors[det],dpi=300)
    plt.clf()

for det in range(len(detectors)):

    maskDet = ((df2['Detector']==detectors[det])) # & mask1Fit & )
    plt.scatter(p2Ealpha[maskDet],p2_res[maskDet],c='b',marker='.')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Resolution (%)')
    plt.title('p2 Resolution %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('resolution/P2/p2_%s.png'%detectors[det],dpi=300)
    plt.clf()

for det in range(len(detectors)):

    maskDet = ((df3['Detector']==detectors[det])) # & mask1Fit & )
    plt.scatter(p1Ealpha[maskDet],a1_res[maskDet],c='b',marker='.')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Resolution (%)')
    plt.title('a1 Resolution %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('resolution/A1/a1_%s.png'%detectors[det],dpi=300)
    plt.clf()





"""
# Check gain match factor at each angle run by run
p1_a = df1['a'].values
p2_a = df2['a'].values
a1_a = df3['a'].values

for det in range(len(detectors)):

    maskDet = ((df1['Detector']==detectors[det])) # & mask1Fit & )
    plt.scatter(p1Ealpha[maskDet],p1_a[maskDet],c='b',marker='.')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Gain Match')
    plt.ylim(.98,1.02)
    plt.title('p1 a %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('gainMatch/P1/p1_%s.png'%detectors[det],dpi=300)
    plt.clf()

for det in range(len(detectors)):

    maskDet = ((df2['Detector']==detectors[det])) # & mask1Fit & )
    plt.scatter(p2Ealpha[maskDet],p2_a[maskDet],c='b',marker='.')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Gain Match')
    plt.ylim(.98,1.02)
    plt.title('p2 a %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('gainMatch/P2/p2_%s.png'%detectors[det],dpi=300)
    plt.clf()

for det in range(len(detectors)):

    maskDet = ((df3['Detector']==detectors[det])) # & mask1Fit & )
    plt.scatter(a1Ealpha[maskDet],a1_a[maskDet],c='b',marker='.')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Gain Match')
    plt.ylim(.98,1.02)
    plt.title('a1 a %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('gainMatch/A1/a1_%s.png'%detectors[det],dpi=300)
    plt.clf()
# """
