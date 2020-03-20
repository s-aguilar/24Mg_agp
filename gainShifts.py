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

runs1 = np.array([int(x[-3:]) for x in df1['Run']])
# runs2 = np.array([int(x[-3:]) for x in df2['Run']])
# runs3 = np.array([int(x[-3:]) for x in df3['Run']])

# print(len(runs1),runs1)

p1Ealpha = df1['Ea'].to_numpy()/1000    # Convert keV to MeV
# p2Ealpha = df2['Ea'].to_numpy()/1000    # Convert keV to MeV
# a1Ealpha = df3['Ea'].to_numpy()/1000    # Convert keV to MeV

# Check gain match factor at each angle run by run
p1_a = df1['a'].to_numpy()
# p2_a = df2['a'].to_numpy()
# a1_a = df3['a'].to_numpy()

print(len(p1_a),len(p1_a[p1_a>0]))


# Angle (deg) for each detector, negative is beam left, positive is beam right
#         00  01 02 03 04 05 06 07  10  11  12   13   14
angle = np.array([120,105,90,45,30,15,0,-15,-30,-45,-90,-105,-120])

for det in range(len(detectors)):

    maskDet = ((df1['Detector']==detectors[det]) & (p1_a > 0.1))
    # plt.scatter(p1Ealpha[maskDet],p1_a[maskDet],c='b',marker='.')
    plt.scatter(runs1[maskDet],p1_a[maskDet],c='b',marker='.')
    plt.plot(runs1[maskDet],p1_a[maskDet],c='k',alpha=.5)
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Gain Match')
    plt.ylim(0.995,1.010)
    plt.title('p1 a %s    %d$^{\circ}$'%(detectors[det],angle[det]))
    plt.savefig('gainMatch/P1/p1_%s.png'%detectors[det],dpi=300)
    plt.clf()

# for det in range(len(detectors)):
#
#     maskDet = ((df2['Detector']==detectors[det])) # & mask1Fit & )
#     # plt.scatter(p2Ealpha[maskDet],p2_a[maskDet],c='b',marker='.')
#     plt.scatter(np.array(runs2)[maskDet],p2_a[maskDet],c='b',marker='.')
#     plt.xlabel('$E_{\\alpha}$ (MeV)')
#     plt.ylabel('Gain Match')
#     plt.ylim(-2,2)
#     plt.title('p2 a %s    %d$^{\circ}$'%(detectors[det],angle[det]))
#     plt.savefig('gainMatch/P2/p2_%s.png'%detectors[det],dpi=300)
#     plt.clf()
#
# for det in range(len(detectors)):
#
#     maskDet = ((df3['Detector']==detectors[det])) # & mask1Fit & )
#     # plt.scatter(a1Ealpha[maskDet],a1_a[maskDet],c='b',marker='.')
#     plt.scatter(np.array(runs3)[maskDet],a1_a[maskDet],c='b',marker='.')
#     plt.xlabel('$E_{\\alpha}$ (MeV)')
#     plt.ylabel('Gain Match')
#     plt.ylim(-2,2)
#     plt.title('a1 a %s    %d$^{\circ}$'%(detectors[det],angle[det]))
#     plt.savefig('gainMatch/A1/a1_%s.png'%detectors[det],dpi=300)
#     plt.clf()
