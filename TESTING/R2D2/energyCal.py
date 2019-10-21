import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# energy calibration
energyCal = [4,5.4]

# Peak positions

# Old calibration
adc0 = [2061.33,2787.95]
adc1 = [2024.12,2751.88]
adc2 =[2039.76,2770.02]
adc3 = [2025.86,2742.33]
adc4 = [2038.6,2764.17]
adc5 = [2015.63,2735.37]
adc6 = [2022.95,2745.21]
adc9 = [2042.8,2779.41]
adc10 = [2066.06,2806.41]
adc11 = [2046.48,2775.62]
adc12 = [2069.07,2802.81]
adc13 = [2060.92,2795.51]
adc14 = [2077.52,2817.09]
adc15 = [2083.46,2821.86]
adc27 = [2085.89,2823.05]
adc28 = [2127.03,2879.41]
adc29 = [2121.81,2871.05]
adc30 = [2185.44,2954.69]
adc31 = [2163.91,2920.11]

# New calibration
adc0 = [2061.33,2787.95]
adc1 = [2024.12,2751.88]
adc2 =[2039.76,2770.02]
adc3 = [2025.86,2742.33]
adc4 = [2038.6,2764.17]
adc5 = [2015.63,2735.37]
adc6 = [2022.95,2745.21]
adc9 = [2042.8,2779.41]
adc10 = [2066.06,2806.41]
adc11 = [2046.48,2775.62]
adc12 = [2069.07,2802.81]
adc13 = [2060.92,2795.51]
adc14 = [2077.52,2817.09]
adc15 = [2083.46,2821.86]
adc27 = [2085.89,2823.05]
adc28 = [2127.03,2879.41]
adc29 = [2121.81,2871.05]
adc30 = [2185.44,2954.69]
adc31 = [2163.91,2920.11]

# Change for each adc you want
channelToEnergy = interp1d(adc29,energyCal,fill_value='extrapolate')
energyToChannel = interp1d(energyCal,adc29,fill_value='extrapolate')

# energy or channel you want to find
energy = 5.2
channel = 2450

print('\nChannel %d'% channel,'is:\t%4.3f' % channelToEnergy(channel),'MeV')
print('\nEnergy %3.2f'% energy,'MeV','is:\t%5.2f' % energyToChannel(energy),'Channel\n')

# plt.plot(adc2,energyCal)
# plt.show()
