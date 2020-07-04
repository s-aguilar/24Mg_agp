import pandas as import pd
import matplotlib.pyplot as plt


angles = ['0','15','30','45','60','75','90','angInt'])
colNames = ['E_cm','E_ex','phi_cm','fit_XS','fit_S','XS','XS_unc','S','S_unc']

chanNames = ['AZUREOUT_aa=1_R=3.out','AZUREOUT_aa=1_R=4.out','AZUREOUT_aa=1_R=5.out',
            'AZUREOUT_aa=2_R=3.out','AZUREOUT_aa=2_R=4.out','AZUREOUT_aa=2_R=5.out',]
chanDict = {'AZUREOUT_aa=1_R=3.out':'27Al_p1',
            'AZUREOUT_aa=1_R=4.out':'27Al_p2',
            'AZUREOUT_aa=1_R=5.out':'27Al_a1',
            'AZUREOUT_aa=1_R=3.out':'24Mg_p1',
            'AZUREOUT_aa=1_R=3.out':'24Mg_p2',
            'AZUREOUT_aa=1_R=3.out':'24Mg_a1'}

colorDict = {'0':'dimgray',
            '15':'darkviolet',
            '30':'red',
            '45':'crimson',
            '60':'orange',
            '75':'green',
            '90':'blue',
            'angInt':'dodgerBlue'}

data = ['24Mg','27Al']

for ang in angles:
    useColor = colorDict[ang]
    for chan in chanNames:
        df = pd.read_table('data/%s/%s'%(ang,chanNames),names=colNames)
        useLabel = chanDict[chan]

        # Plot data
        plt.errorbar(x=df['E_ex'],y='XS',yerr='XS_unc',c='k')
        # Plot fit
        plt.errorbar(x=df['E_ex'],y='XS',yerr='XS_unc',c=useColor,label=useLabel)
        plt.ylim(1e-7,1e-1)
        plt.xlim()
        plt.title('@useLabel_%s'%ang)
        plt.xlabel('Excitation energy (MeV)')
        plt.ylabel('Cross-Section (b)')
        plt.legend()
        plt.savefig()
