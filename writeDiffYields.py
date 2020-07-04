import pandas as pd

channels = ['p1','p2']
colNames = ['Energy','Angle','Cross-section','Error']
angles = [0,15,30,45,60,75,90]

for ch in channels:
    print('Working on %s...'%ch)

    # Read in the data into dataframe
    df = pd.read_table('rMatrix/24Mg_rMatrix_%s.dat'%ch,names=colNames)

    if ch == 'p1':
        # df = df.query('Energy > 3.74')# COM
        df = df.query('Energy > 4.37')

    # Make different sheets for each detector
    with pd.ExcelWriter('Yields/%s/%sDiffYieldsPerAngle.xlsx'%(ch.upper(),ch)) as writer:
        for ang in angles:
            tempdf = df.query('Angle == @ang')
            tempdf.to_excel(writer,sheet_name='%s'%ang,index=False)
