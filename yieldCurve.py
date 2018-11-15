import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Read in the data into dataframe
df1 = pd.read_csv('peakAreasP1.csv')
df2 = pd.read_csv('peakAreasP2.csv')

print(df1.head())


# Extract the columns of the DataFrame as numpy arrays
p1Run = df1['Run'].values
p1Det = df1['Detector'].values
p1Yield = df1['Yield'].values
p1Yield_err = df1['Yield_err'].values
p1Fit = df1['Fit Status'].values


# Fit Status == 0 -> Good Fit
# Fit Status == 1 -> Bad Fit
#
# Mask for which the fit was bad
maskFit = df1['Fit Status'] == 0





# Extract the values into lists


# for index, row in df.iterrows():
#     r1 = dist(row["x1"],row["y1"],row["xc"],row["yc"])
#     r2 = dist(row["x2"],row["y2"],row["xc"],row["yc"])
#     r3 = dist(row["x3"],row["y3"],row["xc"],row["yc"])
#     rc = np.float64(0.)
#     x.append([row["x1"],row["x2"],row["x3"],row["xc"]])
#     y.append([row["y1"],row["y2"],row["y3"],row["yc"]])
#     ranges.append([r1,r2,r3,rc])
