# 24Mg_agp

peakAreasP1.C and peakAreasP2.C fit the p1 and p2 peak respectively. Certain lines need to be un/commented depending on
whether you will be running it locally or on the CRC. For local use you would want to adjust the files that it runs over,
more specifically the for loop would run up to i < 162. This would the run over the 3 included .root files.

p1Areas.sh and p2Areas.sh are used to submit the above macros as a job to be run on the CRC. 

edit_csv.py takes the output from above (peakAreasP1/P2.csv) and also the 24Mg run spreadsheet. It takes the Energy of each run and then appends that to the final P1/P2Yields.csv file

yieldCurve.py is a work in progress. 

