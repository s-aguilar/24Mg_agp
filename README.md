# 24Mg_agp

p1Yields.C and p2Yields.C fit the p1 and p2 peak respectively. There is also the a1Yield.C file which fits the a1 peak. The structure of all the files is the same.


Certain lines need to be un/commented depending on whether you will be running it locally or on the CRC. The directory in which it looks for the runs must be changed (in this case comment out the crc directory)
    For local use you would want to adjust the files that it runs over,
    more specifically the for loop that occurs in the main body would run up to i < 162. This then runs over the 3 included .root files.


In these xxYields.C files the main body, the outer loop - loops over the runs, and the inner loop - loops over the detectors in that run. The peakFitter subroutine is what actually gets run in each loop cycle. There are some global variables that are created because within each function call of peakFitter I need to store some values so that they can then be used in subsequent runs (outer loops) [probably a better/cleaner way of doing this but it works for now]. Also some files are created which will store the peak ranges in the BG and run spectra so that we can track how the gains drift run by run.

The general procedure is as follows: For longer runs, fit the BG peaks that are present in both BG and run spectra. Perform a gain match of the BG to the run spectra and store those gain match calibration spectra to be used for the shorter runs. Subtract out the BG. Then find the peak centroids (p1, p2, a1) and calibrate the run spectra from channel to energy and store those calibration constants to be used for the shorter runs. Once calibrated constrain the fit to the energy of the peak (p1, p2, a1). The first fit parameter is the peak area which is then divided by charge to get the yield and store it into csv file.


edit_csv.py takes the output from above (Yields/P1/_P1.csv) and also the 24Mg run spreadsheet. It takes the Energy of each run and then appends that to the final column of the _P1.csv file


crossSection.py calculates the respective yields to differential cross-sections and then plots it as a function of E_alpha.
