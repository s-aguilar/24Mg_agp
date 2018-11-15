#!/bin/bash

#$ -pe smp 16           # Specify 'parallel env' and number of legally requested cores 'smp 1' NOTE: if more than 10GB memory needed, might nee$
#$ -q debug             # Specify queue type. 'debug'-> 4h/64GB; 'long"-> 15D/256GB; Going shorter means shorter wall-time
#$ -N seabassJob        # Specify job's name
#$ -r y

#export DISPLAY=localhost:0.0
module load root/6.02
module load gcc
#export PATH="/afs/crc.nd.edu/group/nsl/activetarget/users/saguilar/anaconda3/bin:$PATH"
root -b /afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/peakAreas.C
