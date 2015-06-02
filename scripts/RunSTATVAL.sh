#! /bin/bash

#PBS -o results_mblstat_20150529/output${JOBNUM}.txt
#PBS -j oe

numberPE=20

cd $PBS_O_WORKDIR
iter=$(($JOBNUM%$numberPE))
./DoFit --run_number $JOBNUM --bootstrap --fit --masspnt 172.5 --mbl --numPE $numberPE --PE $iter --outdir results_mblstat_20150529
