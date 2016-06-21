#! /bin/bash

#PBS -o results_2Dfit_statval_20160222/output${JOBNUM}.txt
#PBS -j oe

numberPE=20

cd $PBS_O_WORKDIR
iter=$(($JOBNUM%$numberPE))
./DoFit --run_number $JOBNUM --bootstrap --fit --masspnt 172.5 --mbl --mt2_221 --numPE $numberPE --PE $iter --jfactor --outdir results_2Dfit_statval_20160222/
