#! /bin/bash

#PBS -e logs/error_staval_${JOBNUM}.txt
#PBS -o logs/output_statval_${JOBNUM}.txt

numberPE=20

cd $PBS_O_WORKDIR
iter=$(($JOBNUM%$numberPE))
./DoFit --run_number $JOBNUM --bootstrap --fit --masspnt 172.5 --mbl --numPE $numberPE --PE $iter --outfile results_mblstatval/fitresults_${JOBNUM}.root
