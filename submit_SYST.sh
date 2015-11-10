#! /bin/bash

for n in {0..20}
#for n in {0..50}
#for n in {0..100}
do
   qsub -v JOBNUM=${n} -l mem=4gb,vmem=4gb,pmem=4gb scripts/RunSYST.sh
done
