#! /bin/bash

#for n in {0..100}
for n in {0..20}
do
   qsub -v JOBNUM=${n} -l pmem=5gb scripts/RunSYST.sh
done
