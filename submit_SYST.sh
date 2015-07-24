#! /bin/bash

for n in {0..100}
do
   qsub -v JOBNUM=${n} -l pmem=3gb scripts/RunSYST.sh
done
