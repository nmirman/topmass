#! /bin/bash

#cd /uscms/home/nmirman/CMSSW_5_3_9_patch3/src
#eval `scramv1 runtime -sh`

export WORKING_DIR=$_CONDOR_SCRATCH_DIR/CMSSW_5_3_9_patch3

rnum=$((3*$1))

# mc top masses
mass=(161.5 163.5 166.5 169.5 172.5 175.5 178.5 181.5)

#iter=$(($1/100))
iter=$((($1-1)/50))
#iter = 4
#iter=$(($1-1))
#echo ${mass[$iter]} $lbound $rbound
penum=$((($1-1)/100))


echo $1
echo $iter
echo ${mass[$iter]}

#pwd
tar -xvzf forgrid.tar.gz
#mkdir results

#cd /uscms/home/nmirman/topmass/
#./DoFit --run_number $1 --bootstrap --fit --masspnt ${mass[$iter]}
#./DoFit --run_number $1 --bootstrap --fit --masspnt ${mass[$iter]} --gnorm 30000
#./DoFit --run_number $1 --bootstrap --fit --masspnt ${mass[$iter]} --lbnd 200 --rbnd 300
#./DoFit --fit --masspnt ${mass[$iter]} --lbnd 0 --rbnd 30
#./DoFit --fit --masspnt ${mass[$1]} --lmbl 13 --lmt 32

#./DoFit --templates --maos220 --l_maos220 $rnum
#./DoFit --templates --mt2
#./DoFit --learnparams --maos220 --l_maos220 15
#./DoFit --fit --maos220 --l_maos220 15 --masspnt ${mass[$iter]}
#./DoFit --fit --maos210 --masspnt ${mass[$iter]}

#./DoFit --run_number $1 --bootstrap --fit --maos210 --maoscuts210 0 --masspnt ${mass[$iter]}
##./DoFit --run_number $1 --bootstrap --numPE 20 --PE $penum --fit --maos210 --maoscuts210 4 --masspnt ${mass[4]}
##./DoFit --run_number $1 --bootstrap --fit --mbl --mt2_220 --maos210 --maoscuts210 4 --maos220 --maoscuts220 4 --masspnt ${mass[$iter]}
./DoFit --run_number $1 --bootstrap --fit --mt2_220 --masspnt ${mass[4]}
#./DoFit --fit --maos210 --maoscuts210 3 --masspnt ${mass[$iter]}


#./DoFit --templates --maos220 --maoscuts220 5

#ls

#cp plotsTemplates.root ..

#mv results/plotsTemplates.root ./plotsTemplates.root
