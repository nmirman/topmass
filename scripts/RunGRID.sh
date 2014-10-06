#! /bin/bash

set -o verbose

echo "Starting job on " `date`
echo "Running on " `uname -a`
echo "System release " `cat /etc/redhat-release`

# setup ROOT
cd src
eval `scramv1 runtime -sh`
cd ..
echo $LD_LIBRARY_PATH

# unpack tarball
mkdir topmass/
mv topmassforgrid.tar.gz topmass/
cd topmass/
tar -xvzf  topmassforgrid.tar.gz

export WORKING_DIR=$(pwd)
echo $WORKING_DIR

# mc top masses
mass=(161.5 163.5 166.5 169.5 172.5 175.5 178.5 181.5)

job=$(($1-1))
iter=$(($job/50))

./DoFit --run_number $job --bootstrap --fit --masspnt ${mass[$iter]} --maos220

cd ..
mv topmass/fitresults.root .
mv topmass/plotsFitResults.root .

echo "End of job on " `date`
