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

# systematic variations
syst[0]=Central
syst[1]=TotalUP
syst[2]=TotalDN
syst[3]=CorrelationGroupMPFInSituUP
syst[4]=CorrelationGroupMPFInSituDN
syst[5]=CorrelationGroupFlavorUP
syst[6]=CorrelationGroupFlavorDN
syst[7]=CorrelationGroupIntercalibrationUP
syst[8]=CorrelationGroupIntercalibrationDN
syst[9]=CorrelationGroupUncorrelatedUP
syst[10]=CorrelationGroupUncorrelatedDN
syst[11]=CorrelationGroupbJESUP
syst[12]=CorrelationGroupbJESDN
syst[13]=MCscaleup
syst[14]=MCscaledown
syst[15]=MCmatchingup
syst[16]=MCmatchingdown
syst[17]=JetEnergyResolutionUP
syst[18]=JetEnergyResolutionDN
syst[19]=METUnclusteredUP
syst[20]=METUnclusteredDN
syst[21]=PileUpUP
syst[22]=PileUpDN
syst[23]=ElectronEnergyScaleUP
syst[24]=ElectronEnergyScaleDN
syst[25]=MuonMomentumScaleUP
syst[26]=MuonMomentumScaleDN
syst[27]=PtTopReweightingUP

job=$(($1-1))

./DoFit --syst ${syst[$job]} --fit --masspnt 172.5 --mbl

cd ..
mv topmass/fitresults.root .
mv topmass/plotsFitResults.root .

echo "End of job on " `date`
