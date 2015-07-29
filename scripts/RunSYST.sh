#! /bin/bash

#PBS -o results_mbl221syst9_jfact_20150728/output${JOBNUM}.txt
#PBS -j oe

cd $PBS_O_WORKDIR

# systematic variations
syst[0]=Central
syst[1]=JESTotalUP
syst[2]=JESTotalDN
syst[3]=JESCorrelationGroupMPFInSituUP
syst[4]=JESCorrelationGroupMPFInSituDN
syst[5]=JESCorrelationGroupFlavorUP
syst[6]=JESCorrelationGroupFlavorDN
syst[7]=JESCorrelationGroupIntercalibrationUP
syst[8]=JESCorrelationGroupIntercalibrationDN
syst[9]=JESCorrelationGroupUncorrelatedUP
syst[10]=JESCorrelationGroupUncorrelatedDN
syst[11]=JESCorrelationGroupbJESUP
syst[12]=JESCorrelationGroupbJESDN
syst[13]=JESFlavorPureGluonUP
syst[14]=JESFlavorPureGluonDN
syst[15]=JESFlavorPureQuarkUP
syst[16]=JESFlavorPureQuarkDN
syst[17]=JESFlavorPureCharmUP
syst[18]=JESFlavorPureCharmDN
syst[19]=JESFlavorPureBottomUP
syst[20]=JESFlavorPureBottomDN
syst[21]=MCscaleup
syst[22]=MCscaledown
syst[23]=MCmatchingup
syst[24]=MCmatchingdown
syst[25]=JetEnergyResolutionUP
syst[26]=JetEnergyResolutionDN
syst[27]=METUnclusteredUP
syst[28]=METUnclusteredDN
syst[29]=PileUpUP
syst[30]=PileUpDN
syst[31]=ElectronEnergyScaleUP
syst[32]=ElectronEnergyScaleDN
syst[33]=MuonMomentumScaleUP
syst[34]=MuonMomentumScaleDN
syst[35]=PtTopReweightingUP
syst[36]=JESRelativeFSRUP
syst[37]=JESRelativeFSRDN
syst[38]=ElectronIdUP
syst[39]=ElectronIdDN
syst[40]=MuonIdUP
syst[41]=MuonIdDN
syst[42]=BTaggingUP
syst[43]=BTaggingDN
syst[44]=MCTuneP11
syst[45]=MCTuneP11TeV
syst[46]=MCTuneP11mpiHi
syst[47]=MCTuneP11noCR
syst[48]=BFRAGnuUP
syst[49]=BFRAGnuDN
syst[50]=BFRAGrbLEPUP
ipdf=51
for i in {0..49}
do
   syst[ipdf+i]=PDFvar$i
done

./DoFit --run_number ${JOBNUM} --syst ${syst[$JOBNUM]} --fit --masspnt 172.5 --mt2_221 --mbl --jfactor --outdir results_mbl221syst9_jfact_20150728
