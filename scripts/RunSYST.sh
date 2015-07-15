#! /bin/bash

#PBS -o results_mblmt2syst_jfact_20150715/output${JOBNUM}.txt
#PBS -j oe

cd $PBS_O_WORKDIR

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
syst[13]=FlavorPureGluonUP
syst[14]=FlavorPureGluonDN
syst[15]=FlavorPureQuarkUP
syst[16]=FlavorPureQuarkDN
syst[17]=FlavorPureCharmUP
syst[18]=FlavorPureCharmDN
syst[19]=FlavorPureBottomUP
syst[20]=FlavorPureBottomDN
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
syst[36]=RelativeFSRUP
syst[37]=RelativeFSRDN
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

./DoFit --run_number ${JOBNUM} --syst ${syst[$JOBNUM]} --fit --masspnt 172.5 --mt2_220 --mbl --outdir results_mblmt2syst_jfact_20150715

