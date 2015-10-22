#! /usr/bin/env python

from ROOT import *
import sys
from operator import itemgetter, attrgetter, methodcaller

file = TFile( sys.argv[1] )
tree = file.Get('FitResults')

systs = { 'Central': [0 for i in range(5)] }
for i in range( tree.GetEntries() ):
   tree.GetEntry(i)

   if( tree.fitStatus != 0 ):
      continue

   name = str(tree.syst)

   var = 0
   if 'UP' in name:
      var = 1
      name = name.replace('UP','')
   elif 'DN' in name:
      var = 2
      name = name.replace('DN','')

   if ( name not in systs ):
      systs.update( {name: [0 for i in range(5)]} )

   systs[name][var] = tree.mt

cent = systs['Central'][0]
print( "\n*** Central Value = %5.2f ***\n" % cent )

# fill in missing systematics
slist_temp = [ 'MCscaleup', 'MCscaledown', 'MCmatchingup', 'MCmatchingdown', 'PtTopReweighting',
      'MCTuneP11', 'MCTuneP11mpiHi', 'MCTuneP11TeV', 'MCTuneP11noCR', 'BFRAGrbLEP', 'BFRAGnu', 'JetEnergyResolution', 'METUnclustered', 'PileUp', 'ElectronEnergyScale', 'MuonMomentumScale', 'JESRelativeFSR', 'ElectronId', 'MuonId', 'BTagging', 'JESFlavorPureGluon', 'JESFlavorPureQuark', 'JESFlavorPureCharm', 'JESFlavorPureBottom', 'JESCorrelationGroupMPFInSitu', 'JESCorrelationGroupIntercalibration', 'JESCorrelationGroupUncorrelated' ]
for var in range(50):
   slist_temp.append( 'PDFvar'+str(var) )
print 'Missing systematics: '
for n in slist_temp:
   if ( n not in systs ):
      print ' * ', n
      systs.update( {n: [systs['Central'][0] for i in range(5)]} )

# convert up/down into max/min
for n in systs:
   if( systs[n][0] != 0 ):
      continue
   up = systs[n][1]-cent
   dn = systs[n][2]-cent
   smin = min(up,dn)
   smax = max(up,dn)
   systs[n][3] = smax
   systs[n][4] = smin

# flavor systematics
systs.update( {'JESFlavorTotal': [0 for i in range(5)]} )
for i in range(3,5):
   systs['JESFlavorTotal'][i] = systs['JESFlavorPureGluon'][i] + systs['JESFlavorPureQuark'][i] + \
         systs['JESFlavorPureCharm'][i] + systs['JESFlavorPureBottom'][i]

# jes systematics
slist_jes = [ 'JESCorrelationGroupMPFInSitu', 'JESCorrelationGroupIntercalibration', 'JESCorrelationGroupUncorrelated',
      'JESFlavorTotal' ]
      #'JESFlavorPureGluon', 'JESFlavorPureQuark', 'JESFlavorPureCharm', 'JESFlavorPureBottom' ]

# symmetrize jes systematics
for n in systs:
   up = systs[n][3]
   dn = systs[n][4]
   amax = max( abs(up), abs(dn) )
   if( up < 0 or dn > 0 ):
      systs[n][3] = amax
      systs[n][4] = -1.0*amax

temp = [0.0, 0.0]
for n in slist_jes:
   for i in range(2):
      temp[i] += systs[n][i+3]**2
systs.update( {'JESTotal': [0.0, 0.0, 0.0, sqrt(temp[0]), -1.0*sqrt(temp[1])]} )

#
## MC-related systematics -- require special treatment
#

# scale and matching
up = systs['MCscaleup'][0]
dn = systs['MCscaledown'][0]
systs.update( {'MCscale': [ 0.0, up, dn, max(up-cent,dn-cent), min(up-cent,dn-cent) ]} )

up = systs['MCmatchingup'][0]
dn = systs['MCmatchingdown'][0]
systs.update( {'MCmatching': [ 0.0, up, dn, max(up-cent,dn-cent), min(up-cent,dn-cent) ]} )

# pt top reweighting
up = systs['PtTopReweighting'][1]
systs['PtTopReweighting'] = [0.0, up, 0.0, abs(up-cent), -1.0*abs(up-cent)]

# b fragmentation rbLEP
up = systs['BFRAGrbLEP'][1]
systs['BFRAGrbLEP'] = [0.0, up, 0.0, abs(up-cent), -1.0*abs(up-cent)]
systs.update( {'BFRAG': [ cent, 0.0, 0.0, sqrt(systs['BFRAGnu'][3]**2+systs['BFRAGrbLEP'][3]**2),
   sqrt(systs['BFRAGnu'][4]**2+systs['BFRAGrbLEP'][4]**2) ]} )

# underlying event tunes
central = systs['MCTuneP11'][0]
up = systs['MCTuneP11mpiHi'][0]
dn = systs['MCTuneP11TeV'][0]
systs.update( {'MCEventTunes': [central, up, dn, max(up-central,dn-central), min(up-central,dn-central) ]} )

# color reconnection
up = systs['MCTuneP11noCR'][0]
systs.update( {'MCColorReconnection': [central, up, 0.0, abs(up-central), -1.0*abs(up-central) ]} )

# PDFs
up_sum = 0
dn_sum = 0
for i in range(24):
   up_sum += max( systs['PDFvar'+str(2*i)][0]-cent, systs['PDFvar'+str(2*i+1)][0]-cent )**2
   dn_sum += min( systs['PDFvar'+str(2*i)][0]-cent, systs['PDFvar'+str(2*i+1)][0]-cent )**2
up = sqrt(up_sum)
dn = -1.0*sqrt(dn_sum)
systs.update( {'PDF': [central, 0.0, 0.0, up, dn ]} )

# symmetrize all systematics
for n in systs:
   up = systs[n][3]
   dn = systs[n][4]
   amax = max( abs(up), abs(dn) )
   if( up < 0 or dn > 0 ):
      systs[n][3] = amax
      systs[n][4] = -1.0*amax

# sum up other systematics
slist_all = []
temp = [0.0, 0.0]
slist_all = [ 'JESTotal', 'JetEnergyResolution', 'METUnclustered', 'PileUp', 'ElectronEnergyScale', 'MuonMomentumScale',
      'JESRelativeFSR', 'ElectronId', 'MuonId', 'BTagging', 'PtTopReweighting', 'MCscale', 'MCmatching', 'MCEventTunes', 'PDF',
      'MCColorReconnection', 'BFRAG' ]
for n in slist_all:
   for i in range(2):
      temp[i] += systs[n][i+3]**2

# --------------- print systematics -----------------

# print jes systematics
outformat = '   %-35r: %+3.2f: %+3.2f'
print( "\n JES Systematics:\n" )
for n in slist_jes:
   print( outformat % (n, systs[n][3], systs[n][4]) )
print '   --------------------------------------------------'
print( outformat % ('JESTotal', systs['JESTotal'][3], systs['JESTotal'][4]) )

# print other systematics
print( "\n Other Systematics:\n" )
for n in slist_all:
   print( outformat % (n, systs[n][3], systs[n][4]) )

# total
systs.update( {'TotalSum': [0.0, 0.0, 0.0, sqrt(temp[0]), -1.0*sqrt(temp[1])]} )
print '   --------------------------------------------------'
print( outformat % ('TotalSum', systs['TotalSum'][3], systs['TotalSum'][4]) )

#
## print in latex format
#
texformat = '   %-35r: $^{%+3.2f}_{%+3.2f}$'
print( "\n TeX Format:\n" )
for n in slist_jes:
   print( texformat % (n, systs[n][3], systs[n][4]) )
print '   --------------------------------------------------'
for n in slist_all:
   print( texformat % (n, systs[n][3], systs[n][4]) )
print '   --------------------------------------------------'
print( texformat % ('TotalSum', systs['TotalSum'][3], systs['TotalSum'][4]) )

# dump
outformat = '   %-35r: %+3.2f: %+3.2f: %+3.2f: %+3.2f: %+3.2f'
print( "\n Dump:\n" )
for n in systs:
   shift = [0 for i in range(3)]
   for i in range(3):
      if( systs[n][i] != 0 ):
         shift[i] = cent
   print( outformat % (n, systs[n][0]-shift[0], systs[n][1]-shift[1], systs[n][2]-shift[2], systs[n][3], systs[n][4]) )

print ''

