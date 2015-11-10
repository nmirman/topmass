#! /usr/bin/env python

from ROOT import *
import sys

file = TFile( sys.argv[1] )
tree = file.Get('FitResults')

countmt = [0]*7
meanmt = [0]*7
varmt = [0]*7
meanmt_gaus = [0]*7
varmt_gaus = [0]*7

mcmasses = [166.5,169.5,171.5,172.5,173.5,175.5,178.5]

#hmt = []
#for i in range(7):
#   hmt.append( TH1D('hmt'+str(mcmasses[i]),'hmt'+str(mcmasses[i]),50,mcmasses[i]-5,mcmasses[i]+5) )

# set up gaussian fits in ROOFit
topMass = []
topMassMean = []
topMassSigma = []
ds = []
for i in range(7):
   topMass.append( RooRealVar("topMass", "Bootstrap M_{t} @ "+str(mcmasses[i]), mcmasses[i]-2, mcmasses[i]+2, "GeV") )
   topMassMean.append( RooRealVar("M_{t}", "M_{t}", mcmasses[i], mcmasses[i]-2, mcmasses[i]+2, "GeV") )
   topMassSigma.append( RooRealVar("#sigma_{M_{t}}", "#sigma_{M_{t}}", 0.2, 0.0, 5.0, "GeV") )
   ds.append( RooDataSet("topMass", "topMass", RooArgSet(topMass[i])) )

for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   iter=-1
   if tree.mcmass == 166.5:
      iter=0
   elif tree.mcmass == 169.5:
      iter=1
   elif tree.mcmass == 171.5:
      iter=2
   elif tree.mcmass == 172.5:
      iter=3
   elif tree.mcmass == 173.5:
      iter=4
   elif tree.mcmass == 175.5:
      iter=5
   elif tree.mcmass == 178.5:
      iter=6

#   hmt[iter].Fill( tree.mt )
   if tree.fitStatus != 0:
      print 'Fit Status ', tree.fitStatus, ' in PE ', i, ', masspoint ', tree.mcmass
      continue
   if abs(tree.jesfactor-1) > 0.01:
      continue

   countmt[iter] += 1
   meanmt[iter] += tree.mt
   varmt[iter] += (tree.mt)**2

   topMass[iter].setVal(tree.mt)
   ds[iter].add( RooArgSet(topMass[iter]) )

for iter in range(7):
   if countmt[iter] != 0:
      meanmt[iter] = meanmt[iter]/countmt[iter]
      varmt[iter] = varmt[iter]/countmt[iter] - meanmt[iter]**2

# conduct the fit
topgaus = []
fTop = []
#cTop = TCanvas('cbootstrap','cbootstrap',1400,800)
#cTop.Divide(4,2)
for i in range(7):
   topgaus.append( RooGaussian("topgaus"+str(mcmasses[i]), "topgaus"+str(mcmasses[i]), topMass[i], topMassMean[i], topMassSigma[i]) )
   topMassMean[i].setVal(mcmasses[i])
   topMassSigma[i].setVal(0.2)
   topgaus[i].fitTo(ds[i])

   # copy to array of results
   meanmt[i] = topMassMean[i].getVal()
   varmt[i] = topMassSigma[i].getVal()**2

   # display the fit
   #cTop.cd(i+1)
   fTop.append( topMass[i].frame() )
   ds[i].plotOn(fTop[i], RooFit.Binning(40))
   ds[i].statOn(fTop[i], RooFit.Layout( 0.7, 0.97, 0.94),
             RooFit.What("N"))
   topgaus[i].plotOn(fTop[i])
   plotShow = RooArgSet(topMassMean[i], topMassSigma[i])
   topgaus[i].paramOn(fTop[i],RooFit.Parameters(plotShow),
             RooFit.Layout(0.17, 0.48, 0.94),
                 RooFit.Format("NU",RooFit.AutoPrecision(1)))
   if i == 3:
      fTop[i].Draw()

#cTop.Print('results/bootstrap_mt.pdf')

# temp pull plot
hpull = TH1D('hpull', 'Pull for M_{t}(MC) = 172.5;(M_{ti}-#mu)/#sigma_{MIGRAD};Pseudoexperiments', 25, -5, 5)
for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   if tree.mcmass == 172.5:
      hpull.Fill( (tree.mt - meanmt[3])/tree.mt_err )

#cpull = TCanvas('cpull','cpull',800,800)
#hpull.Draw()

#temp
"""
devts = 49000
varmt[0] /= 280000.0/devts
varmt[1] /= 390000.0/devts
varmt[2] /= 260000.0/devts
varmt[3] /= 640000.0/devts
varmt[4] /= 290000.0/devts
varmt[5] /= 440000.0/devts
varmt[6] /= 270000.0/devts
"""

gresults = TGraphErrors()
chi2=0
for i in range(7):
   gresults.SetPoint(i, mcmasses[i], meanmt[i]-mcmasses[i])
   gresults.SetPointError(i, 0.0, sqrt(varmt[i]))
   chi = 0
   if varmt[i] != 0:
      chi = (meanmt[i]-mcmasses[i])/sqrt(varmt[i])
   chi2 += chi*chi
   print str(mcmasses[i])+': '+str(meanmt[i])+' +- '+str(sqrt(varmt[i]))

print 'chi2 = '+str(chi2)
gresults.SetMarkerStyle(20)

fline = TF1('fline','[0]*x',150,200)
fline.SetParameter(0,0.0)
fline.SetLineStyle(7)

canvas = TCanvas('canvas','canvas',800,800)
canvas.cd()

gresults.SetMaximum(1.5)
gresults.SetMinimum(-1.5)
gresults.Draw('AEP')
fline.Draw('same')

raw_input('')
