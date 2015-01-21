#! /usr/bin/env python

from ROOT import *
import sys

file = TFile( sys.argv[1] )
tree = file.Get('FitResults')
title = sys.argv[2]

countmt = [0]*20
meanmt = [0]*20
varmt = [0]*20
meanmt_gaus = [0]*20
varmt_gaus = [0]*20

mcmasses = [161.5,163.5,166.5,169.5,172.5,175.5,178.5,181.5]

#hmt = []
#for i in range(8):
#   hmt.append( TH1D('hmt'+str(mcmasses[i]),'hmt'+str(mcmasses[i]),50,mcmasses[i]-5,mcmasses[i]+5) )

# set up gaussian fits in ROOFit
#topMass = []
#topMassMean = []
#topMassSigma = []
#ds = []
#for i in range(20):
#   topMass.append( RooRealVar("topMass", "Bootstrap M_{t} @ "+str(i), i-2, i+2, "GeV") )
#   topMassMean.append( RooRealVar("M_{t}", "M_{t}", i, i-2, mcmasses[i]+2, "GeV") )
#   topMassSigma.append( RooRealVar("#sigma_{M_{t}}", "#sigma_{M_{t}}", 0.2, 0.0, 5.0, "GeV") )
#   ds.append( RooDataSet("topMass", "topMass", RooArgSet(topMass[i])) )

for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   iter=tree.indexPE

#   hmt[iter].Fill( tree.mt )
   countmt[iter] += 1
   meanmt[iter] += tree.mt
   varmt[iter] += (tree.mt)**2

#   topMass[iter].setVal(tree.mt)
#   ds[iter].add( RooArgSet(topMass[iter]) )

for iter in range(20):
   if countmt[iter] != 0:
      meanmt[iter] = meanmt[iter]/countmt[iter]
      varmt[iter] = varmt[iter]/countmt[iter] - meanmt[iter]**2

# conduct the fit
#topgaus = []
#fTop = []
#cTop = TCanvas('cbootstrap','cbootstrap',1400,800)
#cTop.Divide(4,2)
#for i in range(20):
#   topgaus.append( RooGaussian("topgaus"+str(i), "topgaus"+str(i), topMass[i], topMassMean[i], topMassSigma[i]) )
#   topMassMean[i].setVal(i)
#   topMassSigma[i].setVal(0.2)
#   topgaus[i].fitTo(ds[i])

   # copy to array of results
#   meanmt[i] = topMassMean[i].getVal()
#   varmt[i] = topMassSigma[i].getVal()**2

   # display the fit
#   cTop.cd(i+1)
#   fTop.append( topMass[i].frame() )
#   ds[i].plotOn(fTop[i], RooFit.Binning(40))
#   ds[i].statOn(fTop[i], RooFit.Layout( 0.7, 0.97, 0.94),
#             RooFit.What("N"))
#   topgaus[i].plotOn(fTop[i])
#   plotShow = RooArgSet(topMassMean[i], topMassSigma[i])
#   topgaus[i].paramOn(fTop[i],RooFit.Parameters(plotShow),
#             RooFit.Layout(0.17, 0.48, 0.94),
#                 RooFit.Format("NU",RooFit.AutoPrecision(1)))
   #fTop[i].Draw()

#cTop.Print('results/bootstrap_mt.pdf')

# temp pull plot
hpull = TH1D('hpull', 'Pull for M_{t}(MC) = 172.5, '+title+';(M_{ti}-172.5)/#sigma;Pseudoexperiments', 25, -5, 5)
for i in range( 20 ):
   hpull.Fill( (meanmt[i] - 172.5)/varmt[i] )

cpull = TCanvas('cpull','cpull',800,800)
hpull.Draw()

#spread plot
#mmin = 100
#mmax = -100
#for j in range( tree.GetEntries() ):
#   tree.GetEntry(j)
#   if tree.mcmass == 172.5:
#      if (tree.mt - 172.5) < mmin:
#         mmin = tree.mt - 172.5
#      if (tree.mt - 172.5) > mmax:
#         mmax = tree.mt - 172.5

#hspread = TH1D('hspread', 'Spread for M_{t}(MC) = 172.5, '+title+';(M_{ti}-172.5);Pseudoexperiments', 25, mmin - 0.5, mmax + 0.5)
#for i in range( tree.GetEntries() ):
#   tree.GetEntry(i)
#   if tree.mcmass == 172.5:
#      hspread.Fill( tree.mt - 172.5 )

#cspread = TCanvas('cspread','cspread',800,800)
#hspread.Draw()

#gresults = TGraphErrors()
#chi2=0
#for i in range(8):
#   gresults.SetPoint(i, mcmasses[i], meanmt[i]-mcmasses[i])
#   gresults.SetPointError(i, 0.0, sqrt(varmt[i]))
#   chi = 0
#   if varmt[i] != 0:
#      chi = (meanmt[i]-mcmasses[i])/sqrt(varmt[i])
#   chi2 += chi*chi
#   print str(mcmasses[i])+': '+str(meanmt[i])+' +- '+str(sqrt(varmt[i]))

#print 'chi2 = '+str(chi2)
#gresults.SetMarkerStyle(20)
#gresults.SetTitle(title)
#gresults.GetXaxis().SetTitle('GEN Top Mass (GeV)')
#gresults.GetYaxis().SetTitle('MEAS-GEN Top Mass (GeV)')


#fline = TF1('fline','[0]*x',150,200)
#fline.SetParameter(0,0.0)
#fline.SetLineStyle(7)

#canvas = TCanvas('canvas','canvas',800,800)
#canvas.cd()

#gresults.SetMaximum(1.5)
#gresults.SetMinimum(-1.5)
#gresults.Draw('AEP')
#fline.Draw('same')

raw_input('')
