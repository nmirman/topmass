#! /usr/bin/env python

from ROOT import *
import CMS_lumi, tdrstyle
import sys

tdrstyle.setTDRStyle()

iPeriod = 0
iPos = 11
CMS_lumi.lumi_sqrtS = '8 TeV'
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = 'Simulation Preliminary'

file = TFile( sys.argv[1] )
tree = file.Get('FitResults')

tree.GetEntry(0)
numPE = tree.statval_numPE

print 'Running on ', numPE, ' pseudoexperiments.'

mcmasses = [166.5,169.5,171.5,172.5,173.5,175.5,178.5]
pullmean = [0]*7
pullmeanerr = [0]*7
pullsigma = [0]*7
pullsigmaerr = [0]*7

#for m in range(7):
m = 3
mcmass = mcmasses[m]

countmt = [0]*numPE
meanmt = [0]*numPE
sigmamt = [0]*numPE
meanmt_gaus = [0]*numPE
sigmamt_gaus = [0]*numPE

# set up gaussian fits in ROOFit
topMass = []
topMassMean = []
topMassSigma = []
ds = []
for i in range(numPE):
   topMass.append( RooRealVar("topMass", "Bootstrap M_{t} @ PE "+str(i), mcmass-10, mcmass+10, "GeV") )
   topMassMean.append( RooRealVar("M_{t}", "M_{t}", mcmass, mcmass-10, mcmass+10, "GeV") )
   topMassSigma.append( RooRealVar("#sigma_{M_{t}}", "#sigma_{M_{t}}", 0.5, 0.0, 5.0, "GeV") )
   ds.append( RooDataSet("topMass", "topMass", RooArgSet(topMass[i])) )

for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   if tree.mcmass == mcmass and tree.fitStatus == 0:
      iter = tree.statval_PE
      topMass[iter].setVal(tree.mt)
      ds[iter].add( RooArgSet(topMass[iter]) )


# conduct the fit
topgaus = []*numPE
for i in range(numPE):
   topgaus.append( RooGaussian("topgaus"+str(i), "topgaus"+str(i), topMass[i], topMassMean[i], topMassSigma[i]) )
   topMassMean[i].setVal(mcmass)
   topMassSigma[i].setVal(0.2)
   topgaus[i].fitTo(ds[i])

   # copy to array of results
   meanmt[i] = topMassMean[i].getVal()
   sigmamt[i] = topMassSigma[i].getVal()


# now make a pull plot
pull = RooRealVar("pull", "pull", -5, 5, "")
pullMean = RooRealVar("#mu_{pull}", "pullmean", 0, -3, 3, "")
pullSigma = RooRealVar("#sigma_{pull}", "pullsigma", 1, 0, 3, "")
pulldat = RooDataSet("pull", "pull", RooArgSet(pull) )

for i in range(numPE):
   pull.setVal( (meanmt[i]-mcmass)/sigmamt[i] )
   pulldat.add( RooArgSet(pull) )
   print 'Pull ', i, ' ', meanmt[i], sigmamt[i], ': ', (meanmt[i]-mcmass)/sigmamt[i]

pullgaus = RooGaussian("pullgaus", "pullgaus", pull, pullMean, pullSigma)
pullMean.setVal(0)
pullSigma.setVal(1)
pullgaus.fitTo(pulldat)

# display the fit
fPull = pull.frame()
pulldat.plotOn(fPull, RooFit.Binning(20))
#pulldat.statOn(fPull, RooFit.Layout(0.7, 0.97, 0.94), RooFit.What("N"))
pullgaus.plotOn(fPull)
plotShow = RooArgSet(pullMean, pullSigma)
pullgaus.paramOn(fPull, RooFit.Parameters(plotShow), RooFit.Layout(0.63, 0.95, 0.91),
      RooFit.Format("NELU",RooFit.AutoPrecision(2)))
fPull.getAttLine().SetBorderSize(0)

c = TCanvas('c','c',800,800)
fPull.Draw()

CMS_lumi.CMS_lumi(c, iPeriod, iPos)
c.Draw()
c.Print('pull_mt.pdf')

raw_input('...')
sys.exit()

pullmean[m] = pullMean.getVal()
pullmeanerr[m] = pullMean.getError()
pullsigma[m] = pullSigma.getVal()
pullsigmaerr[m] = pullSigma.getError()

gmean = TGraphErrors()
gsigma = TGraphErrors()
gmean.GetXaxis().SetTitle('MC Mass (GeV)')
gsigma.GetXaxis().SetTitle('MC Mass (GeV)')
gmean.GetYaxis().SetTitle('Pull Mean')
gsigma.GetYaxis().SetTitle('Pull Width')
gmean.SetMaximum(3)
gmean.SetMinimum(-3)
gsigma.SetMaximum(3)
gsigma.SetMinimum(-3)
gmean.SetMarkerStyle(20)
gsigma.SetMarkerStyle(20)
for i in range(7):
   gmean.SetPoint(i, mcmasses[i], pullmean[i])
   gmean.SetPointError(i, 0.0, pullmeanerr[i])
   gsigma.SetPoint(i, mcmasses[i], pullsigma[i])
   gsigma.SetPointError(i, 0.0, pullsigmaerr[i])

canvas = TCanvas('canvas','canvas',1400,700)
canvas.Divide(2,1)
canvas.cd(1)
gmean.Draw('AP')
canvas.cd(2)
gsigma.Draw('AP')

raw_input('')
