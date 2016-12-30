#! /usr/bin/env python

from ROOT import *
import CMS_lumi, tdrstyle
import sys

tdrstyle.setTDRStyle()
gStyle.SetPadTopMargin(0.06)

iPeriod = 0
iPos = 11
CMS_lumi.lumi_sqrtS = '8 TeV'
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = 'Simulation'
fit_type = '1D fit'

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
   topMass.append( RooRealVar("topMass", "M_{t}", mcmasses[i]-2, mcmasses[i]+2, "GeV") )
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
   #if abs(tree.jesfactor-1) > 0.02:
   #   continue
   if tree.fitStatus != 0:
      continue

   mode = sys.argv[2]
   whyb = 0.8
   if mode == '1d':
      mt = tree.mt_fix
   elif mode == '2d':
      mt = tree.mt
      fit_type = '2D fit'
   elif mode == 'hyb':
      mt = (1-whyb)*tree.mt + whyb*tree.mt_fix
      fit_type = 'Hybrid fit'
   elif mode == 'maos':
      mt = tree.mt
      fit_type = 'MAOS fit'
   elif mode == 'jsf':
      mt = tree.jesfactor
      fit_type = '2D fit'
   else:
      print 'ERROR: fit mode not recognized!'
      sys.exit()

   countmt[iter] += 1
   meanmt[iter] += mt
   varmt[iter] += mt**2

   topMass[iter].setVal(mt)
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
   #ds[i].statOn(fTop[i], RooFit.Layout( 0.7, 0.97, 0.94),
   #          RooFit.What("N"))
   topgaus[i].plotOn(fTop[i])
   #plotShow = RooArgSet(topMassMean[i], topMassSigma[i])
   #topgaus[i].paramOn(fTop[i],RooFit.Parameters(plotShow),
   #          RooFit.Layout(0.647, 0.9, 0.85),
   #              RooFit.Format("NU",RooFit.AutoPrecision(1)))
   #fTop[i].getAttLine().SetBorderSize(0)
   fTop[i].SetMinimum(1e-5)
   if i == 3:
      cTop = TCanvas('cboot','cboot',800,800)
      fTop[i].GetYaxis().SetTitle('pseudo-experiments per 0.1 GeV')
      fTop[i].Draw()

CMS_lumi.CMS_lumi(cTop, iPeriod, iPos)
latex = TLatex()
latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.04)
latex.DrawLatex(0.8, 0.887, fit_type)
latex.DrawLatex(0.66, 0.80, 'M_{t} = '+'%.2f'%meanmt[3]+' GeV')
latex.DrawLatex(0.68, 0.75, '#sigma_{M_{t}} = '+'%.2f'%sqrt(varmt[3])+' GeV')
cTop.Print('boot_mt.pdf')

# temp pull plot
hpull = TH1D('hpull', 'Pull for M_{t}(MC) = 172.5;(M_{ti}-#mu)/#sigma_{MIGRAD};Pseudoexperiments', 25, -5, 5)
for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   if tree.mcmass == 172.5:
      hpull.Fill( (tree.mt - meanmt[3])/tree.mt_err )

#cpull = TCanvas('cpull','cpull',800,800)
#hpull.Draw()

#temp
devts = 49000
varmt[0] /= 280000.0/devts
varmt[1] /= 390000.0/devts
varmt[2] /= 260000.0/devts
varmt[3] /= 640000.0/devts
varmt[4] /= 290000.0/devts
varmt[5] /= 440000.0/devts
varmt[6] /= 270000.0/devts

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
gresults.GetXaxis().SetTitle('M_{t}^{MC} [GeV]')
gresults.GetYaxis().SetTitle('M_{t}^{fit} - M_{t}^{MC} [GeV]')

fline = TF1('fline','[0]*x',150,200)
fline.SetParameter(0,0.0)
fline.SetLineStyle(7)
fline.SetLineColor(1)

canvas = TCanvas('canvas','canvas',800,800)
canvas.cd()

gresults.SetMaximum(1.5)
gresults.SetMinimum(-1.5)
gresults.Draw('AEP')
fline.Draw('same')

fcalib = TF1('fcalib','pol1(0)',150,200)
#fcalib = TF1('fcalib','[0] + [1]*(x-172.5)',150,200)
gresults.Fit(fcalib,"N")
fcalib.SetLineWidth(2)
fcalib.Draw('same')
gresults.Draw('EP')

#latex = TLatex()
#latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.04)
latex.DrawLatex(0.861, 0.887, fit_type)
latex.SetTextSize(0.03)
latex.DrawLatex(0.687, 0.842, 'slope = '+str(round(fcalib.GetParameter(1),3))+' #pm '+str(round(fcalib.GetParError(1),3)))
latex.DrawLatex(0.676, 0.8, 'y-intercept = '+str(round(fcalib.GetParameter(0),2))+' #pm '+str(round(fcalib.GetParError(0),2)))
#latex.DrawLatex(0.687, 0.894, 'slope = '+str(round(fcalib.GetParameter(1),3))+' #pm '+str(round(fcalib.GetParError(1),3)))
#latex.DrawLatex(0.676, 0.852, 'y-intercept = '+str(round(fcalib.GetParameter(0),1))+' #pm '+str(round(fcalib.GetParError(0),1)))

CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)

canvas.Print('calib_mt.pdf')

raw_input('')
