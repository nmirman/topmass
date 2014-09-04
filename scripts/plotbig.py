#! /usr/bin/env python

from ROOT import *
from array import array
import sys
import os

inname = sys.argv[1]
nameprefix = os.path.splitext(inname)[0]

file = TFile( inname )
tree = file.Get('FitResults')

gresults = TGraphErrors()
gchi2 = TGraph ()
overallchi2 = array( 'f', [0] )
chi2tree = TTree ('chi2tree','chi2tree')
chi2tree.Branch('overallchi2', overallchi2, 'overallchi2/F')

chi2=0
for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   gresults.SetPoint(i, tree.mcmass, tree.mt-tree.mcmass)
   gresults.SetPointError(i, 0.0, tree.mt_err)
   gchi2.SetPoint(i, tree.mcmass, tree.fitchi2)
   chi = (tree.mt-tree.mcmass)/tree.mt_err
   chi2 += chi*chi
   print str(tree.mcmass)+': '+str(tree.mt)+' +- '+str(tree.mt_err)

#chi2tree.GetEntry(0)
#chi2tree.overallchi2 = chi2
overallchi2[0] = chi2
chi2tree.Fill()

print 'chi2 = '+str(chi2)

title = sys.argv[2]
gresults.SetMarkerStyle(20)
gresults.SetTitle(title)
gresults.GetXaxis().SetTitle('MC Top Mass (GeV)')
gresults.GetYaxis().SetTitle('Fitted - MC Top Mass (GeV)')
gchi2.SetMarkerStyle(20)
gchi2.SetTitle(title)
gchi2.GetXaxis().SetTitle('MC Top Mass (GeV)')
gchi2.GetYaxis().SetTitle('Fit chi^{2}')


fline = TF1('fline','[0]*x',150,200)
fline.SetParameter(0,0.0)
fline.SetLineStyle(7)

canvas = TCanvas('resultscanvas','resultscanvas',800,800)
canvas.cd()

#gresults.SetMaximum(1.5)
#gresults.SetMinimum(-1.5)
gresults.Draw('AEP')
fline.Draw('same')

#canvas.Print('results/'+sys.argv[1]+'_bigplot.root')

chi2canvas = TCanvas('chi2canvas','chi2canvas',800,800)
chi2canvas.cd()

gchi2.Draw('AP')

#chi2canvas.Print('results/'+sys.argv[1]+'_bigplot.root')

fileout = TFile(nameprefix+'_bigplot.root','RECREATE')
fileout.cd()
canvas.Write()
chi2canvas.Write()
chi2tree.Write()
fileout.Close()

raw_input('')
