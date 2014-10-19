#! /usr/bin/env python

from ROOT import *
import sys

file = TFile( sys.argv[1] )
tree = file.Get('FitResults')

gresults = TGraphErrors()
chi2=0
for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   gresults.SetPoint(i, tree.mcmass, tree.mt-tree.mcmass)
   gresults.SetPointError(i, 0.0, tree.mt_err)
   chi = (tree.mt-tree.mcmass)/tree.mt_err
   chi2 += chi*chi
   print str(tree.mcmass)+': '+str(tree.mt)+' +- '+str(tree.mt_err)

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
