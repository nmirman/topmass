#! /usr/bin/env python

from ROOT import *
from array import array
import sys
import os

inname1 = sys.argv[1]
inname2 = sys.argv[2]
#nameprefix = os.path.splitext(inname)[0]

nameprefix = '/uscms_data/d3/shtan/templateplots/oct20code/corr/plotcorr_'

file1 = TFile( inname1 )
tree1 = file1.Get('FitResults')
file2 = TFile( inname2 )
tree2 = file2.Get('FitResults')

#gresults = TGraphErrors()
gcorr = TGraph ()
overallcorr = array( 'f', [0] )
corrtree = TTree ('corrtree','corrtree')
corrtree.Branch('overallcorr', overallcorr, 'overallcorr/F')

#chi2=0
#summ=0
for i in range( tree1.GetEntries() ):
   tree1.GetEntry(i)
   mt1 = tree1.mt
   tree2.GetEntry(i)
   mt2 = tree2.mt
   gcorr.SetPoint(i, mt1, mt2);
   

#   summ += tree.mt-tree.mcmass
#   gresults.SetPoint(i, tree.mcmass, tree.mt-tree.mcmass)
#   gresults.SetPointError(i, 0.0, tree.mt_err)
#   gchi2.SetPoint(i, tree.mcmass, tree.fitchi2)
#   chi = (tree.mt-tree.mcmass)/tree.mt_err
#   chi2 += chi*chi
#   print str(tree.mcmass)+': '+str(tree.mt)+' +- '+str(tree.mt_err)


#chi2tree.GetEntry(0)
#chi2tree.overallchi2 = chi2
#overallchi2[0] = chi2
#chi2tree.Fill()
#ave=summ/8


overallcorr = gcorr.GetCorrelationFactor()
print 'correlation coefficient = '+str(overallcorr)

title = sys.argv[5]
xlabel = sys.argv[3]
ylabel = sys.argv[4]
gcorr.SetMarkerStyle(20)
gcorr.SetTitle(title)
gcorr.GetXaxis().SetTitle('Top Mass from '+xlabel+' (GeV)')
gcorr.GetYaxis().SetTitle('Top Mass from '+ylabel+' (GeV)')

#fline = TF1('fline','[0]*x',150,200)
#fline.SetParameter(0,0.0)
#fline.SetLineStyle(7)

#aveline = TF1('fline','[0]*x+[1]',150,200)
#aveline.SetParameter(0,0.0)
#aveline.SetParameter(1,ave)
#aveline.SetLineStyle(7)
#aveline.SetLineColor(4)

canvas = TCanvas('resultscanvas','resultscanvas',800,800)
canvas.cd()

#gresults.SetMaximum(1.5)
#gresults.SetMinimum(-1.5)
gcorr.Draw('AEP')
#fline.Draw('same')
#aveline.Draw('same')


#canvas.Print('results/'+sys.argv[1]+'_bigplot.root')

#chi2canvas = TCanvas('chi2canvas','chi2canvas',800,800)
#chi2canvas.cd()

#gchi2.Draw('AP')

#chi2canvas.Print('results/'+sys.argv[1]+'_bigplot.root')

fileout = TFile(nameprefix+'_xlabel'+'_'+'ylabel'+'.root','RECREATE')
fileout.cd()
canvas.Write()
#chi2canvas.Write()
corrtree.Write()
fileout.Close()

raw_input('')
