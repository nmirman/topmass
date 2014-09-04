#include "Bigplotting.h"
//#include "TopMass.h"
//#include "Mt2Calculator.h"

#include "TNtuple.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

Bigplotting::Bigplotting(){
}

Bigplotting::~Bigplotting(){
}

void Bigplotting::BigPlot( fitresultsfile ){

	double mt, mt_err, fitchi2, mcmass;

	TFile filein ( fitresultsfile);
	TTree *treein = (TTree*)filein.Get("FitResults");
	treein->SetBranchAddress("mt",&mt);
	treein->SetBranchAddress("mt_err",&mt_err);
	treein->SetBranchAddress("fitchi2",&fitchi2);
	treein->SetBranchAddress("mcmass",&mcmass);

	vector<double> mtsdif, mt_errs, fitchi2s, mcmasses, yerrs;
	mtsdif.clear();
	mt_errs.clear();
	fitchi2s.clear();
	mcmasses.clear();
	yerrs.clear();

	for (int i=0; i<treein->GetEntries(); i++){
		treein->GetEntry(i);
		mtsdif.push_back(mt-mcmass);
		mt_errs.push_back(mt_err);
		fitchi2s.push_back(fitchi2);
		mcmasses.push_back(mcmass);
		yerrs.push_back(0);
	}

	TGraphErrors *plotbig = new TGraphErrors( mcmasses.size(), &(mcmasses[0]), &(mtsdif[0]), &(yerrs[0]), &(mt_errs[0]) );
	plotbig->SetMarkerColor(1);
	plotbig->SetMarkerStyle(7);
	plotbig->SetTitle("Maos from 220, 20GeV, eta1 only");
	plotbig->GetXaxis()->SetTitle("MC Mass (GeV)");
	plotbig->GetYaxis()->SetTitle("Fitted - MC Mass (GeV)"); 
//		graphMetdif->GetXaxis()->SetTitleSize(100.0);
//		graphMetdif->GetYaxis()->SetTitleSize(100.0);

	TLine *line = new TLine(0,0,200,0);
	line->SetLineColor(2);
	line->SetLineWidth(1);
	line->SetLineStyle(2);

	TGraph *plotchi = new TGraph( mcmasses.size(), &(mcmasses[0]), &(fitchi2s[0]) );
	plotchi->SetMarkerColor(1);
	plotchi->SetMarkerStyle(2);
	plotchi->SetTitle("Maos from 220, 20GeV, eta1 only");
	plotchi->GetXaxis()->SetTitle("MC Mass (GeV)");
	plotchi->GetYaxis()->SetTitle("Fit Chi^2 (GeV)"); 
	
	TCanvas *canvasPlotbig = new TCanvas ( "plotbig", "plotbig", 700, 700);
	TCanvas *canvasPlotchi = new TCanvas ( "plotchi", "plotchi", 700, 700);

	TFile *fileout = new TFile ( (fitresultsfile+"_plotbig.root").c_str(), "RECREATE");
	fileout->cd();

	canvasPlotbig->cd();
	plotbig->Draw("AP");	
	line->Draw("SAME");
	canvasPlotbig->Write();
	canvasPlotbig->ls();

	canvasPlotchi->cd();
	plotchi->Draw("AP");
	canvasPlotchi->Write();
	canvasPlotchi->ls();

	fileout->Close();

	delete canvasPlotbig;
	delete plotbig;
	delete canvasPlotchi;
	delete plotchi;
	delete line;

}


