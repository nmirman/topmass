#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TAxis.h"

#include <stdlib.h>
#include <iostream>

using namespace std;

void slopetest(){

   string dist = "maos210";
   double nevts = 2*50000.0;
   string dim = "mt";

   double p1=0, p2=0, p3=0;
   double dx = 0;
   double eps = 0.2;
   if( dim == "mt" ){
      p1 = 172.0; p2 = 172.5, p3 = 173.0;
      //p1 = 172.5; p2 = 173, p3 = 173.5;
      //p1 = 166.0; p2 = 166.5, p3 = 167.0;
      //p1 = 178.0; p2 = 178.5, p3 = 179.0;
   }
   if( dim == "jes" ){
      p1 = 0.999; p2 = 1.000, p3 = 1.001;
   }
   if( dist == "mbl" ){
      dx = 2.8;
   }
   if( dist == "mt2_221" ){
      dx = 1.9;
   }

   TFile *f = new TFile("../results/plotsTemplates.root");
   TDirectory *d = (TDirectory*)f->Get((dim+"shape_"+dist+"_gp").c_str());

   /*
   TCanvas *ctemp = (TCanvas*)f->Get(("c_"+dist+"_gp_signal_JSF1000").c_str());
   TH1 *h = (TH1*)(ctemp->GetListOfPrimitives()->At(1));
   dx = h->GetBinWidth(1);

   delete ctemp;
   */
   dx = 10;

   TGraph *gslope = new TGraph();
   TGraph *gint = new TGraph();

   double sig = 0;
   double d2sig = 0;

   TIter nextkey(d->GetListOfKeys());
   TKey *key;
   while( (key = (TKey*)nextkey()) ){

      string name = key->GetName();
      string classname = key->GetClassName();

      TCanvas *c = (TCanvas*)d->Get(name.c_str());
      TGraph *g = (TGraph*)(c->GetListOfPrimitives()->At(1));

      double y1=0, y2=0;
      double ycent = 0;
      for(int i=0; i < g->GetN(); i++){
         double x=0, y=0;
         g->GetPoint(i, x, y);
         if( x == p1 ) y1 = y;
         if( x == p2 ) ycent = y;
         if( x == p3 ) y2 = y;
      }
      double slope = y2-y1;

      string sbin = name;
      string ntemp = "c"+dist+"_gp";
      sbin.erase(sbin.begin(), sbin.begin()+ntemp.length());
      if( dist == "mt2_221" ) sbin.erase(sbin.end(), sbin.end()+1);
      double dbin = atof(sbin.c_str());

      double integrand = slope*slope/ycent;

      gint->SetPoint(gint->GetN(), dbin, integrand);
      gslope->SetPoint(gslope->GetN(), dbin, slope);

      sig += slope*slope*eps*eps*-0.5*dx/ycent;
      d2sig += integrand*dx;

   }

   cout << "Projected S(" << eps << ") = " << sig-2.0*nevts << endl;
   cout << "Projected sigma = " << sqrt(1.0/(d2sig*nevts)) << endl;

   TCanvas *c1 = new TCanvas();
   gslope->SetLineWidth(2);
   gslope->GetXaxis()->SetTitle("(GeV)");
   gslope->GetXaxis()->SetTitleSize(0.05);
   gslope->GetXaxis()->SetTitleOffset(0.8);
   gslope->Draw("AC");

   TCanvas *c2 = new TCanvas();
   gint->SetLineWidth(2);
   gint->GetXaxis()->SetTitle("(GeV)");
   gint->GetXaxis()->SetTitleSize(0.05);
   gint->GetXaxis()->SetTitleOffset(0.8);
   gint->Draw("AC");

   return;
}
