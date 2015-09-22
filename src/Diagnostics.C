#include "TopMass.h"
#include "Shapes.h"
#include "Mt2Calculator.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TRandom3.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;


void Fitter::DeclareHists( map< string, map<string, TH1D*> >& hists_, map< string, map<string, TH2D*> >& hists2d_, string label ){

   vector<string> type;
   type.push_back("data");
   type.push_back("other");
   for(int i=0; i < NMP; i++){
      stringstream dstr;
      dstr << floor(masspoints[i]);
      string dname = "ttbar"+dstr.str();
      type.push_back(dname+"_signal");
      type.push_back(dname+"_mistag");
      type.push_back(dname+"_taus");
      type.push_back(dname+"_hadronic");
   }
   type.push_back("fitevts"); // for displaying fit results (could be mc or data)

   for(vector<string>::iterator t = type.begin(); t < type.end(); t++){

      string name = *t;
      string namel = *t+label;

      //
      // variables
      //
      hists_["mbl_gp"][name] = new TH1D( ("hmbl_gp_"+namel).c_str(),
            "M_{bl};M_{bl} (GeV);Events", 100, dists["mbl_gp"].lbnd, dists["mbl_gp"].rbnd );
      hists_["mt2_220_gp"][name] = new TH1D( ("hmt2_220_gp_"+namel).c_str(),
            "M_{T2} 220;M_{T2} 220 (GeV);Events", 100, dists["mt2_220_gp"].lbnd, dists["mt2_220_gp"].rbnd );
      hists_["mt2_221_gp"][name] = new TH1D( ("hmt2_221_gp_"+namel).c_str(),
            "M_{T2} 221;M_{T2} 221 (GeV);Events", 100, dists["mt2_221_gp"].lbnd, dists["mt2_221_gp"].rbnd );
      hists_["maos220_gp"][name] = new TH1D( ("hmaos220_gp_"+namel).c_str(),
            "MAOS from M_{T2} 220;blv mass(GeV);Events", 100, dists["maos220_gp"].lbnd, dists["maos220_gp"].rbnd );
      hists_["maos210_gp"][name] = new TH1D( ("hmaos210_gp_"+namel).c_str(),
            "MAOS from M_{T2} 210;blv mass(GeV);Events", 100, dists["maos210_gp"].lbnd, dists["maos210_gp"].rbnd );

      hists_["mt2_220"][name] = new TH1D( ("hmt2_220_"+namel).c_str(),
            "M_{T2} 220;M_{T2} 220 (GeV);Events/3 GeV", 100, 0, dists["mt2_220_nomatchmbl"].rbnd );
      hists_["mt2_221"][name] = new TH1D( ("hmt2_221_"+namel).c_str(),
            "M_{T2} 221;M_{T2} 221 (GeV);Events/2 GeV", 100, 50, 250 );
      hists_["mbl"][name] = new TH1D( ("hmbl_"+namel).c_str(),
            "M_{bl};M_{bl} (GeV);Events/3 GeV", 100, 0, dists["mbl"].rbnd );
      hists_["mt2_210"][name] = new TH1D( ("hmt2_210_"+namel).c_str(),
            "M_{T2} 210;M_{T2} 210 (GeV);Events/1.5 GeV", 100, 0, 150 );
      hists_["mt2_220_matchmbl"][name] = new TH1D( ("hmt2_220_matchmbl_"+namel).c_str(),
            "M_{T2} 220;M_{T2} 220 (GeV);Events/2.5 GeV", 100, 0, dists["mt2_220_nomatchmbl"].rbnd );
      hists_["mt2_220_nomatchmbl"][name] = new TH1D( ("hmt2_220_nomatchmbl_"+namel).c_str(),
            "M_{T2} 220;M_{T2} 220 (GeV);Events/3 GeV", 100, 0, dists["mt2_220_nomatchmbl"].rbnd );
      hists_["maos220blv"][name] = new TH1D( ("hmaos220blv_"+namel).c_str(),
            "MAOS from M_{T2} 220;blv mass(GeV);Events/5GeV", 100, 0, dists["maos220blv"].rbnd );
      hists_["maos210blv"][name] = new TH1D( ("hmaos210blv_"+namel).c_str(),
            "MAOS from M_{T2} 210;blv mass(GeV);Events/5GeV", 100, 0, dists["maos210blv"].rbnd );


      //
      // kinematic distributions
      //
      hists2d_["220Vmbl"][name] = new TH2D( ("h220Vmbl_"+namel).c_str(),
            "MT2 220 vs. mbl;M_{bl} (GeV);M_{T2} (GeV) ", 50, 0, dists["mbl"].rbnd, 50, 0, dists["mt2_220_nomatchmbl"].rbnd );
      hists2d_["maos220blvVmbl"][name] = new TH2D( ("hmaos220blvVmbl_"+namel).c_str(),
            "MAOS 220 bl#nu vs. mbl;M_{bl} (GeV);MAOS bl#nu (GeV)", 50, 0, dists["mbl"].rbnd, 50, 0, dists["maos220blv"].rbnd );
      hists2d_["maos210blvVmbl"][name] = new TH2D( ("hmaos210blvVmbl_"+namel).c_str(),
            "MAOS 210 bl#nu vs. mbl;M_{bl} (GeV);MAOS bl#nu (GeV)", 50, 0, dists["mbl"].rbnd, 50, 0, dists["maos210blv"].rbnd );
      hists2d_["maos220blvV220"][name] = new TH2D( ("hmaos220blvV220_"+namel).c_str(),
            "MAOS 220 bl#nu vs. M_{T2} 220;M_{T2} (GeV);MAOS bl#nu (GeV)", 50, 0, dists["mt2_220_nomatchmbl"].rbnd,
            50, 0, dists["maos220blv"].rbnd );
      hists2d_["maos210blvV220"][name] = new TH2D( ("hmaos210blvV220_"+namel).c_str(),
            "MAOS 210 bl#nu vs. M_{T2} 220;M_{T2} (GeV);MAOS bl#nu (GeV)", 50, 0, dists["mt2_220_nomatchmbl"].rbnd,
            50, 0, dists["maos210blv"].rbnd );
      hists2d_["maos220blvVmaos210blv"][name] = new TH2D( ("hmaos220blvVmaos210blv_"+namel).c_str(),
            "MAOS 220 bl#nu vs. MAOS 210 bl#nu;MAOS bl#nu (GeV);MAOS bl#nu (GeV)", 50, 0, dists["maos210blv"].rbnd,
            50, 0, dists["maos220blv"].rbnd );

      // pt, eta, phi
      hists_["b_pt"][name] = new TH1D( ("hb_pt_"+namel).c_str(),
            "b quark p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["b_eta"][name] = new TH1D( ("hb_eta_"+namel).c_str(),
            "b quark #eta;#eta;Events", 100, -3, 3 );
      hists_["b_phi"][name] = new TH1D( ("hb_phi_"+namel).c_str(),
            "b quark #phi;#phi;Events", 100, -3.5, 3.5 );
      hists_["l_pt"][name] = new TH1D( ("hl_pt_"+namel).c_str(),
            "lepton p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["l_eta"][name] = new TH1D( ("hl_eta_"+namel).c_str(),
            "lepton #eta;#eta;Events", 100, -3, 3 );
      hists_["l_phi"][name] = new TH1D( ("hl_phi_"+namel).c_str(),
            "lepton #phi;#phi;Events", 100, -3.5, 3.5 );
      hists_["bb_pt"][name] = new TH1D( ("hbb_pt_"+namel).c_str(),
            "bb p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["bb_eta"][name] = new TH1D( ("hbb_eta_"+namel).c_str(),
            "bb Rapidity;Rapidity;Events", 100, -3, 3 );
      hists_["bb_phi"][name] = new TH1D( ("hbb_phi_"+namel).c_str(),
            "bb #phi;#phi;Events", 100, -3.5, 3.5 );
      hists_["ll_pt"][name] = new TH1D( ("hll_pt_"+namel).c_str(),
            "ll p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["ll_eta"][name] = new TH1D( ("hll_eta_"+namel).c_str(),
            "ll Rapidity;Rapidity;Events", 100, -3, 3 );
      hists_["ll_phi"][name] = new TH1D( ("hll_phi_"+namel).c_str(),
            "ll #phi;#phi;Events", 100, -3.5, 3.5 );
      hists_["bbll_pt"][name] = new TH1D( ("hbbll_pt_"+namel).c_str(),
            "bbll p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["bbll_eta"][name] = new TH1D( ("hbbll_eta_"+namel).c_str(),
            "bbll Rapidity;Rapidity;Events", 100, -3, 3 );
      hists_["bbll_phi"][name] = new TH1D( ("hbbll_phi_"+namel).c_str(),
            "bbll #phi;#phi;Events", 100, -3.5, 3.5 );

      // invariant mass
      hists_["bb_m"][name] = new TH1D( ("hbb_m_"+namel).c_str(),
            "bb invariant mass;m_{bb} (GeV);Events/4 GeV", 100, 0, 400 );
      hists_["ll_m"][name] = new TH1D( ("hll_m_"+namel).c_str(),
            "ll invariant mass;m_{ll} (GeV);Events/4 GeV", 100, 0, 400 );
      hists_["bbl_m"][name] = new TH1D( ("hbbl_m_"+namel).c_str(),
            "bbl invariant mass;m_{bbl} (GeV);Events/6 GeV", 100, 0, 600 );
      hists_["bll_m"][name] = new TH1D( ("hbll_m_"+namel).c_str(),
            "bll invariant mass;m_{bll} (GeV);Events/6 GeV", 100, 0, 600 );
      hists_["bbll_m"][name] = new TH1D( ("hbbll_m_"+namel).c_str(),
            "bbll invariant mass;m_{bbll} (GeV);Events/10 GeV", 100, 0, 1000 );

      // angles
      hists_["bb_dR"][name] = new TH1D( ("hbb_dR_"+namel).c_str(),
            "#DeltaR between bb;#DeltaR;Events", 100, 0, 5.0 );
      hists_["bl_dR"][name] = new TH1D( ("hbl_dR_"+namel).c_str(),
            "#DeltaR between bl;#DeltaR;Events", 100, 0, 5.0 );
      hists_["ll_dR"][name] = new TH1D( ("hll_dR_"+namel).c_str(),
            "#DeltaR between ll;#DeltaR;Events", 100, 0, 5.0 );

      // met
      hists_["met_x"][name] = new TH1D( ("hmet_x_"+namel).c_str(),
            "MET_{x};MET_{x} (GeV);Events/4 GeV", 100, -200, 200 );
      hists_["met_y"][name] = new TH1D( ("hmet_y_"+namel).c_str(),
            "MET_{y};MET_{y} (GeV);Events/4 GeV", 100, -200, 200 );
      hists_["met"][name] = new TH1D( ("hmet_"+namel).c_str(),
            "total MET;MET (GeV);Events/3 GeV", 100, 0, 300 );
   }

   return;
}

void Fitter::DeleteHists( map< string, map<string, TH1D*> >& hists_, map< string, map<string, TH2D*> >& hists2d_ ){

   typedef map<string, TH1D*> tmap;
   typedef map<string, tmap> hmap;

   for( hmap::iterator h = hists_.begin(); h != hists_.end(); h++){
      for( tmap::iterator t = h->second.begin(); t != h->second.end(); t++){
         delete (t->second);
      }
   }

   typedef map<string, TH2D*> t2map;
   typedef map<string, t2map> h2map;

   for( h2map::iterator h = hists2d_.begin(); h != hists2d_.end(); h++){
      for( t2map::iterator t = h->second.begin(); t != h->second.end(); t++){
         delete (t->second);
      }
   }


}

void Fitter::FillHists( map< string, map<string, TH1D*> >& hists_, map< string, map<string, TH2D*> >& hists2d_,
      vector<Event>& eventvec, bool fit_events ){

   // event loop
   for( vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){

      string type = fit_events ? "fitevts" : ev->type;

      if( fit_events and !(ev->fit_event) ){
         cout << "ERROR: FILL WITH FIT EVENTS" << endl;
         return;
      }

      // define objects
      TLorentzVector jet1 = ev->jet1;
      TLorentzVector jet2 = ev->jet2;
      TLorentzVector lep1 = ev->lep1;
      TLorentzVector lep2 = ev->lep2;
      TLorentzVector met = ev->met;
      TLorentzVector up221 = -((jet1)+(jet2)+(lep1)+(lep2)+(met));
      TLorentzVector bb = jet1 + jet2;
      TLorentzVector ll = lep1 + lep2;
      TLorentzVector bbll = jet1 + jet2 + lep1 + lep2;
      TLorentzVector bll1 = jet1 + lep1 + lep2;
      TLorentzVector bll2 = jet2 + lep1 + lep2;
      TLorentzVector bbl1 = jet1 + jet2 + lep1;
      TLorentzVector bbl2 = jet1 + jet2 + lep2;

      if (sin((jet1).DeltaPhi(up221))*sin((jet2).DeltaPhi(up221)) > 0){
         hists_["mt2_221"][type]->Fill( ev->mt2_221, ev->weight );
         hists_["mt2_221_gp"][type]->Fill( ev->mt2_221, ev->weight );
      }

      //MAOS
      double blv210array [] = { ev->maos210_blvmass1ap, ev->maos210_blvmass1am, ev->maos210_blvmass2ap, ev->maos210_blvmass2am, ev->maos210_blvmass1bp, ev->maos210_blvmass1bm, ev->maos210_blvmass2bp, ev->maos210_blvmass2bm };
      double blv220array [] = { ev->maos220_blvmass1ap, ev->maos220_blvmass1am, ev->maos220_blvmass2ap, ev->maos220_blvmass2am, ev->maos220_blvmass1bp, ev->maos220_blvmass1bm, ev->maos220_blvmass2bp, ev->maos220_blvmass2bm };

      vector<bool> useMaos220 = ev->maoscut220;
      for (int i=0; i<8; i++){
         if (useMaos220[i]){
            hists_["maos220blv"][type]->Fill( blv220array[i], ev->weight );
            hists_["maos220_gp"][type]->Fill( blv220array[i], ev->weight );
         }
      }
      
      vector<bool> useMaos210 = ev->maoscut210;
      for (int i=0; i<8; i++){
         if (useMaos210[i]){
            hists_["maos210blv"][type]->Fill( blv210array[i], ev->weight );
            hists_["maos210_gp"][type]->Fill( blv210array[i], ev->weight );
         }
      }

      // mbl
      bool matchmbl = false;
      for( unsigned int m=0; m < ev->mbls.size(); m++ ){

         if( ev->mbls[m] == ev->mt2_220 ) matchmbl = true;
         else hists2d_["220Vmbl"][type]->Fill( ev->mbls[m], ev->mt2_220, ev->weight );

         hists_["mbl"][type]->Fill( ev->mbls[m], ev->weight );
         hists_["mbl_gp"][type]->Fill( ev->mbls[m], ev->weight );
         for(int i=0; i<8; i++){
            if (useMaos220[i]){
               hists2d_["maos220blvVmbl"][type]->Fill( ev->mbls[m], blv220array[i], ev->weight );
            }
            if (useMaos210[i]){
               hists2d_["maos210blvVmbl"][type]->Fill( ev->mbls[m], blv210array[i], ev->weight );
            }
         }
      }

      // mt2
      hists_["mt2_220"][type]->Fill( ev->mt2_220, ev->weight );
      if( ev->mt2_210 > 1 ) hists_["mt2_210"][type]->Fill( ev->mt2_210, ev->weight );

      if( matchmbl){
         hists_["mt2_220_matchmbl"][type]->Fill( ev->mt2_220, ev->weight ); 
      }
      else{
         hists_["mt2_220_nomatchmbl"][type]->Fill( ev->mt2_220, ev->weight );
         hists_["mt2_220_gp"][type]->Fill( ev->mt2_220, ev->weight );
         for(int i=0; i<8; i++){
            if (useMaos220[i]){
               hists2d_["maos220blvV220"][type]->Fill( ev->mt2_220, blv220array[i], ev->weight );
               for(int j=0; j<8; j++){
                  if (useMaos210[j]){
                     hists2d_["maos220blvVmaos210blv"][type]->Fill( blv210array[j], blv220array[i], ev->weight );
                  }
               }
            }
            if (useMaos210[i]){
               hists2d_["maos210blvV220"][type]->Fill( ev->mt2_220, blv210array[i], ev->weight );
            }
         }
      }

      //
      // kinematic distributions
      //

      double bb_dR = fabs(jet1.DeltaR(jet2));
      double ll_dR = fabs(lep1.DeltaR(lep2));
      double bl_dR1 = fabs(jet1.DeltaR(lep1));
      double bl_dR2 = fabs(jet2.DeltaR(lep1));
      double bl_dR3 = fabs(jet1.DeltaR(lep2));
      double bl_dR4 = fabs(jet2.DeltaR(lep2));

      // pt, eta, phi
      hists_["b_pt"][type]->Fill( jet1.Pt(), ev->weight );
      hists_["b_pt"][type]->Fill( jet2.Pt(), ev->weight );
      hists_["b_eta"][type]->Fill( jet1.Eta(), ev->weight );
      hists_["b_eta"][type]->Fill( jet2.Eta(), ev->weight );
      hists_["b_phi"][type]->Fill( jet1.Phi(), ev->weight );
      hists_["b_phi"][type]->Fill( jet2.Phi(), ev->weight );
      hists_["l_pt"][type]->Fill( lep1.Pt(), ev->weight );
      hists_["l_pt"][type]->Fill( lep2.Pt(), ev->weight );
      hists_["l_eta"][type]->Fill( lep1.Eta(), ev->weight );
      hists_["l_eta"][type]->Fill( lep2.Eta(), ev->weight );
      hists_["l_phi"][type]->Fill( lep1.Phi(), ev->weight );
      hists_["l_phi"][type]->Fill( lep2.Phi(), ev->weight );
      hists_["bb_pt"][type]->Fill( bb.Pt(), ev->weight );
      hists_["bb_eta"][type]->Fill( bb.Rapidity(), ev->weight );
      hists_["bb_phi"][type]->Fill( bb.Phi(), ev->weight );
      hists_["ll_pt"][type]->Fill( ll.Pt(), ev->weight );
      hists_["ll_eta"][type]->Fill( ll.Rapidity(), ev->weight );
      hists_["ll_phi"][type]->Fill( ll.Phi(), ev->weight );
      hists_["bbll_pt"][type]->Fill( bbll.Pt(), ev->weight );
      hists_["bbll_eta"][type]->Fill( bbll.Rapidity(), ev->weight );
      hists_["bbll_phi"][type]->Fill( bbll.Phi(), ev->weight );

      // invariant mass
      hists_["bb_m"][type]->Fill( bb.M(), ev->weight );
      hists_["ll_m"][type]->Fill( ll.M(), ev->weight );
      hists_["bbll_m"][type]->Fill( bbll.M(), ev->weight );
      hists_["bll_m"][type]->Fill( bll1.M(), ev->weight );
      hists_["bll_m"][type]->Fill( bll2.M(), ev->weight );
      hists_["bbl_m"][type]->Fill( bbl1.M(), ev->weight );
      hists_["bbl_m"][type]->Fill( bbl2.M(), ev->weight );

      // dR
      hists_["bb_dR"][type]->Fill( bb_dR, ev->weight );
      hists_["ll_dR"][type]->Fill( ll_dR, ev->weight );
      hists_["bl_dR"][type]->Fill( bl_dR1, ev->weight );
      hists_["bl_dR"][type]->Fill( bl_dR2, ev->weight );
      hists_["bl_dR"][type]->Fill( bl_dR3, ev->weight );
      hists_["bl_dR"][type]->Fill( bl_dR4, ev->weight );

      // met
      hists_["met"][type]->Fill( met.Pt(), ev->weight );
      hists_["met_x"][type]->Fill( met.Px(), ev->weight );
      hists_["met_y"][type]->Fill( met.Py(), ev->weight );

   }

   return;
}

void Fitter::PrintHists( map< string, map<string, TH1D*> >& hists_, map< string, map<string, TH2D*> >& hists2d_, string outfile ){
   cout << "Printing data/mc plots." << endl;

   typedef map<string, TH1D*> tmap;
   typedef map<string, tmap> hmap;

   typedef map<string, TH2D*> t2map;
   typedef map<string, t2map> h2map;

   TFile *fileout = new TFile( outfile.c_str(), "RECREATE" );
   fileout->cd();

   TDirectory *dall = fileout->mkdir( "all" );
   TDirectory *dkin = fileout->mkdir( "kinematics" );
   TDirectory *dmass = fileout->mkdir( "masses" );
   TDirectory *d2dim = fileout->mkdir( "2dim" );

   // print all hists, all channels
   dall->cd();
   for( hmap::iterator h = hists_.begin(); h != hists_.end(); h++){
      TDirectory *dtemp = dall->mkdir( (h->first).c_str() );
      dtemp->cd();
      for( tmap::iterator t = h->second.begin(); t != h->second.end(); t++){
         t->second->Write();
      }
   }

   // print all hists, all channels
   d2dim->cd();
   for( h2map::iterator h = hists2d_.begin(); h != hists2d_.end(); h++){
      TDirectory *dtemp = d2dim->mkdir( (h->first).c_str() );
      dtemp->cd();
      for( t2map::iterator t = h->second.begin(); t != h->second.end(); t++){
         t->second->Write();
      }
   }
      

   // data/mc stack plots w/ mt = 172.5
   for(hmap::iterator h = hists_.begin(); h != hists_.end(); h++){

      string name = h->first;
      tmap hist = h->second;
      TH1D* hdata = (TH1D*)hist["data"]->Clone("data_cl");

      // formatting
      TH1D* httbar_signal = (TH1D*)hist["ttbar172_signal"]->Clone("ttbar172_signal_cl");
      TH1D* httbar_mistag = (TH1D*)hist["ttbar172_mistag"]->Clone("ttbar172_mistag_cl");
      TH1D* httbar_taus = (TH1D*)hist["ttbar172_taus"]->Clone("ttbar172_taus_cl");
      TH1D* httbar_hadronic = (TH1D*)hist["ttbar172_hadronic"]->Clone("ttbar172_hadronic_cl");
      TH1D* hother = (TH1D*)hist["other"]->Clone("other_cl");

      httbar_signal->SetFillColor( kAzure-9 );
      httbar_mistag->SetFillColor( kYellow-9 );
      httbar_taus->SetFillColor( kRed-10 );
      httbar_hadronic->SetFillColor( kCyan-10 );
      hother->SetFillColor( kGreen-10 );

      TH1D* hallmc = (TH1D*)httbar_signal->Clone("hallmc");
      hallmc->Add( httbar_mistag );
      hallmc->Add( httbar_taus );
      hallmc->Add( httbar_hadronic );
      hallmc->Add( hother );

      double norm = 1.0;//hdata->Integral("width") / hallmc->Integral("width");

      // renormalize distribution in mc for all processes
      httbar_signal->Scale( norm );
      httbar_mistag->Scale( norm );
      httbar_taus->Scale( norm );
      httbar_hadronic->Scale( norm );
      hother->Scale( norm );
      hallmc->Scale( norm );

      TCanvas * canvas = new TCanvas( ("c"+name).c_str(), ("c"+name).c_str(), 800, 800 );
      canvas->SetFillColor(0);
      TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
      pad1->SetTopMargin(0.1);
      pad1->SetBottomMargin(0.01);
      pad1->SetRightMargin(0.1);
      pad1->SetFillColor(0);
      pad2->SetTopMargin(0.01);
      pad2->SetBottomMargin(0.3);
      pad2->SetRightMargin(0.1);
      pad2->SetFillColor(0);
      pad1->Draw();
      pad2->Draw();

      pad1->cd();

      THStack * hstack = new THStack( ("hs"+name).c_str(), ("hs"+name).c_str() );
      hstack->Add( hother );
      hstack->Add( httbar_hadronic );
      hstack->Add( httbar_taus );
      hstack->Add( httbar_mistag );
      hstack->Add( httbar_signal );

      hstack->Draw("HIST");

      hstack->GetYaxis()->SetTitle( hdata->GetYaxis()->GetTitle() );
      //hstack->SetMinimum( 0.5 );
      hstack->GetXaxis()->SetTitleSize(0.00);
      hstack->GetYaxis()->SetLabelSize(0.07);
      hstack->GetYaxis()->SetTitleSize(0.08);
      hstack->GetYaxis()->SetTitleOffset(1.0);
      hstack->GetXaxis()->SetLabelFont(42);
      hstack->GetYaxis()->SetLabelFont(42);
      hstack->GetXaxis()->SetTitleFont(42);
      hstack->GetYaxis()->SetTitleFont(42);

      //hdata->SetMarkerStyle(20);
      //hdata->Draw( "same EP" );

      TLegend * legend = new TLegend(0.717,0.650,0.874,0.870);
      //legend->AddEntry( hdata, "data" );
      legend->AddEntry( httbar_signal, "signal", "f" );
      legend->AddEntry( httbar_mistag, "mistag bkg", "f" );
      legend->AddEntry( httbar_taus, "tau decays", "f" );
      legend->AddEntry( httbar_hadronic, "hadronic decays", "f" );
      legend->AddEntry( hother, "non-ttbar bkg", "f" );

      legend->SetFillStyle(0);
      legend->SetBorderSize(0);
      legend->Draw( "same" );

      pad2->cd();

      TH1D *hratio = (TH1D*)hdata->Clone("hratio");
      hratio->Divide( hallmc );
      hratio->SetTitle(";"+TString(hdata->GetXaxis()->GetTitle())+";data/mc");
      hratio->GetYaxis()->CenterTitle();
      hratio->SetStats(0);

      hratio->GetXaxis()->SetTitleSize(0.14);
      hratio->GetXaxis()->SetLabelSize(0.14);
      hratio->GetYaxis()->SetLabelSize(0.11);
      hratio->GetYaxis()->SetTitleSize(0.14);
      hratio->GetYaxis()->SetTitleOffset(0.28);
      hratio->GetXaxis()->SetLabelFont(42);
      hratio->GetYaxis()->SetLabelFont(42);
      hratio->GetXaxis()->SetTitleFont(42);
      hratio->GetYaxis()->SetTitleFont(42);
      hratio->SetMaximum( 1.6 );
      hratio->SetMinimum( 0.4 );
      hratio->GetYaxis()->SetNdivisions(505);
      hratio->Draw("EP");

      TF1 *func = new TF1("func","[0]",-10E6,10E6);
      func->SetParameter(0,1.0);
      func->SetLineWidth(1);
      func->SetLineStyle(7);
      func->SetLineColor(1);
      func->Draw("same");

      fileout->cd();
      if( name.find("220")==string::npos and name.find("221")==string::npos
            and name.find("mbl")==string::npos and name.find("210")==string::npos  ){
         dkin->cd();
      }
      canvas->Write();

      pad1->cd();
      TLatex *text = new TLatex();
      text->SetTextSize(0.05);
      text->SetNDC();
      text->DrawLatex(0.4, 0.9, hdata->GetTitle() );
      canvas->Print( ("results/plots/stack_"+name+".pdf").c_str() );

      delete canvas;
      delete legend;
      delete func;
      delete text;
   }

   TCanvas *c220_mblstack = new TCanvas("c220_mblstack","c220_mblstack",700,700);
   c220_mblstack->cd();

   TH1D* h220_matchmbl = (TH1D*)hists_["mt2_220_matchmbl"]["data"]->Clone("mt2_220_matchmbl_cl");
   TH1D* h220_nomatchmbl = (TH1D*)hists_["mt2_220_nomatchmbl"]["data"]->Clone("mt2_220_nomatchmbl_cl");

   h220_matchmbl->Add( h220_nomatchmbl);

   h220_matchmbl->SetFillColor( kRed-9 );
   h220_nomatchmbl->SetFillColor( kAzure-9 );

   h220_matchmbl->Draw( "HIST" );
   h220_nomatchmbl->Draw( "same HIST" );

   TLegend *l220_mblstack = new TLegend(0.605528,0.655866,0.866834,0.816333);
   l220_mblstack->AddEntry( h220_matchmbl, "220 = mbl", "f" );
   l220_mblstack->AddEntry( h220_nomatchmbl, "220 #neq mbl", "f" );
   l220_mblstack->SetFillStyle(0);
   l220_mblstack->SetBorderSize(0);
   l220_mblstack->Draw("same");

   c220_mblstack->Write();

   //
   // chi2 plots
   //
   gErrorIgnoreLevel = 2000;
   for(hmap::iterator h = hists_.begin(); h != hists_.end(); h++){

      string hname = h->first;
      tmap hist = h->second;
      TH1D* hdata = (TH1D*)hist["data"]->Clone("data_cl");

      TGraph *gchi2 = new TGraph();

      // do this for all top masses
      for(int i=0; i < NMP; i++){

         stringstream dstr;
         dstr << floor(masspoints[i]);
         string dname = "ttbar"+dstr.str();

         TH1D *hallmc = (TH1D*)hist[dname+"_signal"]->Clone("hallmc");
         hallmc->Add( hist[dname+"_mistag"] );
         hallmc->Add( hist[dname+"_taus"] );
         hallmc->Add( hist[dname+"_hadronic"] );
         hallmc->Add( hist["other"] );

         double chi2 = 0.0;//hdata->Chi2Test(hallmc,"UW CHI2");

         gchi2->SetPoint(i, masspoints[i], chi2);

         // write out mass hists
         if( hname.find("mbl") != string::npos or hname.find("mt2") != string::npos ){
            dmass->cd();
            hallmc->SetName( (hname+"_"+dname).c_str() );
            hallmc->Write();
            fileout->cd();
         }

      }

      TCanvas * canvas = new TCanvas( ("g"+hname).c_str(), ("g"+hname).c_str(), 800, 800 );
      canvas->SetFillColor(0);
      canvas->cd();

      gchi2->SetName(("g"+hname).c_str());
      gchi2->SetMarkerStyle(20);
      gchi2->SetTitle(hname.c_str());
      gchi2->GetXaxis()->SetTitle("Top Mass (GeV)");
      gchi2->GetYaxis()->SetTitle("#chi^{2}");

      bool yrange = false;
      string names [] = {"b_pt","bb_pt","bbl_m", "bbll_m", "mbl", 
         "mt2_220","mt2_221","mt2_220_nomatchmbl"};
      for(unsigned int i=0; i < sizeof(names)/sizeof(names[0]); i++){
         if( hname.compare(names[i]) == 0 ) yrange = true;
      }

      if( yrange ){
         gchi2->SetMaximum(4000);
         gchi2->SetMinimum(0);
      }

      gchi2->Draw("AP");

      TF1 *func = new TF1("func","pol2",160,182);
      gchi2->Fit(func,"Q");

      gStyle->SetOptFit();

      //double pb = func->GetParameter(1);
      //double pa = func->GetParameter(2);
      double minimum = 0.0;//-pb/(2*pa);

      //cout << hname << " minimum = " << minimum << endl;

      stringstream minstr;
      minstr << minimum;
      TLatex *text = new TLatex();
      text->SetTextSize(0.04);
      text->SetNDC();
      text->DrawLatex(0.5, 0.7, ("min @ "+minstr.str()+" GeV").c_str() );

      fileout->cd();
      if( hname.find("220")==string::npos and hname.find("221")==string::npos
            and hname.find("mbl")==string::npos ){
         dkin->cd();
      }
      canvas->Write();

      text->SetTextSize(0.05);
      text->DrawLatex(0.4, 0.9, hdata->GetTitle() );
      canvas->Print( ("results/plots/chi2_"+hname+".pdf").c_str() );


      delete gchi2;
      delete func;
      delete canvas;
      delete text;
   }
   fileout->Close();

   return;
}

vector<bool> Fitter::MaosCut220( vector<Event>::iterator ev, TLorentzVector &nu1ap, TLorentzVector &nu1am, TLorentzVector &nu2ap, TLorentzVector &nu2am, TLorentzVector &nu1bp, TLorentzVector &nu1bm, TLorentzVector &nu2bp, TLorentzVector &nu2bm ){
 
   bool distcut = 0;
   bool etadisamb = 0;
   bool blmatch = 0;
   double endpt220 = 172.5;

   if (maoscuts220 == 1){ distcut = 1; }
   else if (maoscuts220 == 2){ etadisamb = 1; }
   else if (maoscuts220 == 3){distcut = 1; etadisamb = 1; }
   else if (maoscuts220 == 4){blmatch = 1; }
   else if (maoscuts220 == 5){distcut = 1; blmatch = 1; }
   else if (maoscuts220 == 6){etadisamb = 1; blmatch = 1; }
   else if (maoscuts220 == 7){distcut = 1; etadisamb = 1; blmatch = 1; }

   vector<bool> toreturn;
   toreturn.clear();
   for (int i=0; i<8; i++){
      toreturn.push_back(true);
   }

   TLorentzVector neuarray [] = { nu1ap, nu1am, nu2ap, nu2am, nu1bp, nu1bm, nu2bp, nu2bm };

   //Nans
   for (int i=0; i<8; i++){

      if ( isnan( (neuarray[i]).Pz() ) ){ toreturn[i] = 0; }
   }

   //Near Endpoint
   if (distcut){
      if ( (ev->mt2_220grida < endpt220-20) or (ev->mt2_220grida > endpt220+20) ){
         for (int i=0; i<4; i++){ toreturn[i] = 0; }
      }
      if ( (ev->mt2_220gridb < endpt220-20) or (ev->mt2_220gridb > endpt220+20) ){
         for (int i=4; i<8; i++){ toreturn[i] = 0; }
      }
   }

   //Eta 1
   if (etadisamb){
      for (int i=0; i<4; i++){
         if ( isnan( (neuarray[2*i]).Pz() ) or isnan( (neuarray[2*i+1]).Pz() ) ) continue;
         if ( abs( (neuarray[2*i]).Eta() ) < abs( (neuarray[2*i+1]).Eta() ) ){
            toreturn[2*i+1] = 0;
         } else {
            toreturn[2*i] = 0;
         }
      }
   }

   //Bl matching
   if (blmatch){
      if ( ev->mt2_220grida < ev->mt2_220gridb ){
         for (int i=4; i<8; i++){ toreturn[i] = 0; }
      }

      if ( ev->mt2_220grida > ev->mt2_220gridb ){
         for (int i=0; i<4; i++){ toreturn[i] = 0; }
      }

   } 

   return toreturn;
}

vector<bool> Fitter::MaosCut210( vector<Event>::iterator ev, TLorentzVector &nu1p, TLorentzVector &nu1m, TLorentzVector &nu2p, TLorentzVector &nu2m ){
  
   bool distcut = 0;
   bool etadisamb = 0;
   bool blmatch = 0;
   double endpt210 = 80.4;

   if (maoscuts210 == 1){ distcut = 1; }
   else if (maoscuts210 == 2){ etadisamb = 1; }
   else if (maoscuts210 == 3){distcut = 1; etadisamb = 1; }
   else if (maoscuts210 == 4){blmatch = 1; }
   else if (maoscuts210 == 5){distcut = 1; blmatch = 1; }
   else if (maoscuts210 == 6){etadisamb = 1; blmatch = 1; }
   else if (maoscuts210 == 7){distcut = 1; etadisamb = 1; blmatch = 1; }
 
   vector<bool> toreturn;
   toreturn.clear();
   for (int i=0; i<8; i++){
      toreturn.push_back(true);
   }

   TLorentzVector neuarray [] = { nu1p, nu1m, nu2p, nu2m };

   //Nans
   for (int i=0; i<4; i++){
      if ( isnan( (neuarray[i]).Pz() ) ){ toreturn[i] = 0; toreturn[i+4] = 0; }
   }

   //Near Endpoint
   if (distcut){
      if ( (ev->mt2_210grid < endpt210-30) or (ev->mt2_210grid > endpt210+20) ){
         for (int i=0; i<8; i++){ toreturn[i] = 0; }
      }
   }

   //Eta 1
   if (etadisamb){
      for (int i=0; i<2; i++){
         if ( isnan( (neuarray[2*i]).Pz() ) or isnan( (neuarray[2*i+1]).Pz() ) ) continue;
         if ( abs( (neuarray[2*i]).Eta() ) < abs( (neuarray[2*i+1]).Eta() ) ){
            toreturn[2*i+1] = 0;
            toreturn[2*i+1+4] = 0;
         } else {
            toreturn[2*i] = 0;
            toreturn[2*i+4] = 0;
         }
      }
   }

   //Bl matching
   if (blmatch){
      if ( ev->mt2_220grida < ev->mt2_220gridb ){
         for (int i=4; i<8; i++){ toreturn[i] = 0; }
      }

      if ( ev->mt2_220grida > ev->mt2_220gridb ){
         for (int i=0; i<4; i++){ toreturn[i] = 0; }
      }

   } 

   return toreturn;
}

void Fitter::PlotTemplates( vector< map< string, map<string, TH1D*> > >& hists_ ){

   std::string pathstr;
   char* path = std::getenv("WORKING_DIR");
   if (path==NULL) {
      pathstr = "./results";
   }else {
      pathstr = path;
   }

   TFile *fileout = new TFile( (pathstr+"/plotsTemplates.root").c_str(), "RECREATE" );
   fileout->cd();


   // templates for all masses
   string sb[] = {"sig","bkg"};

   for( map<string, Distribution>::iterator it = dists.begin(); it != dists.end(); it++ ){

      string name = it->first;
      Distribution *dist = &(it->second);

      if( dist->activate ){// only do this if we're fitting the variable in question

         TDirectory *dir = fileout->mkdir( name.c_str() );
         dir->cd();

         for(unsigned int k=0; k < sizeof(sb)/sizeof(sb[0]); k++){ // sig,bkg
            for(int j=0; j < NMP; j++){ // masses

               stringstream ssmass;
               ssmass << floor(masspoints[j]);
               string smass = ssmass.str();

               TCanvas *canvas = new TCanvas( ("c"+sb[k]+"_"+name+smass).c_str(),
                     ("c"+sb[k]+"_"+name+smass).c_str(), 800, 800);
               canvas->SetFillColor(0);
               canvas->cd();

               TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
               TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
               pad1->SetTopMargin(0.1);
               pad1->SetBottomMargin(0.01);
               pad1->SetRightMargin(0.1);
               pad1->SetFillColor(0);
               pad2->SetTopMargin(0.01);
               pad2->SetBottomMargin(0.3);
               pad2->SetRightMargin(0.1);
               pad2->SetFillColor(0);
               pad1->Draw();
               pad2->Draw();

               // line for ratio plot
               TF1 *func = new TF1("func","[0]",-10E6,10E6);
               func->SetParameter(0,1.0);
               func->SetLineWidth(1);
               func->SetLineStyle(7);
               func->SetLineColor(1);

               // pad 1
               pad1->cd();

               TH1D *hmc;
               if( sb[k] == "sig" ){
                  hmc = (TH1D*)hists_[1][name]["ttbar"+smass+"_signal"]->Clone("hmc");
               }else{
                  hmc = (TH1D*)hists_[1][name]["ttbar"+smass+"_mistag"]->Clone("hmc");
                  hmc->Add( hists_[1][name]["ttbar"+smass+"_taus"] );
                  hmc->Add( hists_[1][name]["ttbar"+smass+"_hadronic"] );
                  hmc->Add( hists_[1][name]["other"] );
               }

               hmc->SetTitle( hmc->GetTitle()+TString(" "+sb[k]+" shape @ "+smass+".5") );
               hmc->GetXaxis()->SetTitleSize(0.00);
               hmc->GetYaxis()->SetLabelSize(0.07);
               hmc->GetYaxis()->SetTitleSize(0.08);
               hmc->GetYaxis()->SetTitleOffset(1.0);
               hmc->GetXaxis()->SetLabelFont(42);
               hmc->GetYaxis()->SetLabelFont(42);
               hmc->GetXaxis()->SetTitleFont(42);
               hmc->GetYaxis()->SetTitleFont(42);

               hmc->Scale(1.0/hmc->Integral("width"));

               hmc->SetMarkerStyle(20);
               hmc->DrawCopy();

               Shapes * fptr = new Shapes( name, dist->ptrain,
                     dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->lbnd, dist->rbnd );
               fptr->aGPsig.ResizeTo( dist->aGPsig.GetNoElements() );
               fptr->aGPsig = dist->aGPsig;
               fptr->aGPbkg.ResizeTo( dist->aGPbkg.GetNoElements() );
               fptr->aGPbkg = dist->aGPbkg;

               fptr->Ainv_sig.ResizeTo( dist->aGPsig.GetNoElements(), dist->aGPsig.GetNoElements() );
               fptr->Ainv_sig = dist->Ainv_sig;
               fptr->Ainv_bkg.ResizeTo( dist->aGPbkg.GetNoElements(), dist->aGPbkg.GetNoElements() );
               fptr->Ainv_bkg = dist->Ainv_bkg;

               TF1 *ftemplate = new TF1("ftemplate", fptr, &Shapes::Ftot, dist->lbnd, dist->rbnd, 6);
               ftemplate->SetNpx(500);

               // normalization inside likelihood function (temp)
               ftemplate->SetParameters( masspoints[j], 1.0, 1-k, 1.0, 1.0, 1.0 );
               double integralsig = (sb[k] == "sig") ? ftemplate->Integral(dist->lbnd, dist->rbnd) : 1.0;
               double integralbkg = (sb[k] == "bkg") ? ftemplate->Integral(dist->lbnd, dist->rbnd) : 1.0;
               ftemplate->SetParameters( masspoints[j], 1.0, 1-k,
                     1.0, integralsig, integralbkg );
               ftemplate->SetLineWidth(2);

               // TGraph with GP covariance
               TGraphErrors *gpvar = new TGraphErrors();
               for(int x=dist->lbnd; x < dist->rbnd; x++){
                  gpvar->SetPoint(x, x, ftemplate->Eval(x));
                  gpvar->SetPointError(x, 0, sqrt(fptr->Fmbl_gp_var(x,masspoints[j],1.0,sb[k])));
               }
               gpvar->SetLineColor(2);
               gpvar->SetFillColor(5);
               gpvar->Draw("E3");

               ftemplate->DrawCopy("same");
               hmc->DrawCopy("same"); // redraw points

               // pad 2
               pad2->cd();
               TH1D *hratio = (TH1D*)hmc->Clone("hratio");
               hratio->Divide( ftemplate );

               hratio->SetTitle("");
               hratio->GetYaxis()->SetTitle("data/mc");
               hratio->GetYaxis()->CenterTitle();
               hratio->SetStats(0);

               hratio->GetXaxis()->SetTitleSize(0.14);
               hratio->GetXaxis()->SetLabelSize(0.14);
               hratio->GetYaxis()->SetLabelSize(0.11);
               hratio->GetYaxis()->SetTitleSize(0.14);
               hratio->GetYaxis()->SetTitleOffset(0.28);
               hratio->GetXaxis()->SetLabelFont(42);
               hratio->GetYaxis()->SetLabelFont(42);
               hratio->GetXaxis()->SetTitleFont(42);
               hratio->GetYaxis()->SetTitleFont(42);
               hratio->SetMaximum( 1.6 );
               hratio->SetMinimum( 0.4 );
               hratio->GetYaxis()->SetNdivisions(505);

               hratio->Draw("EP");
               func->Draw("same");

               canvas->Write();

               delete canvas;
               delete func;
               delete ftemplate;
               delete fptr;
            }
         }
      }
   }
   fileout->cd();

   /*
   // plot template as a function of top mass
   for( map<string, Distribution>::iterator it = dists.begin(); it != dists.end(); it++ ){

      string name = it->first;
      Distribution *dist = &(it->second);

      if( dist->activate ){// only do this if we're fitting the variable in question

         TDirectory *dir = fileout->mkdir( ("mtshape_"+name).c_str() );
         dir->cd();

         for(unsigned int k=0; k < sizeof(sb)/sizeof(sb[0]); k++){ // sig,bkg
            for(double x=0; x <= dist->range; x+=10){ // bin of mbl

               stringstream ssx;
               ssx << x;
               string sx = ssx.str();

               TCanvas *canvas = new TCanvas( ("c"+sb[k]+"_"+name+sx).c_str(),
                     ("c"+sb[k]+"_"+name+sx).c_str(), 800, 800);
               canvas->SetFillColor(0);
               canvas->cd();

               // graph with template value at mbl = x
               TGraph *gtemplate = new TGraph();
               Shapes * fptr = new Shapes( name, dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->range );
               fptr->aGPsig.ResizeTo( dist->aGPsig.GetNoElements() );
               fptr->aGPsig = dist->aGPsig;
               fptr->aGPbkg.ResizeTo( dist->aGPbkg.GetNoElements() );
               fptr->aGPbkg = dist->aGPbkg;
               TF1 *ftemplate = new TF1("ftemplate", fptr, &Shapes::Ftot, 0, dist->range, 5);
               int count=0;
               for(double m=160.0; m <= 183.0; m+=0.5){ // value of mt
                  // normalization inside likelihood function (temp)
                  ftemplate->SetParameters( m, 1-k, 1.0, 1.0, 1.0 );
                  double integralsig = (sb[k] == "sig") ? ftemplate->Integral(0,dist->range) : 1.0;
                  double integralbkg = (sb[k] == "bkg") ? ftemplate->Integral(0,dist->range) : 1.0;
                  ftemplate->SetParameters( m, 1-k, 1.0, integralsig, integralbkg );

                  gtemplate->SetPoint(count, m, ftemplate->Eval(x));
                  count++;
               }
               gtemplate->SetTitle( hists_[name]["ttbar172_signal"]->GetTitle()
                     + TString(" "+sb[k]+" shape @ "+dist->title+" = "+sx) );
               gtemplate->SetLineColor(2);
               gtemplate->SetLineWidth(2);

               // now do the same at mc masspoints
               TGraphErrors *gmc = new TGraphErrors();
               count=0;
               for(int j=0; j < NMP; j++){
                  stringstream ssmass;
                  ssmass << floor(masspoints[j]);
                  string smass = ssmass.str();

                  TH1D *hmc;
                  if( sb[k] == "sig" ){
                     hmc = (TH1D*)hists_[name]["ttbar"+smass+"_signal"]->Clone("hmc");
                  }else{
                     hmc = (TH1D*)hists_[name]["ttbar"+smass+"_mistag"]->Clone("hmc");
                     hmc->Add( hists_[name]["ttbar"+smass+"_taus"] );
                     hmc->Add( hists_[name]["ttbar"+smass+"_hadronic"] );
                     hmc->Add( hists_[name]["other"] );
                  }
                  hmc->Scale( 1.0/hmc->Integral("width") );

                  gmc->SetPoint(count, masspoints[j], hmc->GetBinContent(hmc->FindBin(x)) );
                  gmc->SetPointError(count, 0.0, hmc->GetBinError(hmc->FindBin(x)) );
                  count++;
               }

               gmc->SetMarkerStyle(20);

               gmc->SetMinimum( min(gtemplate->GetMinimum(),gmc->GetMinimum()) );
               gmc->SetMaximum( max(gtemplate->GetMaximum(),gmc->GetMaximum()) );
               gmc->Draw("AEP");
               gtemplate->Draw("same C");
               gmc->Draw("EP");

               canvas->Write();

               delete canvas;
               delete ftemplate;
               delete fptr;
            }
         }
      }
   }
   */

   // TODO
   // need to figure out variance band for this plot
   /*
      TDirectory *dir = fileout->mkdir( "mtshape" );
      dir->cd();
      for(unsigned int k=0; k < sizeof(sb)/sizeof(sb[0]); k++){ // sig,bkg
      for(double x=0; x <= rangembl; x+=10){ // bin of mbl
      cout << "mbl: " << x << endl;

      stringstream ssx;
      ssx << x;
      string sx = ssx.str();

      TCanvas *canvas = new TCanvas( ("c"+sb[k]+"_mbl"+sx).c_str(),
      ("c"+sb[k]+"_mbl"+sx).c_str(), 800, 800);
      canvas->SetFillColor(0);
      canvas->cd();

   // graph with template value at mbl = x
   TGraphErrors *gtemplate = new TGraphErrors();
   Shapes * fptr = new Shapes( gplength_mbl, gplength_mt, lbnd, rbnd, gnorm1, gnorm2 );
   fptr->aGPsig.ResizeTo( aGPsig.GetNoElements() );
   fptr->aGPsig = aGPsig;
   fptr->aGPbkg.ResizeTo( aGPbkg.GetNoElements() );
   fptr->aGPbkg = aGPbkg;

   fptr->Ainv_sig.ResizeTo( aGPsig.GetNoElements(), aGPsig.GetNoElements() );
   fptr->Ainv_sig = Ainv_sig;
   fptr->Ainv_bkg.ResizeTo( aGPbkg.GetNoElements(), aGPbkg.GetNoElements() );
   fptr->Ainv_bkg = Ainv_bkg;

   TF1 *ftemplate = new TF1("ftemplate", fptr, &Shapes::Fmbl_tot, 0, rangembl, 5);
   int count=0;
   for(double m=160.0; m <= 183.0; m+=0.5){ // value of mt
   // normalization inside likelihood function (temp)
   ftemplate->SetParameters( m, 1-k, 1.0, 1.0, 1.0 );
   double integralsig = (sb[k] == "sig") ? ftemplate->Integral(0,rangembl) : 1.0;
   double integralbkg = (sb[k] == "bkg") ? ftemplate->Integral(0,rangembl) : 1.0;
   ftemplate->SetParameters( m, 1-k, 1.0, integralsig, integralbkg );

   gtemplate->SetPoint(count, m, ftemplate->Eval(x));
   gtemplate->SetPointError(count, 0, sqrt(fptr->Fmbl_gp_var(x,m,sb[k])));
   count++;
   }
   gtemplate->SetTitle( TString("M_{bl} "+sb[k]+" shape @ mbl = "+sx) );
   gtemplate->SetLineColor(2);
   gtemplate->SetLineWidth(2);
   gtemplate->SetFillStyle(3004);
   gtemplate->SetFillColor(4);

   // now do the same at mc masspoints
   TGraphErrors *gmc = new TGraphErrors();
   count=0;
   for(int j=0; j < NMP; j++){
   stringstream ssmass;
   ssmass << floor(masspoints[j]);
   string smass = ssmass.str();

   TH1D *hmc;
   if( sb[k] == "sig" ){
   hmc = (TH1D*)hists_["mbl"]["ttbar"+smass+"_signal"]->Clone("hmc");
   }else{
   hmc = (TH1D*)hists_["mbl"]["ttbar"+smass+"_mistag"]->Clone("hmc");
   hmc->Add( hists_["mbl"]["ttbar"+smass+"_taus"] );
   hmc->Add( hists_["mbl"]["ttbar"+smass+"_hadronic"] );
   hmc->Add( hists_["mbl"]["other"] );
   }
   }
   }
   }
   */

   for( map<string, Distribution>::iterator it = dists.begin(); it != dists.end(); it++ ){

      string name = it->first;
      Distribution *dist = &(it->second);

      if( dist->activate ){// only do this if we're fitting the variable in question
         fileout->cd();

         TCanvas *cmbl_signal = new TCanvas( ("c_"+name+"_signal").c_str(),
               (dist->title+" Template").c_str(), 800, 800);
         cmbl_signal->cd();

         // mass points
         TH1D* mbl166 = (TH1D*)hists_[1][name]["ttbar166_signal"]->Clone("mbl166");
         TH1D* mbl172 = (TH1D*)hists_[1][name]["ttbar172_signal"]->Clone("mbl172");
         TH1D* mbl178 = (TH1D*)hists_[1][name]["ttbar178_signal"]->Clone("mbl178");

         mbl166->Scale( 1.0/mbl166->Integral("width") );
         mbl172->Scale( 1.0/mbl172->Integral("width") );
         mbl178->Scale( 1.0/mbl178->Integral("width") );

         mbl166->SetLineColor(2);
         mbl172->SetLineColor(1);
         mbl178->SetLineColor(3);
         mbl166->SetMarkerColor(2);
         mbl172->SetMarkerColor(1);
         mbl178->SetMarkerColor(3);

         mbl166->DrawCopy();
         mbl172->DrawCopy("same");
         mbl178->DrawCopy("same");

         // mbl likelihood

         Shapes * fptr = new Shapes( name, dist->ptrain,
               dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->lbnd, dist->rbnd );
         fptr->aGPsig.ResizeTo( dist->aGPsig.GetNoElements() );
         fptr->aGPsig = dist->aGPsig;
         fptr->aGPbkg.ResizeTo( dist->aGPbkg.GetNoElements() );
         fptr->aGPbkg = dist->aGPbkg;
         TF1 *fmbl_tot = new TF1( ("f"+name+"_tot").c_str(), fptr, &Shapes::Ftot, dist->lbnd, dist->rbnd, 6);

         fmbl_tot->SetParameters( 166.5, 1.0, 1.0, 1.0, 1.0, 1.0 );
         fmbl_tot->SetParameters( 166.5, 1.0, 1.0, 1.0, fmbl_tot->Integral(dist->lbnd, dist->rbnd), 1.0 );
         fmbl_tot->SetLineColor(2);
         fmbl_tot->DrawCopy("same");

         fmbl_tot->SetParameters( 172.5, 1.0, 1.0, 1.0, 1.0, 1.0 );
         fmbl_tot->SetParameters( 172.5, 1.0, 1.0, 1.0, fmbl_tot->Integral(dist->lbnd, dist->rbnd), 1.0 );
         fmbl_tot->SetLineColor(1);
         fmbl_tot->DrawCopy("same");

         fmbl_tot->SetParameters( 178.5, 1.0, 1.0, 1.0, 1.0, 1.0 );
         fmbl_tot->SetParameters( 178.5, 1.0, 1.0, 1.0, fmbl_tot->Integral(dist->lbnd, dist->rbnd), 1.0 );
         fmbl_tot->SetLineColor(3);
         fmbl_tot->DrawCopy("same");

         TLegend *lm = new TLegend(0.76,0.59,0.98,0.77);
         lm->SetFillStyle(0);
         lm->SetBorderSize(0);
         lm->AddEntry( mbl166, "166.5" );
         lm->AddEntry( mbl172, "172.5" );
         lm->AddEntry( mbl178, "178.5" );
         lm->Draw("same");

         cmbl_signal->Write();

         delete fptr;
         delete cmbl_signal;
         delete fmbl_tot;
         delete lm;

         TCanvas *cmbl_signal_jfact = new TCanvas( ("c_"+name+"_signal_jfact").c_str(),
               (dist->title+" Template").c_str(), 800, 800);
         cmbl_signal_jfact->cd();

         // jfactor points
         TH1D* mblDN = (TH1D*)hists_[0][name]["ttbar172_signal"]->Clone("mblDN");
         TH1D* mblCENT = (TH1D*)hists_[1][name]["ttbar172_signal"]->Clone("mblCENT");
         TH1D* mblUP = (TH1D*)hists_[2][name]["ttbar172_signal"]->Clone("mblUP");

         mblDN->Scale( 1.0/mblDN->Integral("width") );
         mblCENT->Scale( 1.0/mblCENT->Integral("width") );
         mblUP->Scale( 1.0/mblUP->Integral("width") );

         mblDN->SetLineColor(2);
         mblCENT->SetLineColor(1);
         mblUP->SetLineColor(3);
         mblDN->SetMarkerColor(2);
         mblCENT->SetMarkerColor(1);
         mblUP->SetMarkerColor(3);

         mblDN->DrawCopy();
         mblCENT->DrawCopy("same");
         mblUP->DrawCopy("same");

         // likelihood
         Shapes * fptr_jfact = new Shapes( name, dist->ptrain,
               dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->lbnd, dist->rbnd );
         fptr_jfact->aGPsig.ResizeTo( dist->aGPsig.GetNoElements() );
         fptr_jfact->aGPsig = dist->aGPsig;
         fptr_jfact->aGPbkg.ResizeTo( dist->aGPbkg.GetNoElements() );
         fptr_jfact->aGPbkg = dist->aGPbkg;
         TF1 *fmbl_tot_jfact = new TF1( ("f"+name+"_tot").c_str(), fptr_jfact, &Shapes::Ftot, dist->lbnd, dist->rbnd, 6);

         fmbl_tot_jfact->SetParameters( 172.5, jfactpoints[0], 1.0, 1.0, 1.0, 1.0 );
         fmbl_tot_jfact->SetParameters( 172.5, jfactpoints[0], 1.0, 1.0, fmbl_tot_jfact->Integral(dist->lbnd, dist->rbnd), 1.0 );
         fmbl_tot_jfact->SetLineColor(2);
         fmbl_tot_jfact->DrawCopy("same");

         fmbl_tot_jfact->SetParameters( 172.5, jfactpoints[1], 1.0, 1.0, 1.0, 1.0 );
         fmbl_tot_jfact->SetParameters( 172.5, jfactpoints[1], 1.0, 1.0, fmbl_tot_jfact->Integral(dist->lbnd, dist->rbnd), 1.0 );
         fmbl_tot_jfact->SetLineColor(1);
         fmbl_tot_jfact->DrawCopy("same");

         fmbl_tot_jfact->SetParameters( 172.5, jfactpoints[2], 1.0, 1.0, 1.0, 1.0 );
         fmbl_tot_jfact->SetParameters( 172.5, jfactpoints[2], 1.0, 1.0, fmbl_tot_jfact->Integral(dist->lbnd, dist->rbnd), 1.0 );
         fmbl_tot_jfact->SetLineColor(3);
         fmbl_tot_jfact->DrawCopy("same");

         TLegend *lm_jfact = new TLegend(0.76,0.59,0.98,0.77);
         lm_jfact->SetFillStyle(0);
         lm_jfact->SetBorderSize(0);
         lm_jfact->AddEntry( mblDN, "j=0.98" );
         lm_jfact->AddEntry( mblCENT, "j=1.00" );
         lm_jfact->AddEntry( mblUP, "j=1.02" );
         lm_jfact->Draw("same");

         cmbl_signal_jfact->Write();

         delete fptr_jfact;
         delete cmbl_signal_jfact;
         delete fmbl_tot_jfact;
         delete lm_jfact;

      }
   }


   fileout->Close();

   return;

}


