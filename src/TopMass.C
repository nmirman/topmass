#include "TopMass.h"
#include "Shapes.h"
#include "Mt2Calculator.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TRandom3.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"

using namespace std;


//
// constructor and destructor
//

Fitter::Fitter(){

   // MINUIT variables
   gMinuit = 0;
   fFunc = 0;

   // compute weighted errors
   TH1::SetDefaultSumw2();

   // plot formatting
   gROOT->ProcessLineSync(".L scripts/tdrstyle.C");
   gROOT->ProcessLineSync("setTDRStyle()");

   compute_profile = false;

}

Fitter::~Fitter(){
      if (gMinuit) delete gMinuit;
      if (fFunc) delete fFunc;
}

void Fitter::InitializeDists(){

   // gaussian process length scales
   dists[ "mbl" ] = Distribution( "mbl", "M_{bl}", 5.8, 10.7, 11.8, 26.1, 300 );
   dists[ "mt2_220_nomatchmbl" ] = Distribution( "mt2_220_nomatchmbl", "M_{T2} 220", 3.9, 8.4, 8.5, 12.5, 300 );
   dists[ "maos220blv" ] = Distribution( "maos220blv","blv mass from Maos neutrinos from M_{T2} 220", 1.6, 6.4, 19.2, 19.2, 500 );
   dists[ "maos210blv" ] = Distribution( "maos210blv","blv mass from Maos neutrinos from M_{T2} 210", 0.60, 8.5, 19.1, 18.9, 500 );

}

//
// member definitions
//

void Fitter::LoadDatasets( map<string, Dataset>& datasets ){

   // check machine location
   char hostname[1024];
   gethostname(hostname, 1024);
   char *pch;
   pch = strstr( hostname, "lpc" );

   // file path
   string path;
   if( pch != NULL ) path = "root://cmseos:1094//eos/uscms/store/user/nmirman/Ntuples/TopMass/20141030/";
   if( pch == NULL ) path = "/afs/cern.ch/work/n/nmirman/public/Ntuples/TopMass/20150226/";

   // filenames
   datasets[ "data" ]      = Dataset( path, "ntuple_data.root" );
   datasets[ "t_tw" ]      = Dataset( path, "ntuple_T_tW.root" );
   datasets[ "tbar_tw" ]   = Dataset( path, "ntuple_Tbar_tW.root" );
   datasets[ "dy" ]        = Dataset( path, "ntuple_DYJetsToLL.root" );
   datasets[ "wjets" ]     = Dataset( path, "ntuple_WJetsToLNu.root" );
   datasets[ "ww" ]        = Dataset( path, "ntuple_WW.root" );
   datasets[ "wz" ]        = Dataset( path, "ntuple_WZ.root" );
   datasets[ "zz" ]        = Dataset( path, "ntuple_ZZ.root" );
   datasets[ "ttbar161" ]  = Dataset( path, "ntuple_TTJets_mass161_5.root" ); 
   datasets[ "ttbar163" ]  = Dataset( path, "ntuple_TTJets_mass163_5.root" ); 
   datasets[ "ttbar166" ]  = Dataset( path, "ntuple_TTJets_mass166_5.root" ); 
   datasets[ "ttbar169" ]  = Dataset( path, "ntuple_TTJets_mass169_5.root" ); 
   datasets[ "ttbar172" ]  = Dataset( path, "ntuple_TTJets.root" );
   datasets[ "ttbar175" ]  = Dataset( path, "ntuple_TTJets_mass175_5.root" ); 
   datasets[ "ttbar178" ]  = Dataset( path, "ntuple_TTJets_mass178_5.root" ); 
   datasets[ "ttbar181" ]  = Dataset( path, "ntuple_TTJets_mass181_5.root" ); 
   datasets[ "ttbarsyst_scaleup" ]        = Dataset( path, "ntuple_TTJets_scaleup.root" );
   datasets[ "ttbarsyst_scaledown" ]      = Dataset( path, "ntuple_TTJets_scaledown.root" );
   datasets[ "ttbarsyst_matchingup" ]     = Dataset( path, "ntuple_TTJets_matchingup.root" );
   datasets[ "ttbarsyst_matchingdown" ]   = Dataset( path, "ntuple_TTJets_matchingdown.root" );

   // for mc weights
   datasets[ "dy" ].mc_nevts         = 30459503;
   datasets[ "ww" ].mc_nevts         = 10000431;
   datasets[ "wz" ].mc_nevts         = 10000283;
   datasets[ "zz" ].mc_nevts         = 9799908;
   datasets[ "tbar_tw" ].mc_nevts    = 493460;
   datasets[ "t_tw" ].mc_nevts       = 497658;
   datasets[ "wjets" ].mc_nevts      = 57709905;
   datasets[ "ttbar161" ].mc_nevts   = 5369214;
   datasets[ "ttbar163" ].mc_nevts   = 5365348;
   datasets[ "ttbar166" ].mc_nevts   = 4469095;
   datasets[ "ttbar169" ].mc_nevts   = 5202817;
   datasets[ "ttbar172" ].mc_nevts   = 6923750;
   datasets[ "ttbar175" ].mc_nevts   = 5186494;
   datasets[ "ttbar178" ].mc_nevts   = 4733483;
   datasets[ "ttbar181" ].mc_nevts   = 5145140;
   datasets[ "ttbarsyst_scaleup" ].mc_nevts        = 5009488;
   datasets[ "ttbarsyst_scaledown" ].mc_nevts      = 5387181;
   datasets[ "ttbarsyst_matchingup" ].mc_nevts     = 5415010;
   datasets[ "ttbarsyst_matchingdown" ].mc_nevts   = 5476728;

   datasets[ "dy" ].mc_xsec          = 3351.97;
   datasets[ "ww" ].mc_xsec          = 54.838;
   datasets[ "wz" ].mc_xsec          = 33.21;
   datasets[ "zz" ].mc_xsec          = 8.059;
   datasets[ "tbar_tw" ].mc_xsec     = 11.1;
   datasets[ "t_tw" ].mc_xsec        = 11.1;
   datasets[ "wjets" ].mc_xsec       = 37509.0;
   datasets[ "ttbar161" ].mc_xsec    = 234;
   datasets[ "ttbar163" ].mc_xsec    = 234;
   datasets[ "ttbar166" ].mc_xsec    = 234;
   datasets[ "ttbar169" ].mc_xsec    = 234;
   datasets[ "ttbar172" ].mc_xsec    = 234;
   datasets[ "ttbar175" ].mc_xsec    = 234;
   datasets[ "ttbar178" ].mc_xsec    = 234;
   datasets[ "ttbar181" ].mc_xsec    = 234;
   datasets[ "ttbarsyst_scaleup" ].mc_xsec         = 234;
   datasets[ "ttbarsyst_scaledown" ].mc_xsec       = 234;
   datasets[ "ttbarsyst_matchingup" ].mc_xsec      = 234;
   datasets[ "ttbarsyst_matchingdown" ].mc_xsec    = 234;

   return;
}

void Fitter::ReadNtuple( string path, string process, double mcweight, 
      string selection, vector<Event>& eventvec, int opt, double fracevts,
      int statval_numPE, int statval_PE ){
   
   // declare variables
   TLorentzVector *jet1 = new TLorentzVector();
   TLorentzVector *jet2 = new TLorentzVector();
   TLorentzVector *lep1 = new TLorentzVector();
   TLorentzVector *lep2 = new TLorentzVector();
   TLorentzVector *met = new TLorentzVector();

   int nvert;
   double jet1PtRes, jet1PhiRes, jet1EtaRes, jet2PtRes, jet2PhiRes, jet2EtaRes;
   int lpPdgIdGEN, lmPdgIdGEN, nPdgIdGEN, nbPdgIdGEN;
   int jet1GenId, jet2GenId;
   float weight_pu=1, weight_toppt=1, weight_btag=1, weight_mu=1, weight_elec=1;
   int nmuons, nelectrons;
   vector<float> *pdf_weights = 0;

   // open ntuple
   TFile file( path.c_str() );
   TTree *tree = (TTree*)file.Get(selection.c_str());
   TBranch *branch = 0;

   tree->SetBranchAddress("jet1FourVector", &jet1);
   tree->SetBranchAddress("jet2FourVector", &jet2);
   tree->SetBranchAddress("lep1FourVector", &lep1);
   tree->SetBranchAddress("lep2FourVector", &lep2);
   tree->SetBranchAddress("metFourVector", &met);
   tree->SetBranchAddress("jet1PtResolution", &jet1PtRes);
   tree->SetBranchAddress("jet1PhiResolution", &jet1PhiRes);
   tree->SetBranchAddress("jet1EtaResolution", &jet1EtaRes);
   tree->SetBranchAddress("jet2PtResolution", &jet2PtRes);
   tree->SetBranchAddress("jet2PhiResolution", &jet2PhiRes);
   tree->SetBranchAddress("jet2EtaResolution", &jet2EtaRes);
   tree->SetBranchAddress("vertices", &nvert);

   if( process.find("ttbar") != string::npos ){
      tree->SetBranchAddress("lpPdgIdGEN", &lpPdgIdGEN);
      tree->SetBranchAddress("lmPdgIdGEN", &lmPdgIdGEN);
      tree->SetBranchAddress("nPdgIdGEN", &nPdgIdGEN);
      tree->SetBranchAddress("nbPdgIdGEN", &nbPdgIdGEN);
      tree->SetBranchAddress("jet1GenId", &jet1GenId);
      tree->SetBranchAddress("jet2GenId", &jet2GenId);
   }

   if( !(process.find("data") != string::npos) ){
      tree->SetBranchAddress("weight_pu", &weight_pu);
      tree->SetBranchAddress("weight_toppt", &weight_toppt);
      tree->SetBranchAddress("weight_btag", &weight_btag);
      tree->SetBranchAddress("weight_mu", &weight_mu);
      tree->SetBranchAddress("weight_elec", &weight_elec);
   }

   tree->SetBranchAddress("nmuons", &nmuons);
   tree->SetBranchAddress("nelectrons", &nelectrons);
   tree->SetBranchAddress("pdf_weights", &pdf_weights);

   // subset of events to use (for training, testing)
   int start = 0;
   int treesize = tree->GetEntries();
   int end = treesize;
   if( opt == 1 ){
      start = 0;
      end = treesize/3;
   }
   if( opt == 2 ){
      start = treesize/3;
      end = 2*treesize/3;
   }
   if( statval_numPE != -1 ){
      int nevts = end-start;
      int begin = start;
      start = begin + statval_PE*nevts/statval_numPE;
      end = begin + (statval_PE+1)*nevts/statval_numPE;
   }

   // run on fraction of total events
   if( fracevts != -1 ){
      end = start + fracevts*(end-start);
   }

   // fill event vector
   for(int ev=start; ev < end; ev++){
      tree->GetEntry(ev);

      Event evtemp;

      // global quantities
      evtemp.process = process;
      evtemp.weight = mcweight * weight_pu * weight_toppt * weight_btag * weight_mu * weight_elec;
      evtemp.nvertices = nvert;

      // jets, leptons, met
      evtemp.jet1 = *jet1;
      evtemp.jet2 = *jet2;

      evtemp.lep1 = *lep1;
      evtemp.lep2 = *lep2;

      evtemp.met = *met;

      evtemp.isemu = nmuons==1 and nelectrons==1;

      evtemp.pdf_weights = *pdf_weights;

      //
      // classify events
      // 
      if( process.find("data") != string::npos ){
         evtemp.type = process;
      }
      // physics backgrounds
      else if( process.find("ttbar") == string::npos ){
         evtemp.type = "other";
      }
      // hadronic decays
      else if( fabs(lpPdgIdGEN) < 5 or fabs(lmPdgIdGEN) < 5
            or fabs(nPdgIdGEN) < 5 or fabs(nbPdgIdGEN) < 5 ){
         evtemp.type = process+"_hadronic";
      }
      // tau decays
      else if( fabs(lpPdgIdGEN) == 15 or fabs(lmPdgIdGEN) == 15
            or fabs(nPdgIdGEN) == 15 or fabs(nbPdgIdGEN) == 15 ){
         evtemp.type = process+"_taus";
      }
      // mistag bkg
      else if( fabs(jet1GenId) != 5 or fabs(jet2GenId) != 5 ){
         evtemp.type = process+"_mistag";
      }
      // signal ttbar
      else if( (lpPdgIdGEN == -11 or lpPdgIdGEN == -13)
            and (lmPdgIdGEN == 11 or lmPdgIdGEN == 13)
            and (nPdgIdGEN == 12 or nPdgIdGEN == 14)
            and (nbPdgIdGEN == -12 or nbPdgIdGEN == -14)
            and (fabs(jet1GenId) == 5 and fabs(jet2GenId) == 5) ){
         evtemp.type = process+"_signal";
      }
      // check for unclassified events
      else{
         cout << "ERROR IN TTBAR EVENT CLASSIFICATION" << endl;
         break;
      }

      // push back event
      if ( (jet1->M() < 40 and jet2->M() < 40)/* and evtemp.isemu*/ )
         eventvec.push_back( evtemp );

   }

   return;
}


void Fitter::GetVariables( vector<Event>& eventvec ){

   Mt2Calculator::Calculator Calc;
   for( vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){

      Calc.SetParticles( ev->jet1, ev->jet2, ev->lep1, ev->lep2, ev->met );

      ev->mt2_221 = Calc.GetMt2(2,1);
      ev->mt2_220 = Calc.GetMt2(2,0);
      ev->mt2_210 = Calc.GetMt2(1,0);
      ev->mbls = Calc.GetBlInvariantMasses();

      //Maos
      //Declare variables to accept Maos values
      //neutrino 4-vectors
      TLorentzVector maos210_neu1p, maos210_neu1m, maos210_neu2p, maos210_neu2m;
      TLorentzVector maos220_neu1ap, maos220_neu1am, maos220_neu2ap, maos220_neu2am, maos220_neu1bp, maos220_neu1bm, maos220_neu2bp, maos220_neu2bm;

      //blv masses
      double maos210_blv1ap, maos210_blv1am, maos210_blv2ap, maos210_blv2am, maos210_blv1bp, maos210_blv1bm, maos210_blv2bp, maos210_blv2bm;
      double maos220_blv1ap, maos220_blv1am, maos220_blv2ap, maos220_blv2am, maos220_blv1bp, maos220_blv1bm, maos220_blv2bp, maos220_blv2bm;

      //calculate values and set above variables to them; then set event variables to these
      ev->mt2_210grid = Calc.MaosReturn210( maos210_neu1p, maos210_neu1m, maos210_neu2p, maos210_neu2m, maos210_blv1ap, maos210_blv1am, maos210_blv2ap, maos210_blv2am, maos210_blv1bp, maos210_blv1bm, maos210_blv2bp, maos210_blv2bm );
      std::vector<double> Mt2_220grid = Calc.MaosReturn220( maos220_neu1ap, maos220_neu1am, maos220_neu2ap, maos220_neu2am, maos220_neu1bp, maos220_neu1bm, maos220_neu2bp, maos220_neu2bm, maos220_blv1ap, maos220_blv1am, maos220_blv2ap, maos220_blv2am, maos220_blv1bp, maos220_blv1bm, maos220_blv2bp, maos220_blv2bm);
      ev->mt2_220grida = Mt2_220grid.at(0);
      ev->mt2_220gridb = Mt2_220grid.at(1);

      ev->maos210_neutrino1p = maos210_neu1p;
      ev->maos210_neutrino1m = maos210_neu1m;
      ev->maos210_neutrino2p = maos210_neu2p;
      ev->maos210_neutrino2m = maos210_neu2m;

      ev->maos220_neutrino1ap = maos220_neu1ap;
      ev->maos220_neutrino1am = maos220_neu1am;
      ev->maos220_neutrino2ap = maos220_neu2ap;
      ev->maos220_neutrino2am = maos220_neu2am;
      ev->maos220_neutrino1bp = maos220_neu1bp;
      ev->maos220_neutrino1bm = maos220_neu1bm;
      ev->maos220_neutrino2bp = maos220_neu2bp;
      ev->maos220_neutrino2bm = maos220_neu2bm;

      ev->maos210_blvmass1ap = maos210_blv1ap;
      ev->maos210_blvmass1am = maos210_blv1am;
      ev->maos210_blvmass2ap = maos210_blv2ap;
      ev->maos210_blvmass2am = maos210_blv2am;
      ev->maos210_blvmass1bp = maos210_blv1bp;
      ev->maos210_blvmass1bm = maos210_blv1bm;
      ev->maos210_blvmass2bp = maos210_blv2bp;
      ev->maos210_blvmass2bm = maos210_blv2bm;

      ev->maos220_blvmass1ap = maos220_blv1ap;
      ev->maos220_blvmass1am = maos220_blv1am;
      ev->maos220_blvmass2ap = maos220_blv2ap;
      ev->maos220_blvmass2am = maos220_blv2am;
      ev->maos220_blvmass1bp = maos220_blv1bp;
      ev->maos220_blvmass1bm = maos220_blv1bm;
      ev->maos220_blvmass2bp = maos220_blv2bp;
      ev->maos220_blvmass2bm = maos220_blv2bm;

   }

   return;
}

void Fitter::ReweightMC( vector<Event>& eventvec, string name ){
   cout << "Reweighting MC events." << endl;

   double weight_norm = 0;
   double nevts_ttbar = 0;
   for(vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){
      if( ev->type.find(name) != string::npos ){
         weight_norm += ev->weight;
         nevts_ttbar++;
      }
   }
   cout << "---> nevts = " << nevts_ttbar << " weight_norm = " << weight_norm << endl;
   // reweight
   double wgt=0;
   for(vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){
      ev->weight *= 1.0*nevts_ttbar/weight_norm;
      wgt += ev->weight;
   }
   cout << "---> sum of weights = " << wgt << endl;

}

void Fitter::Resample( vector<Event>& eventvec, int randseed ){
   cout << "Resampling events." << endl;

   // put resampled events into a new vector
   vector<Event> eventvec_resampled;

   // initialize random number engine
   TRandom3 *rand = new TRandom3(randseed);

   double maxweight = 0;
   for(vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){
      if( ev->weight > maxweight ) maxweight = ev->weight;
   }

   int numevts_data = 49243;
   // resample with replacement, taking into account event weights
   int count = 0;
   while( count < numevts_data ){
      // get event
      int ev = rand->Uniform( eventvec.size() );
      // if event weight is large enough, take event
      if( rand->Rndm() < eventvec[ev].weight/maxweight ){
         Event temp = eventvec[ev];
         eventvec_resampled.push_back( temp );
         count++;
      }
   }

   // set all event weights to unity
   for(vector<Event>::iterator ev = eventvec_resampled.begin(); ev < eventvec_resampled.end(); ev++){
      ev->weight = 1.0;
   }

   // replace eventvec with new vector
   eventvec = eventvec_resampled;
}

void Fitter::RunMinimizer( vector<Event>& eventvec ){


   gMinuit = new ROOT::Minuit2::Minuit2Minimizer ( ROOT::Minuit2::kMigrad );
   gMinuit->SetPrintLevel(3);

   // Dimension of fFunc needs to be changed if adding more variables
   fFunc = new ROOT::Math::Functor ( this, &Fitter::Min2LL, 5 );
   gMinuit->SetFunction( *fFunc );
   gMinuit->SetVariable(0, "topMass", 175.0, 0.1);

   // If we're fitting mbl, set mbl background as a limited variable, otherwise set it as a fixed variable
   if (dists["mbl"].activate){
      gMinuit->SetLimitedVariable(1, "norm", 0.7, 0.1, 0, 1.0);
   } else {
      gMinuit->SetFixedVariable(1, "norm", 0.70712);
   }

   // If we're fitting 220, set 220 background as a limited variable, otherwise set it as a fixed variable
   if (dists["mt2_220_nomatchmbl"].activate){
      gMinuit->SetLimitedVariable(2, "norm220", 0.7, 0.1, 0, 1.0);
   } else {
      gMinuit->SetFixedVariable(2, "norm220", 0.70712);
   }

   //MAOS 220
   if (dists["maos220blv"].activate){
      gMinuit->SetLimitedVariable(3, "norm_maos220", 0.7, 0.1, 0, 1.0);
   } else {
      gMinuit->SetFixedVariable(3, "norm_maos220", 0.70712);
   }

   //MAOS 210
   if (dists["maos210blv"].activate){
      gMinuit->SetLimitedVariable(4, "norm_maos210", 0.7, 0.1, 0, 1.0);
   } else {
      gMinuit->SetFixedVariable(4, "norm_maos210", 0.70712);
   }

   // set event vector and minimize
   eventvec_fit = &eventvec;

   cout << "\nFitting " << eventvec_fit->size() << " events." << endl;
   gMinuit->Minimize();
   //cout << "\nComputing Hessian." << endl;
   //gMinuit->Hesse();
   double emtLow=0, emtUp=0;
   gMinuit->GetMinosError(0,emtLow,emtUp);
   cout << "MINOS ERROR: -" << emtLow << " +" << emtUp << endl;

   return;
}

double Fitter::Min2LL(const double *x){

   double m2ll = 0;

   for( map<string, Distribution>::iterator it = dists.begin(); it != dists.end(); it++ ){

      string name = it->first;
      Distribution *dist = &(it->second);
      int iparam = -1;
      if( name.compare("mbl") == 0 ) iparam = 1;
      if( name.compare("mt2_220_nomatchmbl") == 0 ) iparam = 2;
      if( name.compare("maos220blv") == 0 ) iparam = 3;
      if( name.compare("maos210blv") == 0 ) iparam = 4;

      if( dist->activate ){// only do this if we're fitting the variable in question

         // normalization inside likelihood function (temp)
         Shapes * fptr = new Shapes( name, dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->range );
         fptr->aGPsig.ResizeTo( dist->aGPsig.GetNoElements() );
         fptr->aGPsig = dist->aGPsig;
         fptr->aGPbkg.ResizeTo( dist->aGPbkg.GetNoElements() );
         fptr->aGPbkg = dist->aGPbkg;

         TF1 *fshape_tot = new TF1( ("f"+name+"_tot").c_str(), fptr, &Shapes::Ftot, 0, dist->range, 5);
         fshape_tot->SetParameters( x[0], 1.0, 1.0, 1.0, 1.0 );
         double integralsig = fshape_tot->Integral(0,dist->range);
         fshape_tot->SetParameters( x[0], 0.0, 1.0, 1.0, 1.0 );
         double integralbkg = fshape_tot->Integral(0,dist->range);
         delete fshape_tot;
         delete fptr;

         Shapes shape( name, dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->range );
         shape.aGPsig.ResizeTo( dist->aGPsig.GetNoElements() );
         shape.aGPsig = dist->aGPsig;
         shape.aGPbkg.ResizeTo( dist->aGPbkg.GetNoElements() );
         shape.aGPbkg = dist->aGPbkg;

         double pfit [] = {x[0], x[iparam], 1.0, integralsig, integralbkg};

         // evaluate likelihood
         for( vector<Event>::iterator ev = eventvec_fit->begin(); ev < eventvec_fit->end(); ev++ ){
            if( !(ev->fit_event) ) continue;

            if ( name.compare("mbl") == 0 ){ // for mbl
               for( unsigned int j=0; j < ev->mbls.size(); j++ ){
                  if( ev->mbls[j] > dist->range ) continue;
                  //if( ev->mbls[j] > lbnd and ev->mbls[j] < rbnd ) continue;
                  double val = shape.Ftot( &(ev->mbls[j]), pfit );
                  m2ll -= 2.0*ev->weight*log( val );
               }
            }
            else if ( name.compare("mt2_220_nomatchmbl") == 0 ){ // for 220
               bool matchmbl = false;
               for ( unsigned int j=0; j < ev->mbls.size(); j++){
                  if ( ev->mbls[j] == ev->mt2_220 ) matchmbl = true;
               }
               if( matchmbl ) continue;
               if( ev->mt2_220 > dist->range ) continue;
               //if( ev->mt2_220 > lbnd and ev->mt2_220 < rbnd ) continue;
               double val = shape.Ftot( &(ev->mt2_220), pfit );
               m2ll -= 2.0*ev->weight*log( val );
            }
            else if ( name.compare("maos220blv") == 0 ){ // for Maos 220
               double blv220array [] = { ev->maos220_blvmass1ap, ev->maos220_blvmass1am, ev->maos220_blvmass2ap, ev->maos220_blvmass2am, ev->maos220_blvmass1bp, ev->maos220_blvmass1bm, ev->maos220_blvmass2bp, ev->maos220_blvmass2bm };

               vector<bool> useMaos220 = MaosCut220( ev );
               for ( unsigned int j=0; j < sizeof(blv220array)/sizeof(blv220array[0]); j++){           
                  if( blv220array[j] > dist->range ) continue;
                  //if( blv_array[j] > lbnd and blv_array[j] < rbnd ) continue;
                  if (useMaos220[j]){
                     double val = shape.Ftot( &(blv220array[j]), pfit );
                     m2ll -= 2.0*ev->weight*log( val );
                  }

               }
            }
            else if ( name.compare("maos210blv") == 0 ){ // for Maos 210
               double blv210array [] = { ev->maos210_blvmass1ap, ev->maos210_blvmass1am, ev->maos210_blvmass2ap, ev->maos210_blvmass2am, ev->maos210_blvmass1bp, ev->maos210_blvmass1bm, ev->maos210_blvmass2bp, ev->maos210_blvmass2bm };

               vector<bool> useMaos210 = MaosCut210( ev );
               for ( unsigned int j=0; j < sizeof(blv210array)/sizeof(blv210array[0]); j++){           
                  if( blv210array[j] > dist->range ) continue;
                  //if( blv_array[j] > lbnd and blv_array[j] < rbnd ) continue;
                  if (useMaos210[j]){
                     double val = shape.Ftot( &(blv210array[j]), pfit );
                     m2ll -= 2.0*ev->weight*log( val );
                  }

               }
            }

         }
      }
   }

   return m2ll;
}

void Fitter::PlotResults( map< string, map<string, TH1D*> >& hists_ ){
   cout << "Plotting fit results." << endl;

   const double *xmin = gMinuit->X();
   const double *xerr = gMinuit->Errors();
   double minvalue = gMinuit->MinValue();

   // plot fit results
   std::string pathstr;
   char* path = std::getenv("WORKING_DIR");
   if (path==NULL) {
      pathstr = "./results";
   }else {
      pathstr = path;
   }

   TFile *fileout = new TFile( (pathstr+"/plotsFitResults.root").c_str() , "RECREATE" );
   fileout->cd();

   double xmin1s [] = {xmin[1], xmin[1], xmin[2], xmin[3], xmin[4]}; // the background for each variable, as it is positioned in minuit's parameter vector
   double xerr1s [] = {xerr[1], xerr[1], xerr[2], xerr[3], xerr[4]}; // the background for each variable, as it is positioned in minuit's parameter vector

   // loop over distributions
   double chi2 = 0;
   for( map<string, Distribution>::iterator it = dists.begin(); it != dists.end(); it++ ){

      string name = it->first;
      Distribution *dist = &(it->second);
      int iparam = -1;
      if( name.compare("mbl") == 0 ) iparam = 1;
      if( name.compare("mt2_220_nomatchmbl") == 0 ) iparam = 2;
      if( name.compare("maos220blv") == 0 ) iparam = 3;
      if( name.compare("maos210blv") == 0 ) iparam = 4;

      if( dist->activate ){// only do this if we're fitting the variable in question

         // normalization inside likelihood function (temp)
         Shapes * fptr = new Shapes( name, dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->range );
         fptr->aGPsig.ResizeTo( dist->aGPsig.GetNoElements() );
         fptr->aGPsig = dist->aGPsig;
         fptr->aGPbkg.ResizeTo( dist->aGPbkg.GetNoElements() );
         fptr->aGPbkg = dist->aGPbkg;

         TF1 *ftemplate = new TF1( "ftemplate", fptr, &Shapes::Ftot, 0, dist->range, 5);

         TCanvas *canvas = new TCanvas( ("c_"+name).c_str(), ("c_"+name).c_str(), 800, 800);
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

         TH1D *hdata = (TH1D*)hists_[name]["fitevts"]->Clone("fitevts");
         hdata->GetXaxis()->SetTitleSize(0.00);
         hdata->GetYaxis()->SetLabelSize(0.07);
         hdata->GetYaxis()->SetTitleSize(0.08);
         hdata->GetYaxis()->SetTitleOffset(1.0);
         hdata->GetXaxis()->SetLabelFont(42);
         hdata->GetYaxis()->SetLabelFont(42);
         hdata->GetXaxis()->SetTitleFont(42);
         hdata->GetYaxis()->SetTitleFont(42);

         hdata->SetMarkerStyle(20);
         hdata->Draw();

         // normalization inside likelihood function (temp)
         ftemplate->SetParameters( xmin[0], 1.0, 1.0, 1.0, 1.0 );
         double integralsig = ftemplate->Integral(0,dist->range);
         ftemplate->SetParameters( xmin[0], 0.0, 1.0, 1.0, 1.0 );
         double integralbkg = ftemplate->Integral(0,dist->range);
         ftemplate->SetParameters( xmin[0], xmin1s[iparam],
               hdata->Integral("width"), integralsig, integralbkg );

         ftemplate->SetLineWidth(2);
         ftemplate->DrawCopy("same");
         hdata->DrawCopy("same"); // redraw points

         for(int n=0; n <= hdata->GetNbinsX(); n++){
            double bincontent = hdata->GetBinContent(n);
            double binerr = hdata->GetBinError(n);
            double feval = ftemplate->Eval(hdata->GetBinCenter(n));
            if( binerr == 0 ) binerr = 1;
            chi2 += pow( (bincontent-feval)/binerr, 2);
         }

         // pad 2
         pad2->cd();
         TH1D *hratio = (TH1D*)hists_[name]["fitevts"]->Clone("hratio");
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

         // This part (the profile computation may not yet be working correctly, or may not make sense.
         if ( compute_profile ){
            //
            // plot likelihood near minimum
            //
            unsigned int npnts_mt = 10;
            unsigned int npnts_kmbl = 10;
            double mt_lrange = xmin[0]-3*xerr[0];
            double mt_rrange = xmin[0]+3*xerr[0];
            double kmbl_lrange = xmin1s[iparam]-3*xerr1s[iparam];
            double kmbl_rrange = xmin1s[iparam]+3*xerr1s[iparam];

            // mt profile
            TGraph *gLmt = new TGraph();

            for(unsigned int k=0; k <= npnts_mt; k++){
               cout << "mt profile, pnt " << k << endl;
               double mt = mt_lrange + (mt_rrange-mt_lrange)*k/npnts_mt;

               // normalization inside likelihood function (temp)
               ftemplate->SetParameters( mt, 1.0, 1.0, 1.0, 1.0 );
               double integralsigk = ftemplate->Integral(0,dist->range);
               ftemplate->SetParameters( mt, 0.0, 1.0, 1.0, 1.0 );
               double integralbkgk = ftemplate->Integral(0,dist->range);
               ftemplate->SetParameters( mt, xmin1s[iparam], 1.0, integralsigk, integralbkgk );

               const double par [] = {mt, xmin1s[iparam], 1.0, integralsigk, integralbkgk};
               gLmt->SetPoint(k, mt, Min2LL(par) - minvalue);
            }

            // kvariable profile
            TGraph *gLkmbl = new TGraph();

            for(unsigned int k=0; k <= npnts_kmbl; k++){
               cout << "kmbl profile, pnt " << k << endl;
               double kmbl = kmbl_lrange + (kmbl_rrange-kmbl_lrange)*k/npnts_kmbl;

               // normalization inside likelihood function (temp)
               ftemplate->SetParameters( xmin[0], 1.0, 1.0, 1.0, 1.0 );
               double integralsigk = ftemplate->Integral(0,dist->range);
               ftemplate->SetParameters( xmin[0], 0.0, 1.0, 1.0, 1.0 );
               double integralbkgk = ftemplate->Integral(0,dist->range);
               ftemplate->SetParameters( xmin[0], kmbl, 1.0, integralsigk, integralbkgk );

               const double par [] = {xmin[0], kmbl, 1.0, integralsigk, integralbkgk};
               gLkmbl->SetPoint(k, kmbl, Min2LL(par) - minvalue);
            }


            // kmbl vs mt
            TH2D *hLmbl = new TH2D("hLmbl", "hLmbl", npnts_mt, mt_lrange, mt_rrange,
                  npnts_kmbl, kmbl_lrange, kmbl_rrange);

            if( false ){ 
               cout << "Generating 2d profile." << endl;
               for(unsigned int k=0; k <= npnts_mt; k++){
                  for(unsigned int j=0; j <= npnts_kmbl; j++){
                     double mt = hLmbl->GetXaxis()->GetBinCenter(k);
                     double kmbl = hLmbl->GetYaxis()->GetBinCenter(j);
                     const double par [] = {mt, kmbl, xmin[2], xmin[3]}; //what's going on here? Why is there xmin[2] and xmin[3] in the original code; there should only be 2 things in xmin.
                     hLmbl->SetBinContent(k, j, Min2LL(par) - minvalue);
                  }
               }
            }

            TCanvas *cLmt = new TCanvas( ("cLmt_"+dist->name).c_str(), ("cLmt_"+dist->name).c_str(), 800, 800);
            cLmt->cd();
            gLmt->SetMarkerStyle(20);
            gLmt->Draw("ACP");
            cLmt->Write("cLmt");

            TCanvas *cLkmbl = new TCanvas( ("cLkmbl_"+dist->name).c_str(), ("cLkmbl_"+dist->name).c_str(), 800, 800);
            cLkmbl->cd();
            gLkmbl->SetMarkerStyle(20);
            gLkmbl->Draw("ACP");
            cLkmbl->Write("cLkmbl");

            TCanvas *cLmbl = new TCanvas( ("cLmbl_"+dist->name).c_str(), ("cLmbl_"+dist->name).c_str(), 800, 800);
            cLmbl->cd();
            hLmbl->Draw("colz");
            cLmbl->Write("cLmbl");



         }

         delete fptr;
         delete ftemplate;

      }
   }

   fitchi2 = chi2;

   /*
   // 2 sigma contour
   unsigned int npnts_contour = 100;
   gMinuit->SetErrorDef(4); // note 4 and not 2
   double *x2s = new double[npnts_contour];
   double *y2s = new double[npnts_contour];
   gMinuit->Contour(0, 1, npnts_contour, x2s, y2s); // par1, par2, numpnts, x, y
   TGraph *gs2 = new TGraph(npnts_contour, x2s, y2s);
   gs2->SetLineColor(5);
   gs2->Draw("same C");

   // 1 sigma contour
   gMinuit->SetErrorDef(1); // note 4 and not 2
   double *xs = new double[100];
   double *ys = new double[100];
   gMinuit->Contour(0, 1, npnts_contour, xs, ys); // par1, par2, numpnts, x, y
   TGraph *gs = new TGraph(npnts_contour, xs, ys);
   gs->SetLineColor(5);
   gs->Draw("same C");
   */


   fileout->Close();
   return;

}

void Fitter::PDFReweight( vector<Event>& eventvec, int evalue ){
   cout << "Applying PDF Reweight, Eigenvalue " << evalue << endl;

   int esize = eventvec[0].pdf_weights.size();
   int nevents = eventvec.size();
   vector<float> pdfweight_sum;
   pdfweight_sum.resize( esize );

   // weight sums
   for( vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){
      for( int i=0; i < esize; i++ ){
         pdfweight_sum[i] += ev->pdf_weights[i];
      }
   }

   // apply weights
   for( vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){
      ev->weight *= ev->pdf_weights[evalue] / pdfweight_sum[evalue] * nevents;
   }

}
