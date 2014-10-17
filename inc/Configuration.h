#ifndef CONFIG_H
#define CONFIG_H

#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TVectorD.h"

#include <vector>
#include <cmath>
#include <map>
#include <TH1.h>
#include <TH2.h>
#include <string>

using namespace std;

struct Dataset {

   // dataset info
   string path;
   string file;

   // monte carlo
   double mc_xsec;
   int mc_nevts;

   // constructor
   Dataset( string p="", string f="" ): path(p), file(f) {
      mc_xsec = 1.0;
      mc_nevts = 1;
   }

};

struct Event {

   string process;

   // classify events
   string type;

   int nvertices;
   double weight;

   // kinematic variables
   double mt2_220;
   double mt2_210;
   double mt2_221;
   vector<double> mbls;

   // kinematic variables for Maos
   double mt2_210grid;
   double mt2_220grida;
   double mt2_220gridb;

   // Maos neutrinos
   TLorentzVector maos210_neutrino1p;
   TLorentzVector maos210_neutrino1m;
   TLorentzVector maos210_neutrino2p;
   TLorentzVector maos210_neutrino2m;

   TLorentzVector maos220_neutrino1ap;
   TLorentzVector maos220_neutrino1am;
   TLorentzVector maos220_neutrino2ap;
   TLorentzVector maos220_neutrino2am;
   TLorentzVector maos220_neutrino1bp;
   TLorentzVector maos220_neutrino1bm;
   TLorentzVector maos220_neutrino2bp;
   TLorentzVector maos220_neutrino2bm;

   double maos210_blvmass1ap;
   double maos210_blvmass1am;
   double maos210_blvmass2ap;
   double maos210_blvmass2am;
   double maos210_blvmass1bp;
   double maos210_blvmass1bm;
   double maos210_blvmass2bp;
   double maos210_blvmass2bm;

   double maos220_blvmass1ap;
   double maos220_blvmass1am;
   double maos220_blvmass2ap;
   double maos220_blvmass2am;
   double maos220_blvmass1bp;
   double maos220_blvmass1bm;
   double maos220_blvmass2bp;
   double maos220_blvmass2bm;

   // reconstructed objects
   TLorentzVector jet1, jet2, lep1, lep2, met;

   // for fit
   bool fit_event;

   Event(){
      process = "";
      type = "";

      nvertices = 0;
      weight = 0;

      jet1 = TLorentzVector();
      jet2 = TLorentzVector();
      lep1 = TLorentzVector();
      lep2 = TLorentzVector();
      met = TLorentzVector();

      mt2_220 = 0;
      mt2_210 = 0;
      mt2_221 = 0;

      mt2_210grid = 0;
      mt2_220grida = 0;
      mt2_220gridb = 0;

      maos210_neutrino1p = TLorentzVector();
      maos210_neutrino1m = TLorentzVector();
      maos210_neutrino2p = TLorentzVector();
      maos210_neutrino2m = TLorentzVector();

      maos220_neutrino1ap = TLorentzVector();
      maos220_neutrino1am = TLorentzVector();
      maos220_neutrino2ap = TLorentzVector();
      maos220_neutrino2am = TLorentzVector();
      maos220_neutrino1bp = TLorentzVector();
      maos220_neutrino1bm = TLorentzVector();
      maos220_neutrino2bp = TLorentzVector();
      maos220_neutrino2bm = TLorentzVector();

      maos210_blvmass1ap = 0;
      maos210_blvmass1am = 0;
      maos210_blvmass2ap = 0;
      maos210_blvmass2am = 0;
      maos210_blvmass1bp = 0;
      maos210_blvmass1bm = 0;
      maos210_blvmass2bp = 0;
      maos210_blvmass2bm = 0;

      maos220_blvmass1ap = 0;
      maos220_blvmass1am = 0;
      maos220_blvmass2ap = 0;
      maos220_blvmass2am = 0;
      maos220_blvmass1bp = 0;
      maos220_blvmass1bm = 0;
      maos220_blvmass2bp = 0;
      maos220_blvmass2bm = 0;

      fit_event = false;
   }

};

struct Distribution {

   string name;
   string title;
   bool activate;

   double gnorm1;
   double gnorm2;
   double glx;
   double glmt;

   TMatrixD Ainv_sig;
   TMatrixD Ainv_bkg;
   
   TVectorD aGPsig;
   TVectorD aGPbkg;

   double lbnd;
   double rbnd;
   double range;

   Distribution( string n="", string t="", double n1=1.0, double n2=1.0, double lx=1.0, double lmt=1.0, double r=300 )
      : name(n), title(t), gnorm1(n1), gnorm2(n2), glx(lx), glmt(lmt), range(r) {
         lbnd = 0;
         rbnd = 0;
         activate = false;
   }

};

#endif