#ifndef TOPMASS_H
#define TOPMASS_H

#include "Shapes.h"
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

using namespace std;

#define NMP 7
#define NJP 5
#define NTP 75

struct Dataset {

   // dataset info
   string path;
   string dir;

   // monte carlo
   double mc_xsec;
   int mc_nevts;

   vector<string> filenames;

   // constructor
   Dataset( string d="" ): dir(d) {
      mc_xsec = 1.0;
      mc_nevts = 1;
   }

};

struct Event {

   string process;

   // classify events
   string type;

   int njets;
   int nbjets;

   int nvertices;
   double weight;
   double weightcorr;

   // kinematic variables for Maos
   double mt2_210grid;
   double mt2_220grida;
   double mt2_220gridb;

   // kinematic variables
   double mt2_220;
   double mt2_210;
   double mt2_221;

   // Maos neutrinos
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

   vector<double> mbls;

   // reconstructed objects
   TLorentzVector jet1, jet2, lep1, lep2, met, met_uncl;

   bool isemu;

   // for fit
   bool fit_event;

   vector<float> pdf_weights;

   // pass various maos cuts
   //vector<bool> maoscut210;
   vector<bool> maoscut220;

   Event(){
      process = "";
      type = "";

      nvertices = 0;
      weight = 0;

      njets = 0;
      nbjets = 0;

      jet1 = TLorentzVector();
      jet2 = TLorentzVector();
      lep1 = TLorentzVector();
      lep2 = TLorentzVector();
      met = TLorentzVector();
      met_uncl = TLorentzVector();

      mt2_220 = 0;
      mt2_210 = 0;
      mt2_221 = 0;

      mt2_210grid = 0;
      mt2_220grida = 0;
      mt2_220gridb = 0;

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
   double gljf;
   double grho;

   TMatrixD Ainv_sig;
   TMatrixD Ainv_bkg;
   
   TVectorD aGPsig;
   TVectorD aGPbkg;

   TVectorD aGPsigUP;
   TVectorD aGPbkgUP;
   TVectorD aGPsigDN;
   TVectorD aGPbkgDN;

   double lbnd;
   double rbnd;

   double ptrain [NTP][NMP][NJP];

   Distribution( string n="", string t="", double n1=1.0, double n2=1.0, double lx=1.0,
         double lmt=1.0, double ljf=0.1, double rho=0.5, double lb=0, double rb=300 )
      : name(n), title(t), gnorm1(n1), gnorm2(n2), glx(lx), glmt(lmt), gljf(ljf), grho(rho), lbnd(lb), rbnd(rb) {
         activate = false;
   }

};

class Fitter{

   public:

      Fitter();
      ~Fitter();

      void InitializeDists();
      void ReadNtuple(Dataset, string, double, string, vector<Event>&, int, double, int, int, double);
      void LoadDatasets(map<string, Dataset>&);
      void ReadDatasets(map<string, Dataset>&, vector<Event>&, string, string, double, double, double, double);
      void GetVariables(vector<Event>&);
      void JShift(Event&, double=1.0);
      void JShift_test(Event&, double=1.0);
      double uncertainty(double=0.0);
      void ReweightMC(vector<Event>&, string);
      vector<int> Resample(vector<Event>&, int, bool);
      void PDFReweight(vector<Event>&, int);

      void RunMinimizer(vector<Event>&);
      void PlotResults(map< string, map<string, TH1D*> >&, string);

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;

      // diagnostics
      void DeclareHists( map< string, map<string, TH1D*> >&, map< string, map<string, TH2D*> >&, string label );
      void DeleteHists( map< string, map<string, TH1D*> >&, map< string, map<string, TH2D*> >& );
      void FillHists( map< string, map<string, TH1D*> >&,
            map< string, map<string, TH2D*> >&, vector<Event>&, bool=false );
      void PrintHists( map< string, map<string, TH1D*> >&, map< string, map<string, TH2D*> >&, string );
      vector<bool> MaosCut210( vector<Event>::iterator ev, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector& );
      vector<bool> MaosCut220( vector<Event>::iterator ev, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector& );
      void PlotTemplates( vector< map< string, map<string, TH1D*> > >& );

      void FindPTrain( map< string, map<string, TH1D*> >&, vector<Event>&, int );

      void ComputeProfile(int, vector<Event>&, string);

      double linapprox( double, double, double, double );

      map<string, Distribution> dists;

      bool compute_profile;
      double fitchi2;

      double tsig_mbl_chi2 [NMP];
      double tbkg_mbl_chi2 [NMP];
      
      int maoscuts220;
      int maoscuts210;

      static const double masspoints [NMP];
      static const double jfactpoints [NJP];

      bool compute_maos220;
      bool compute_maos210;

      bool fit_jfactor;

      double jtest;
      double gplength_jfact;
      double whyb;

      double mt_fix;

      static int clocks [100];

      bool use_data;
      int NDATA;


   private:

      double Min2LL(const double*);

      ROOT::Math::IMultiGenFunction* fFunc;

      vector<Event>* eventvec_fit;

};

#endif

