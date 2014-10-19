#ifndef TOPMASS_H
#define TOPMASS_H

#include "Shapes.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TVectorD.h"
#include "Configuration.h"

#include <vector>
#include <cmath>
#include <map>
#include <TH1.h>
#include <TH2.h>

using namespace std;

class Fitter{

   public:

      Fitter();
      ~Fitter();

      void InitializeDists();
      void ReadNtuple(string, string, double, string, vector<Event>&, int=0, int=0, double=-1, int=-1, int=-1);
      void LoadDatasets(map<string, Dataset>&);
      void GetVariables(vector<Event>&);
      void ReweightMC(vector<Event>&, string);

      void RunMinimizer(vector<Event>&, vector<Event>&, map< string, map<string, TH1D*> >&);
      void PlotResults(map< string, map<string, TH1D*> >&);

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;

      // diagnostics
      void DeclareHists( map< string, map<string, TH1D*> >&, map< string, map<string, TH2D*> >&, string label );
      void DeleteHists( map< string, map<string, TH1D*> >&, map< string, map<string, TH2D*> >& );
      void FillHists( map< string, map<string, TH1D*> >&, map< string, map<string, TH2D*> >&, vector<Event>&, bool=false );
      void PrintHists( map< string, map<string, TH1D*> >&, map< string, map<string, TH2D*> >& );
      vector<bool> MaosCut220( vector<Event>::iterator ev );
      vector<bool> MaosCut210( vector<Event>::iterator ev );
      void PlotTemplates( map< string, map<string, TH1D*> >& );

      map<string, Distribution> dists;

      bool compute_profile;
      double fitchi2;

      double tsig_mbl_chi2 [8];
      double tbkg_mbl_chi2 [8];
      
      int maoscuts220;
      int maoscuts210;

   private:

      double Min2LL(const double*);

      ROOT::Math::IMultiGenFunction* fFunc;

      vector<Event>* eventvec_fit;
      vector<Event>* eventvec_train;
      map< string, map<string, TH1D*> >* hists_fit;

};

#endif

