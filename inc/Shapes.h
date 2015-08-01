#include "TH1D.h"
#include "TVectorD.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

#include <map>

#ifndef SHAPES_H
#define SHAPES_H

#define NMP 7
#define NJP 5
#define NGLB = NMP * NJP

using namespace std;

class Shapes{

   public:

      Shapes( string, vector<double>&, double, double, double, double, double, double );
      ~Shapes();

      double Ftot(double*, double*);
      double Fsig_param(double, double);
      double Fmbl_gp(double, double, double, string);
      double Fmbl_gp_var(double, double, double, string);

      string name;
      double lx, lmass, gnorm1, gnorm2, ljfact;
      double lbx, rbx;
      double ltrain;
      double rtrain;

      vector<double> ptrain;
      TMatrixD Ainv_sig;
      TMatrixD Ainv_bkg;
      TVectorD aGPsig;
      TVectorD aGPbkg;
      TMatrixD Kinv;
      void SetGPopts();
      void TrainGP( vector< map< string, map<string, TH1D*> > >&, double&, double& );
      double GPkern(double, double, double, double, double, double, double, double, double);

      void LearnGPparams( vector< map< string, map<string, TH1D*> > >& );

      void iGP( int, int&, int& );

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;

      double masspnts [NMP];
      double jfactpnts [NJP];

   private:
      double GPm2ll(const double*);
      double GPm2llX(const double*);
      double GPm2llLOOCV(const double*);
      ROOT::Math::IMultiGenFunction* fFunc;
      vector< map< string, map<string, TH1D*> > >* hists_train_;

      bool do_gpvar;

};

#endif

