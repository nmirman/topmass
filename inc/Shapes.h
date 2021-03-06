#include "TH1D.h"
#include "TVectorD.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

#include <map>

#ifndef SHAPES_H
#define SHAPES_H

#define NMP 7
#define NJP 5
#define NTP 75
#define NGLB = NMP * NJP

using namespace std;

class Shapes{

   public:

      Shapes( string, double[NTP][NMP][NJP] , double, double, double, double, double, double, double, double );
      ~Shapes();

      double Ftot(double*, double*);
      double Fmbl_gp(double, double, double);
      double Fmbl_gp_var(double, double, double);

      string name;
      double lx, lmass, gnorm1, gnorm2, ljfact, grho;
      double lbx, rbx;
      double ltrain;
      double rtrain;

      //vector<double> ptrain;
      double ptrain [NTP][NMP][NJP];
      TMatrixD Ainv_sig;
      TVectorD aGPsig;
      void SetGPopts();
      void TrainGP( vector< map< string, map<string, TH1D*> > >&, double& );
      double GPkern(double, double, double, double, double, double);

      void LearnGPparams( vector< map< string, map<string, TH1D*> > >& );

      void iGP( int, int&, int& );

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;

      double masspnts [NMP];
      double jfactpnts [NJP];

      unsigned int ptrainsize;
      double ilx2, ilmass2, iljfact2, mjfcorr, irho2;

   private:
      double GPm2ll(const double*);
      double GPm2llX(const double*);
      double GPm2llLOOCV(const double*);
      ROOT::Math::IMultiGenFunction* fFunc;
      vector< map< string, map<string, TH1D*> > >* hists_train_;

      bool do_gpvar;

};

#endif

