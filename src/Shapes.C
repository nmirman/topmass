#include "TopMass.h"
#include "Shapes.h"

#include "TH1.h"
#include "TMath.h"
#include "TDecompLU.h"
#include "TDecompChol.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;


//
// constructor and destructor
//

Shapes::Shapes( string var, double ptraintmp[NTP][NMP][NJP], double gplength_x, double gplength_mt, double gplength_jfact, double norm1, double norm2, double lbound, double rbound ){

   name = var;
   // GP options
   lx = gplength_x;
   lmass = gplength_mt;
   gnorm1 = norm1;
   gnorm2 = norm2;
   ltrain = lbound;
   rtrain = rbound;
   //for(int i=0; i < ntrain; i++) ptrain.push_back( ltrain + (i+0.5)*(rtrain-ltrain)/ntrain );
   //for(unsigned int i=0; i < ptraintmp.size(); i++) ptrain.push_back( ptraintmp[i] );
   for(int i=0; i < NTP; i++){
      for(int j=0; j < NMP; j++){
         for(int k=0; k < NJP; k++){
            ptrain[i][j][k] = ptraintmp[i][j][k];
         }
      }
   }
   ptrainsize = NTP;

   ljfact = gplength_jfact;

   // right and left bounds -- set to zero unless needed
   lbx = 0.0;
   rbx = 0.0;

   double rho = 0.0;
   // JES 5 point
   //if( name == "mbl_gp" ) rho = 0.698;
   //else if( name == "mt2_221_gp" ) rho = 0.690;
   if( name == "mbl_gp" ) rho = 0.892;
   else if( name == "mt2_221_gp" ) rho = 0.884;
   else cout << "ERROR IN GP CORR: DIST NOT FOUND" << endl;
   // to avoid divisions in GPkern
   irho2 = 1.0/(1-rho*rho);
   mjfcorr = rho/(lmass*ljfact);
   ilx2 = 1.0/(lx*lx);
   ilmass2 = 1.0/(lmass*lmass);
   iljfact2 = 1.0/(ljfact*ljfact);

   // flag for cross validation two-stage fit
   do_gpvar = false;

   for(int i=0; i < NMP; i++) masspnts[i] = Fitter::masspoints[i];
   for(int i=0; i < NJP; i++) jfactpnts[i] = Fitter::jfactpoints[i];

}

Shapes::~Shapes(){
}


//
// member definitions
//

double Shapes::Ftot(double *px, double *pp){
   double x = px[0];
   double mt = pp[0];
   double jfact = pp[1];
   double norm = pp[2];
   double integralsig = pp[3];

   double val = norm*Fmbl_gp(x, mt, jfact)/integralsig;

   if( val <= 0 or (x > lbx and x < rbx) ) return 1E-10;
   else return val;
}

double Shapes::Fmbl_gp(double x, double mt, double jfact){

   double fgp = 0;
   for(unsigned int i=0; i < ptrainsize; i++){
      for(int j=0; j < NMP; j++){
         for(int k=0; k < NJP; k++){
            fgp += aGPsig[i+j*ptrainsize+k*ptrainsize*NMP]
               * GPkern( x, ptrain[i][j][k], mt, masspnts[j], jfact, jfactpnts[k] );
         }
      }
   }

   return fgp;
}

double Shapes::Fmbl_gp_var(double x, double mt, double jfact){

   int ntrain = ptrainsize;

   // vector of covariances
   TVectorD k(ntrain*NMP);
   for(int i=0; i < ntrain*NMP*NJP; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      int ijfact = i / (ntrain*NMP);
      k[i] = GPkern( x, ptrain[im][imass][ijfact], mt, masspnts[imass], jfact, jfactpnts[ijfact] );
   }
   TVectorD kT = k;

   double c1=0, c2=0;
   c1 = GPkern( x, x, mt, mt, jfact, jfact );
   k *= Ainv_sig;
   c2 = kT*k;

   return (c1-c2);
}

double Shapes::GPkern(double x1, double x2, double m1, double m2, double j1, double j2){

   return 1E-06*gnorm2*gnorm1*exp( -0.5*( ilx2*(x1-x2)*(x1-x2) + ilmass2*irho2*(m1-m2)*(m1-m2)
            + iljfact2*irho2*(j1-j2)*(j1-j2) - 2*mjfcorr*irho2*(x1-x2)*(j1-j2) ));
}

void Shapes::TrainGP( vector< map< string, map<string, TH1D*> > >& hists_,
     double &m2llsig ){

   int ntrain = ptrainsize;

   // histograms
   vector<TH1D*> hgp_sig;
   for(int i=0; i < NMP*NJP; i++){

      int im, ijfact;
      iGP(i,im,ijfact);

      stringstream ssmass;
      ssmass << floor(masspnts[im]);
      string smass = ssmass.str();

      stringstream sjfact;
      sjfact << jfactpnts[ijfact];
      string jfact = sjfact.str();

      // signal shape
      hgp_sig.push_back( (TH1D*)hists_[ijfact][name]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp_sig"+smass).c_str()) );
      hgp_sig[i]->Add( hists_[ijfact][name]["ttbar"+smass+"_mistag"] );
      hgp_sig[i]->Add( hists_[ijfact][name]["ttbar"+smass+"_hadronic"] );
      hgp_sig[i]->Add( hists_[ijfact][name]["ttbar"+smass+"_taus"] );
      hgp_sig[i]->Add( hists_[ijfact][name]["other"] );
      hgp_sig[i]->Scale( 1.0/hgp_sig[i]->Integral("width") );

      for(int n=0; n < hgp_sig[i]->GetNbinsX(); n++){
         if( hgp_sig[i]->GetBinError(n) < 5E-06 ) hgp_sig[i]->SetBinError(n, 5E-06);
      }

   }

   // compute covariance matrix
   TMatrixD K(ntrain*NMP*NJP,ntrain*NMP*NJP);
   for(int i=0; i < ntrain*NMP*NJP; i++){
      for(int j=0; j < ntrain*NMP*NJP; j++){
         int im = i % ntrain;
         int jm = j % ntrain;

         int i_index = i / ntrain;
         int j_index = j / ntrain;
         int imass, jmass, ijfact, jjfact;
         iGP( i_index, imass, ijfact );
         iGP( j_index, jmass, jjfact );

         K[i][j] = GPkern( ptrain[im][imass][ijfact], ptrain[jm][jmass][jjfact],
               masspnts[imass], masspnts[jmass], jfactpnts[ijfact], jfactpnts[jjfact] );
     }
   }
   // compute noise matrix
   TMatrixD Nsig(ntrain*NMP*NJP,ntrain*NMP*NJP);
   for(int i=0; i < ntrain*NMP*NJP; i++){
      int im = i % ntrain;
      int index = i / ntrain;
      int imass, ijfact;
      iGP( index, imass, ijfact );
      double binerr_sig = hgp_sig[index]->GetBinError( hgp_sig[index]->FindBin(ptrain[im][imass][ijfact]) );
      binerr_sig *= sqrt(gnorm2);
      for(int j=0; j < ntrain*NMP*NJP; j++){
         if( i==j ){
            Nsig[i][j] = binerr_sig*binerr_sig;//pow( max(binerr_sig,0.001), 2 );
         }else{
            Nsig[i][j] = 0;
         }
     }
   }

   // inverse of sum
   TMatrixD Asig = K + Nsig;
   cout << "---> cholesky decomposition of signal matrix... "; fflush(stdout);
   TDecompChol Cholsig(Asig);
   cout << "done!" << endl;
   TMatrixD AsigU = Cholsig.GetU();
   bool status = 0;
   cout << "---> invert signal matrix... "; fflush(stdout);
   TMatrixDSym Asinv_sig = Cholsig.Invert(status);
   cout << "done!" << endl;
   Ainv_sig.Clear();
   Ainv_sig.ResizeTo( ntrain*NMP*NJP, ntrain*NMP*NJP );
   Ainv_sig = (TMatrixD)Asinv_sig;

   TMatrixD Ainv_sigtemp = Ainv_sig;

   // vector of training points
   TVectorD ysig(ntrain*NMP*NJP);
   for(int i=0; i < ntrain*NMP*NJP; i++){
      int im = i % ntrain;
      int index = i / ntrain;
      int imass, ijfact;
      iGP( index, imass, ijfact );
      ysig[i] = hgp_sig[index]->GetBinContent( hgp_sig[index]->FindBin(ptrain[im][imass][ijfact]) );
   }

   // alpha vector
   aGPsig.Clear();

   aGPsig.ResizeTo( ntrain*NMP*NJP );
   aGPsig = Ainv_sigtemp*ysig;

   // compute marginal likelihood
   double ldetsig = 0.0;
   for(int i=0; i < Cholsig.GetNrows(); i++){
      ldetsig += 2*log(AsigU[i][i]);
   }

   double term1sig = -0.5*ysig*aGPsig;
   double term2sig = -0.5*ldetsig;
   double term3sig = -0.5*ntrain*log(2*TMath::Pi());

   m2llsig = -2.0*(term1sig+term2sig+term3sig);

   return;
}

void Shapes::LearnGPparams( vector< map< string, map<string, TH1D*> > >& hists_ ){

   gMinuit = new ROOT::Minuit2::Minuit2Minimizer ( ROOT::Minuit2::kMigrad );
   gMinuit->SetPrintLevel(3);

   // set training hist
   hists_train_ = &hists_;

   //fFunc = new ROOT::Math::Functor ( this, &Shapes::GPm2llX, 5 );
   fFunc = new ROOT::Math::Functor ( this, &Shapes::GPm2llLOOCV, 6 );
   gMinuit->SetFunction( *fFunc );

   /*
   do_gpvar = false;
   gMinuit->SetLowerLimitedVariable(0, "gpnorm1", 0.5, 1, 0.0);
   gMinuit->SetLowerLimitedVariable(1, "gpnorm2", 30, 1, 0.0);
   gMinuit->SetLowerLimitedVariable(2, "lx", 20, 1, 0.0);
   gMinuit->SetLowerLimitedVariable(3, "lmass", 30, 1, 0.0);
   
   gMinuit->SetTolerance(1);
   gMinuit->Minimize();

   const double *xs = gMinuit->X();
   gnorm1 = xs[0];
   gnorm2 = xs[1];
   lx = xs[2];
   lmass = xs[3];
   */
   /*
   do_gpvar = false;
   // stage 1
   gMinuit->SetFixedVariable(0, "gpnorm1", 1);
   gMinuit->SetFixedVariable(1, "gpnorm2", 1);
   gMinuit->SetLowerLimitedVariable(2, "lx", 10, 1, 0.0);
   gMinuit->SetLowerLimitedVariable(3, "lmass", 10, 1, 0.0);

   gMinuit->SetTolerance(1);
   gMinuit->Minimize();

   const double *xs = gMinuit->X();
   gnorm1 = xs[0];
   gnorm2 = xs[1];
   lx = xs[2];
   lmass = xs[3];

   // stage 2
   gMinuit->SetLowerLimitedVariable(0, "gpnorm1", 1, 0.1, 0.0);
   gMinuit->SetLowerLimitedVariable(1, "gpnorm2", 10, 1, 0.0);
   gMinuit->SetFixedVariable(2, "lx", lx);
   gMinuit->SetFixedVariable(3, "lmass", lmass);

   gMinuit->SetTolerance(1);
   gMinuit->Minimize();

   const double *xs2 = gMinuit->X();
   gnorm1 = xs2[0];
   gnorm2 = xs2[1];
   lx = xs2[2];
   lmass = xs2[3];
   */

   // name(n), title(t), gnorm1(n1), gnorm2(n2), glx(lx), glmt(lmt), gljf(lfj), lbnd(lb), rbnd(rb)
   // name(n), title(t), theta1, theta0, theta2, theta3
   //dists[ "mbl_gp" ] = Distribution( "mbl_gp", "M_{bl}", 0.54, 2.1, 7.2, 8.1, ljf, 20, 300 );

   do_gpvar = false;
   gMinuit->SetFixedVariable(0, "gpnorm1", gnorm1);
   gMinuit->SetFixedVariable(1, "gpnorm2", gnorm2);
   gMinuit->SetFixedVariable(2, "lx", sqrt(1.0/ilx2));
   gMinuit->SetFixedVariable(3, "lmass", sqrt(1.0/ilmass2));
   gMinuit->SetVariable(4, "ljf", sqrt(1.0/iljfact2), 0.1);
   gMinuit->SetLimitedVariable(5, "rho", mjfcorr/sqrt(iljfact2*ilmass2), 0.1, -1 , 1);

   cout << "############ STAGE 1 FIT ###########" << endl;
   gMinuit->Minimize();
   const double *xs = gMinuit->X();
   gnorm1 = xs[0];
   gnorm2 = xs[1];
   double rho = xs[5];
   irho2 = 1.0/(1-rho*rho);
   ilx2 = 1.0/(xs[2]*xs[2]);
   ilmass2 = 1.0/(xs[3]*xs[3]);
   iljfact2 = 1.0/(xs[4]*xs[4]);
   mjfcorr = sqrt(iljfact2*ilmass2)*rho;

   gMinuit->SetFixedVariable(0, "gpnorm1", gnorm1);
   gMinuit->SetFixedVariable(1, "gpnorm2", gnorm2);
   gMinuit->SetVariable(2, "lx", sqrt(1.0/ilx2), 1.0);
   //gMinuit->SetVariable(3, "lmass", sqrt(1.0/ilmass2), 1.0);
   gMinuit->SetFixedVariable(3, "lmass", sqrt(1.0/ilmass2));
   gMinuit->SetFixedVariable(4, "ljf", sqrt(1.0/iljfact2));
   gMinuit->SetFixedVariable(5, "rho", mjfcorr/sqrt(iljfact2*ilmass2));

   cout << "############ STAGE 2 FIT ###########" << endl;
   gMinuit->Minimize();
   const double *xs2 = gMinuit->X();
   gnorm1 = xs2[0];
   gnorm2 = xs2[1];
   rho = xs2[5];
   irho2 = 1.0/(1-rho*rho);
   ilx2 = 1.0/(xs2[2]*xs2[2]);
   ilmass2 = 1.0/(xs2[3]*xs2[3]);
   iljfact2 = 1.0/(xs2[4]*xs2[4]);
   mjfcorr = sqrt(iljfact2*ilmass2)*rho;

   return;
}

double Shapes::GPm2ll( const double *x ){
   cout << "gnorm: " << x[0] << ", " << x[1] << endl;
   cout << "lx, lmt, ljf: " << x[1] << ", " << x[2] << ", " << x[3] << endl;

   gnorm1 = x[0];
   gnorm2 = x[1];
   lx = x[2];
   lmass = x[3];
   double m2llsig;
   TrainGP( *hists_train_, m2llsig );

   return m2llsig;
}

double Shapes::GPm2llX( const double *x ){
   /*
   cout << "gnorm: " << x[0] << ", " << x[1] << endl;
   cout << "lx, lmt, ljf: " << x[2] << ", " << x[3] << ", " << x[4] << endl;

   gnorm1 = x[0];
   gnorm2 = x[1];
   lx = x[2];
   lmass = x[3];
   ljfact = x[4];

   // histograms
   vector<TH1D*> hgp_sig;
   for(int j=0; j < NMP; j++){

      stringstream ssmass;
      ssmass << floor(masspnts[j]);
      string smass = ssmass.str();

      // signal shape
      hgp_sig.push_back( (TH1D*)(*hists_train_)[1][name]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp_sigx"+smass).c_str()) );
      hgp_sig[j]->Add( (*hists_train_)[1][name]["ttbar"+smass+"_mistag"] );
      hgp_sig[j]->Add( (*hists_train_)[1][name]["ttbar"+smass+"_hadronic"] );
      hgp_sig[j]->Add( (*hists_train_)[1][name]["ttbar"+smass+"_taus"] );
      hgp_sig[j]->Add( (*hists_train_)[1][name]["other"] );
      hgp_sig[j]->Scale( 1.0/hgp_sig[j]->Integral("width") );

   }

   int nval = 3;
   double m2ll_tot = 0;
   for(int c=0; c < nval; c++){ // n-fold cross validation

      ptrain.clear();
      vector<double> ptrainX;
      for(int i=1; i <= hgp_sig[0]->GetNbinsX(); i++){ // exclude every nth point
         // check for zero bins
         bool zerobin = false;
         for(int j=0; j < NMP; j++){
            double val = hgp_sig[j]->GetBinContent(i);
            if( val == 0 ) zerobin = true;
         }
         if( i%nval != c or zerobin ){
            ptrain.push_back( hgp_sig[0]->GetBinCenter(i) );
         }else{
            ptrainX.push_back( hgp_sig[0]->GetBinCenter(i) );
         }
      }
      double m2llsig, m2llbkg;
      TrainGP( *hists_train_, m2llsig, m2llbkg );

      // evaluate shape at excluded points
      double m2ll = 0;
      for(unsigned int i=0; i < ptrainX.size(); i++){
         for(int j=0; j < NMP; j++){

            double mean = Fmbl_gp(ptrainX[i],masspnts[j],1.0,"sig");
            double yi = hgp_sig[j]->GetBinContent( hgp_sig[j]->FindBin(ptrainX[i]) );
            double var = -1;
            if( do_gpvar ) var = Fmbl_gp_var(ptrainX[i],masspnts[j],1.0,"sig");
            else var = pow(hgp_sig[j]->GetBinError( hgp_sig[j]->FindBin(ptrainX[i]) ), 2);

            if( var <= 0 ){
               //cout << "NEGATIVE VARIANCE IN GP: " << var << "  ----> setting to minimum" << endl;
               var = 1E-10;
            }

            m2ll += log(var) + pow(mean-yi,2)/var + log(2*TMath::Pi());

         }
      }

      m2ll_tot += m2ll;
   }

   return m2ll_tot;

   */
      return x[0]; // temporary placeholder
}

double Shapes::GPm2llLOOCV( const double *x ){
   cout << "gnorm: " << x[0] << ", " << x[1] << endl;
   cout << "lx, lmt, ljf, corr: " << x[2] << ", " << x[3] << ", " << x[4] << ", " << x[5] << endl;

   double rho = x[5];
   irho2 = 1.0/(1-rho*rho);
   gnorm1 = x[0];
   gnorm2 = x[1];
   ilx2 = 1.0/(x[2]*x[2]);
   ilmass2 = 1.0/(x[3]*x[3]);
   iljfact2 = 1.0/(x[4]*x[4]);
   mjfcorr = sqrt(iljfact2*ilmass2)*rho;

   int ntrain = ptrainsize;

   // histograms
   vector<TH1D*> hgp_sig;
   for(int i=0; i < NMP*NJP; i++){

      int im, ijfact;
      iGP(i,im,ijfact);

      stringstream ssmass;
      ssmass << floor(masspnts[im]);
      string smass = ssmass.str();

      stringstream sjfact;
      sjfact << jfactpnts[ijfact];
      string jfact = sjfact.str();

      // signal shape
      hgp_sig.push_back( (TH1D*)(*hists_train_)[ijfact][name]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp_sigx"+smass).c_str()) );
      hgp_sig[i]->Add( (*hists_train_)[ijfact][name]["ttbar"+smass+"_mistag"] );
      hgp_sig[i]->Add( (*hists_train_)[ijfact][name]["ttbar"+smass+"_hadronic"] );
      hgp_sig[i]->Add( (*hists_train_)[ijfact][name]["ttbar"+smass+"_taus"] );
      hgp_sig[i]->Add( (*hists_train_)[ijfact][name]["other"] );
      hgp_sig[i]->Scale( 1.0/hgp_sig[i]->Integral("width") );

   }

   // precompute Ky vector
   double m2llsig;
   TrainGP( *hists_train_, m2llsig );

   // vector of training points
   TVectorD ysig(ntrain*NMP*NJP);
   for(int i=0; i < ntrain*NMP*NJP; i++){
      int im = i % ntrain;
      int index = i / ntrain;
      int imass, ijfact;
      iGP( index, imass, ijfact );
      ysig[i] = hgp_sig[index]->GetBinContent( hgp_sig[index]->FindBin(ptrain[im][imass][ijfact]) );
   }

   TVectorD Ky = Ainv_sig*ysig;

   double m2ll = 0;
   for(int i=0; i < ntrain*NMP*NJP; i++){
      if( ysig[i] == 0 ) continue;
      int im = i % ntrain;
      int index = i / ntrain;
      int imass, ijfact;
      iGP( index, imass, ijfact );

      double ui = ysig[i] - Ky[i]/Ainv_sig[i][i];
      double vi1 = 1.0/Ainv_sig[i][i];
      double vi2 = pow(hgp_sig[index]->GetBinContent( hgp_sig[index]->FindBin(ptrain[im][imass][ijfact]) ),2);

      double vi = do_gpvar ? vi1 : vi2;

      m2ll += log(vi) + pow(ui-ysig[i],2)/vi + log(2*TMath::Pi());
   }
   cout << "m2ll = " << m2ll << endl;

   return m2ll;

}

void Shapes::iGP( int index, int &mass, int &jfactor ){

   mass = index % NMP;
   jfactor = index / NMP;
   
   return;
}
