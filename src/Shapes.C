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

Shapes::Shapes( string var, vector<double>& ptraintmp, double gplength_x, double gplength_mt, double norm1, double norm2, double lbound, double rbound ){

   name = var;
   // GP options
   lx = gplength_x;
   lmass = gplength_mt;
   gnorm1 = norm1;
   gnorm2 = norm2;
   //ntrain = 100;
   ltrain = lbound;
   rtrain = rbound;
   //for(int i=0; i < ntrain; i++) ptrain.push_back( ltrain + (i+0.5)*(rtrain-ltrain)/ntrain );
   for(unsigned int i=0; i < ptraintmp.size(); i++) ptrain.push_back( ptraintmp[i] );

   // right and left bounds -- set to zero unless needed
   lbx = 0.0;
   rbx = 0.0;

   // flag for cross validation two-stage fit
   do_gpvar = false;

   for(int i=0; i < NMP; i++) masspnts[i] = Fitter::masspoints[i];

}

Shapes::~Shapes(){
}


//
// member definitions
//

double Shapes::Ftot(double *px, double *pp){

   double x = px[0];
   double mt = pp[0];
   double k = pp[1];
   double norm = pp[2];
   double integralsig = pp[3];
   double integralbkg = pp[4];

   double val = norm*(k*Fmbl_gp(x, mt, "sig")/integralsig + (1-k)*Fmbl_gp(x, mt, "bkg")/integralbkg);
   if( val <= 0 or (x > lbx and x < rbx) ) return 1E-10;
   else return val;
}

double Shapes::Fsig_param(double x, double mt){

   double par [9];
   par[0] = 0.04314 - 0.0001872*mt;
   par[1] = -89.26 + 1.139*mt;
   par[2] = 20.86 + 0.004963*mt;
   par[3] = -0.005268 + 7.442E-05*mt;
   par[4] = -90.66 + 0.9241*mt;
   par[5] = -26.65 + 0.2591*mt;
   par[6] = 0.02192 + -6.172E-05*mt;
   par[7] = -23.11 + 0.3386*mt;
   par[8] = -15.22 + 0.1252*mt;

   double gaus1 = par[0]*exp( -0.5*pow((x-par[1])/par[2],2) );
   double gaus2 = par[3]*exp( -0.5*pow((x-par[4])/par[5],2) );
   double landau = par[6]*TMath::Landau(x,par[7],par[8],false);

   return gaus1+gaus2+landau;
}

double Shapes::Fmbl_gp(double x, double mt, string sb){

   double fgp = 0;
   for(unsigned int i=0; i < ptrain.size(); i++){
      for(int j=0; j < NMP; j++){
         double agp = 0;
         if( sb.compare("sig") == 0 ) agp = aGPsig[i+j*ptrain.size()];
         else if( sb.compare("bkg") == 0 ) agp = aGPbkg[i+j*ptrain.size()];
         else{
            cout << "ERROR in GP shape." << endl;
            return -1;
         }
         fgp += agp*GPkern( x, ptrain[i], lx, mt, masspnts[j], lmass );
      }
   }

   return fgp;
}

double Shapes::Fmbl_gp_var(double x, double mt, string sb){

   int ntrain = ptrain.size();

   // vector of covariances
   TVectorD k(ntrain*NMP);
   for(int i=0; i < ntrain*NMP; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      k[i] = GPkern( x, ptrain[im], lx, mt, masspnts[imass], lmass );
   }
   TVectorD kT = k;

   double c1=0, c2=0;
   c1 = GPkern( x, x, lx, mt, mt, lmass );
   if( sb.compare("sig") == 0 ) k *= Ainv_sig;
   else if( sb.compare("bkg") == 0 ) k *= Ainv_bkg;
   else{
      cout << "ERROR in GP shape." << endl;
      return -1;
   }
   c2 = kT*k;

   return (c1-c2);
}

double Shapes::GPkern(double x1, double x2, double lsx, double m1, double m2, double lsm ){

   double kernel = 1E-06*gnorm2*gnorm1*exp(-0.5*(pow( (x1-x2)/lsx, 2)+pow( (m1-m2)/lsm, 2)));
   return kernel;
}

void Shapes::TrainGP( map< string, map<string, TH1D*> > & hists_,
     double &m2llsig, double &m2llbkg ){

   int ntrain = ptrain.size();

   // histograms
   vector<TH1D*> hgp_sig;
   vector<TH1D*> hgp_bkg;
   for(int i=0; i < NMP; i++){

      stringstream ssmass;
      ssmass << floor(masspnts[i]);
      string smass = ssmass.str();

      // signal shape
      hgp_sig.push_back( (TH1D*)hists_[name]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp_sig"+smass).c_str()) );
      hgp_sig[i]->Scale( 1.0/hgp_sig[i]->Integral("width") );

      // background shape
      hgp_bkg.push_back( (TH1D*)hists_[name]["ttbar"+smass+"_mistag"]
            ->Clone( ("hgp_bkg"+smass).c_str()) );
      hgp_bkg[i]->Add( hists_[name]["ttbar"+smass+"_hadronic"] );
      hgp_bkg[i]->Add( hists_[name]["ttbar"+smass+"_taus"] );
      hgp_bkg[i]->Add( hists_[name]["other"] );
      hgp_bkg[i]->Scale( 1.0/hgp_bkg[i]->Integral("width") );

      for(int n=0; n < hgp_sig[i]->GetNbinsX(); n++){
         if( hgp_sig[i]->GetBinError(n) < 5E-06 ) hgp_sig[i]->SetBinError(n, 5E-06);
         if( hgp_bkg[i]->GetBinError(n) < 5E-06 ) hgp_bkg[i]->SetBinError(n, 5E-06);
      }
      

   }

   // compute covariance matrix
   TMatrixD K(ntrain*NMP,ntrain*NMP);
   for(int i=0; i < ntrain*NMP; i++){
      for(int j=0; j < ntrain*NMP; j++){
         int im = i % ntrain;
         int jm = j % ntrain;
         int imass = i / ntrain;
         int jmass = j / ntrain;
         K[i][j] = GPkern( ptrain[im], ptrain[jm], lx, masspnts[imass], masspnts[jmass], lmass );
     }
   }
   // compute noise matrix
   TMatrixD Nsig(ntrain*NMP,ntrain*NMP);
   TMatrixD Nbkg(ntrain*NMP,ntrain*NMP);
   for(int i=0; i < ntrain*NMP; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      double binerr_sig = hgp_sig[imass]->GetBinError( hgp_sig[imass]->FindBin(ptrain[im]) );
      double binerr_bkg = hgp_bkg[imass]->GetBinError( hgp_bkg[imass]->FindBin(ptrain[im]) );
      binerr_sig *= sqrt(gnorm2);
      binerr_bkg *= sqrt(gnorm2);
      for(int j=0; j < ntrain*NMP; j++){
         if( i==j ){
            Nsig[i][j] = binerr_sig*binerr_sig;//pow( max(binerr_sig,0.001), 2 );
            Nbkg[i][j] = binerr_bkg*binerr_bkg;//pow( max(binerr_bkg,0.001), 2 );
         }else{
            Nsig[i][j] = 0;
            Nbkg[i][j] = 0;
         }
     }
   }

   // inverse of sum
   TMatrixD Asig = K + Nsig;
   TMatrixD Abkg = K + Nbkg;
   TDecompChol Cholsig(Asig);
   TDecompChol Cholbkg(Abkg);
   TMatrixD AsigU = Cholsig.GetU();
   TMatrixD AbkgU = Cholbkg.GetU();
   bool status = 0;
   TMatrixDSym Asinv_sig = Cholsig.Invert(status);
   TMatrixDSym Asinv_bkg = Cholbkg.Invert(status);
   Ainv_sig.Clear();
   Ainv_bkg.Clear();
   Ainv_sig.ResizeTo( ntrain*NMP, ntrain*NMP );
   Ainv_bkg.ResizeTo( ntrain*NMP, ntrain*NMP );
   Ainv_sig = (TMatrixD)Asinv_sig;
   Ainv_bkg = (TMatrixD)Asinv_bkg;

   TMatrixD Ktmp = K;
   for(int i=0; i < ntrain*NMP; i++) Ktmp[i][i] += 10E-9;
   TDecompChol CholK(Ktmp);
   status = 0;
   TMatrixDSym Ksinv = CholK.Invert(status);
   Kinv.Clear();
   Kinv.ResizeTo( ntrain*NMP, ntrain*NMP );
   Kinv = (TMatrixD)Ksinv;

   TMatrixD Ainv_sigtemp = Ainv_sig;
   TMatrixD Ainv_bkgtemp = Ainv_bkg;

   // vector of training points
   TVectorD ysig(ntrain*NMP);
   TVectorD ybkg(ntrain*NMP);
   for(int i=0; i < ntrain*NMP; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      ysig[i] = hgp_sig[imass]->GetBinContent( hgp_sig[imass]->FindBin(ptrain[im]) );
      ybkg[i] = hgp_bkg[imass]->GetBinContent( hgp_bkg[imass]->FindBin(ptrain[im]) );
   }

   // alpha vector
   aGPsig.Clear();
   aGPbkg.Clear();

   aGPsig.ResizeTo( ntrain*NMP );
   aGPsig = Ainv_sigtemp*ysig;

   aGPbkg.ResizeTo( ntrain*NMP );
   aGPbkg = Ainv_bkgtemp*ybkg;

   // compute marginal likelihood
   //TDecompLU lusig(Asig);
   //TDecompLU lubkg(Abkg);
   //TMatrixD AsigLU = lusig.GetLU();
   //TMatrixD AbkgLU = lubkg.GetLU();

   double ldetsig = 0.0;
   double ldetbkg = 0.0;
   for(int i=0; i < Cholsig.GetNrows(); i++){
      ldetsig += 2*log(AsigU[i][i]);
      ldetbkg += 2*log(AbkgU[i][i]);
   }

   double term1sig = -0.5*ysig*aGPsig;
   double term2sig = -0.5*ldetsig;
   double term3sig = -0.5*ntrain*log(2*TMath::Pi());

   double term1bkg = -0.5*ybkg*aGPbkg;
   double term2bkg = -0.5*ldetbkg;
   double term3bkg = -0.5*ntrain*log(2*TMath::Pi());

   m2llsig = -2.0*(term1sig+term2sig+term3sig);
   m2llbkg = -2.0*(term1bkg+term2bkg+term3bkg);

   return;
}

void Shapes::LearnGPparams( map< string, map<string, TH1D*> > & hists_ ){

   gMinuit = new ROOT::Minuit2::Minuit2Minimizer ( ROOT::Minuit2::kMigrad );
   gMinuit->SetPrintLevel(3);

   // set training hist
   hists_train_ = &hists_;

   fFunc = new ROOT::Math::Functor ( this, &Shapes::GPm2llX, 4 );
   //fFunc = new ROOT::Math::Functor ( this, &Shapes::GPm2llLOOCV, 4 );
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

   return;
}

double Shapes::GPm2ll( const double *x ){
   cout << "gnorm: " << x[0] << ", " << x[1] << endl;
   cout << "lx, lmt: " << x[1] << ", " << x[2] << endl;

   gnorm1 = x[0];
   gnorm2 = x[1];
   lx = x[2];
   lmass = x[3];
   double m2llsig, m2llbkg;
   TrainGP( *hists_train_, m2llsig, m2llbkg );

   return m2llsig;
}

double Shapes::GPm2llX( const double *x ){
   cout << "gnorm: " << x[0] << ", " << x[1] << endl;
   cout << "lx, lmt: " << x[2] << ", " << x[3] << endl;

   gnorm1 = x[0];
   gnorm2 = x[1];
   lx = x[2];
   lmass = x[3];

   // histograms
   vector<TH1D*> hgp_sig;
   for(int j=0; j < NMP; j++){

      stringstream ssmass;
      ssmass << floor(masspnts[j]);
      string smass = ssmass.str();

      // signal shape
      hgp_sig.push_back( (TH1D*)(*hists_train_)[name]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp_sigx"+smass).c_str()) );
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

            double mean = Fmbl_gp(ptrainX[i],masspnts[j],"sig");
            double yi = hgp_sig[j]->GetBinContent( hgp_sig[j]->FindBin(ptrainX[i]) );
            double var = -1;
            if( do_gpvar ) var = Fmbl_gp_var(ptrainX[i],masspnts[j],"sig");
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

}

double Shapes::GPm2llLOOCV( const double *x ){
   cout << "gnorm: " << x[0] << ", " << x[1] << endl;
   cout << "lx, lmt: " << x[2] << ", " << x[3] << endl;

   gnorm1 = x[0];
   gnorm2 = x[1];
   lx = x[2];
   lmass = x[3];

   int ntrain = 100;

   // histograms
   vector<TH1D*> hgp_sig;
   for(int j=0; j < NMP; j++){

      stringstream ssmass;
      ssmass << floor(masspnts[j]);
      string smass = ssmass.str();

      // signal shape
      hgp_sig.push_back( (TH1D*)(*hists_train_)[name]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp_sigx"+smass).c_str()) );
      hgp_sig[j]->Scale( 1.0/hgp_sig[j]->Integral("width") );

   }

   // precompute Ky vector
   double m2llsig, m2llbkg;
   TrainGP( *hists_train_, m2llsig, m2llbkg );

   // vector of training points
   TVectorD ysig(ntrain*NMP);
   for(int i=0; i < ntrain*NMP; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      ysig[i] = hgp_sig[imass]->GetBinContent( hgp_sig[imass]->FindBin(ptrain[im]) );
   }

   TVectorD Ky = Kinv*ysig;

   double m2ll = 0;
   for(int i=0; i < ntrain*NMP; i++){
      if( ysig[i] == 0 ) continue;
      int im = i % ntrain;
      int imass = i / ntrain;

      double ui = ysig[i] - Ky[i]/Kinv[i][i];
      double vi1 = 1/Kinv[i][i];
      double vi2 = pow(hgp_sig[imass]->GetBinContent( hgp_sig[imass]->FindBin(ptrain[im]) ),2);

      double vi = do_gpvar ? vi1 : vi2;

      m2ll += log(vi) + pow(ui-ysig[i],2)/vi + log(2*TMath::Pi());

   }

   return m2ll;

}
