#include "TopMass.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <string>
#include <getopt.h>
#include <unistd.h>

using namespace std;

void print_usage(){

   cout << "\nUsage: ./DoFit <flags>\n";
   cout << "Flags: \n";
   cout << setiosflags(ios::left);
   cout << setw(25) << "\t-n --run_number" << "Run number.\n";
   cout << setw(25) << "\t-f --fit" << "Turn on fit.\n";
   cout << setw(25) << "\t-d --diagnostics" << "Turn on diagnostics.\n";
   cout << setw(25) << "\t-e --templates" << "Turn on template plots.\n";
   cout << setw(25) << "\t-x --learnparams" << "Do a fit to learn the GP hyperparameters.\n";
   cout << setw(25) << "\t-a --data" << "Run the fit on data (use full mc for training).\n";
   cout << setw(25) << "\t-p --profile" << "Run the likelihood profile.\n";
   cout << setw(25) << "\t-m --masspnt  <value>" << "If running on mc, use masspoint indicated.\n";
   cout << setw(25) << "\t-o --bootstrap" << "Turn on bootstrapping.\n";
   cout << setw(25) << "\t-c --fracevts" << "Fit fraction of events.\n";
   cout << setw(25) << "\t-s --numPE" << "Number of PEs for stat error validation.\n";
   cout << setw(25) << "\t-i --PE" << "PE index (out of numPE total).\n";
   cout << setw(25) << "\t-b --mbl" << "Activate Mbl distribution.\n";
   cout << setw(25) << "\t-t --mt2_220" << "Activate MT2 220 distribution.\n";
   cout << setw(25) << "\t-k --mt2_221" << "Activate MT2 221 distribution.\n";
   cout << setw(25) << "\t-2 --maos220" << "Activate MAOS 220 distribution.\n";
   cout << setw(25) << "\t-1 --maos210" << "Activate MAOS 210 distribution.\n";
   cout << setw(25) << "\t-9 --syst <string>" << "Run with a systematic variation.\n";
   cout << setw(25) << "\t-q --jfactor" << "JES-stable fit.\n";
   cout << setw(25) << "\t-g --outdir <string>" << "Output fit results.\n";
   cout << setw(25) << "\t-h --help" << "Display this menu.\n";
   cout << endl;
   return;

}

int main(int argc, char* argv[]){

   freopen ("stderr.txt","w",stderr);

   // declarations
   Fitter fitter;
   map<string, Dataset> datasets;
   vector<Event> eventvec_datamc;
   //vector<Event> eventvec_train;
   //vector<Event> eventvec_trainUP;
   //vector<Event> eventvec_trainDN;
   vector<Event> eventvec_test;
   map< string, map<string, TH1D*> > hists_all_;
   map< string, map<string, TH1D*> > hists_train_;
   //map< string, map<string, TH1D*> > hists_trainUP_;
   //map< string, map<string, TH1D*> > hists_trainDN_;
   map< string, map<string, TH1D*> > hists_test_;
   map< string, map<string, TH2D*> > hists2d_all_;
   map< string, map<string, TH2D*> > hists2d_train_;
   //map< string, map<string, TH2D*> > hists2d_trainUP_;
   //map< string, map<string, TH2D*> > hists2d_trainDN_;
   map< string, map<string, TH2D*> > hists2d_test_;

   vector<map< string, map<string, TH1D*> > > hists_jvec_train_;
   vector<map< string, map<string, TH2D*> > > hists2d_jvec_train_;

   // output fit results
   int run_number=-1;
   int fitstatus=-1;
   double mt=0, mt_err=0;
   double kmbl=0, kmbl_err=0;
   double k220=0, k220_err=0;
   double k221=0, k221_err=0;
   double kmaos220=0, kmaos220_err=0;
   double kmaos210=0, kmaos210_err=0;
   double jesfactor=0, jesfactor_err=0;
   bool distcut220 = 0;
   bool etadisamb220 = 0;
   bool blmatch220 = 0;
   bool distcut210 = 0;
   bool etadisamb210 = 0;
   bool blmatch210 = 0;
   double mcmass=0;
   double fitchi2=0;
   double tsig_mbl_chi2 [NMP] = {0};
   double tbkg_mbl_chi2 [NMP] = {0};
   string nsyst = "";
   int statval_numPE = -1;
   int statval_PE = -1;
   vector<int> evlist;
   
   TTree *tree = new TTree("FitResults", "FitResults");
   tree->Branch("runNumber", &run_number);
   tree->Branch("fitStatus", &fitstatus);
   tree->Branch("mt", &mt);
   tree->Branch("mt_err", &mt_err);
   tree->Branch("kmbl", &kmbl);
   tree->Branch("kbml_err", &kmbl_err);
   tree->Branch("k220", &k220);
   tree->Branch("k220_err", &k220_err);
   tree->Branch("k221", &k221);
   tree->Branch("k221_err", &k221_err);
   tree->Branch("kmaos220", &kmaos220);
   tree->Branch("kmaos220_err", &kmaos220_err);
   tree->Branch("kmaos210", &kmaos210);
   tree->Branch("kmaos210_err", &kmaos210_err);
   tree->Branch("jesfactor", &jesfactor);
   tree->Branch("jesfactor_err", &jesfactor_err);
   tree->Branch("distcut220", &distcut220);
   tree->Branch("etadisamb220", &etadisamb220);
   tree->Branch("blmatch220", &blmatch220);
   tree->Branch("distcut210", &distcut210);
   tree->Branch("etadisamb210", &etadisamb210);
   tree->Branch("blmatch210", &blmatch210);
   tree->Branch("mcmass", &mcmass);
   tree->Branch("fitchi2", &fitchi2);
   tree->Branch("tsig_mbl_chi2", tsig_mbl_chi2, "tsig_mbl_chi2[7]/D");
   tree->Branch("tbkg_mbl_chi2", tbkg_mbl_chi2, "tbkg_mbl_chi2[7]/D");
   tree->Branch("syst", &nsyst);
   tree->Branch("statval_numPE", &statval_numPE);
   tree->Branch("statval_PE", &statval_PE);
   tree->Branch("evlist", &evlist);

   // option flags
   int c;
   int do_fit = 0;
   int do_diagnostics = 0;
   int use_data = 0;
   float masspnt = 0;
   int do_bootstrap = 0;
   int do_templates = 0;
   int do_learnparams = 0;
   int do_mbl = 0;
   int do_mt2_220 = 0;
   int do_mt2_221 = 0;
   int do_maos220 = 0;
   int do_maos210 = 0;
   int maoscuts220 = 0;
   int maoscuts210 = 0;
   double fracevts = -1;
   string run_number_str = "";
   string outdir = "results";

   struct option longopts[] = {
      { "run_number",   required_argument,   0,                'n' },
      { "fit",          no_argument,         &do_fit,          'f' },
      { "diagnostics",  no_argument,         &do_diagnostics,  'd' },
      { "templates",    no_argument,         &do_templates,    'e' },
      { "learnparams",  no_argument,         &do_learnparams,  'x' },
      { "data",         no_argument,         &use_data,        'a' },
      { "profile",      no_argument,         0,                'p' },
      { "masspnt",      required_argument,   0,                'm' },
      { "bootstrap",    no_argument,         &do_bootstrap,    'o' },
      { "fracevts",     required_argument,   0,                'c' },
      { "numPE",        required_argument,   0,                's' },
      { "PE",           required_argument,   0,                'i' },
      { "mbl",          no_argument,         &do_mbl,          'b' },
      { "mt2_220",      no_argument,         &do_mt2_220,      'k' },
      { "mt2_221",      no_argument,         &do_mt2_221,      't' },
      { "maos220",      no_argument,         &do_maos220,      '2' },
      { "maos210",      no_argument,         &do_maos210,      '1' },
      { "maoscuts220",  required_argument,   0,                'y' },
      { "maoscuts210",  required_argument,   0,                'z' },
      { "jfactor",      no_argument,         0,                'q' },
      // maoscuts220 and maoscuts210 set which cuts to use for maos220 and maos 210 respectively.
      // see lines 237-251 for what number to input.
      { "syst",         required_argument,   0,                '9' },
      { "outdir",      required_argument,   0,                'g' },
      { "help",         no_argument,         NULL,             'h' },
      { 0, 0, 0, 0 }
   };

   while( (c = getopt_long(argc, argv, "fdexahponbkt21qyzm:c:s:i:9:g:", longopts, NULL)) != -1 ) {
      switch(c)
      {
         case 'n' :
            run_number = atoi(optarg);
            run_number_str = "_"+string(optarg);
            break;

         case 'f' :
            do_fit = 1;
            break;

         case 'd' :
            do_diagnostics = 1;
            break;

         case 'e' :
            do_templates = 1;
            break;

         case 'x' :
            do_learnparams = 1;
            break;

         case 'a' :
            use_data = 1;
            break;

         case 'm' :
            masspnt = atof(optarg);
            break;

         case 'p' :
            fitter.compute_profile = true;
            break;

         case 'o' :
            do_bootstrap = 1;
            break;

         case 'c' :
            fracevts = atof(optarg);
            break;

         case 's' :
            statval_numPE = atoi(optarg);
            break;

         case 'i' :
            statval_PE = atoi(optarg);
            break;

         case 'b' :
            do_mbl = true;
            break;

         case 't' :
            do_mt2_220 = true;
            break;

         case 'k' :
            do_mt2_221 = true;
            break;

         case '2' :
            do_maos220 = true;
            break;

         case '1' :
            do_maos210 = true;
            break;

         case 'y' :
            maoscuts220 = atoi(optarg);
            break;

         case 'z' :
            maoscuts210 = atoi(optarg);
            break;

         case '9' :
            nsyst = optarg;
            break;

         case 'g' :
            outdir = optarg;
            break;

         case 'q' :
            fitter.fit_jfactor = true;
            break;

         case 'h' :
            print_usage();
            return -1;
            break;

         case 0:     /* getopt_long() set a variable, just keep going */
            break;

         case ':':   /* missing option argument */
            fprintf(stderr, "%s: option `-%c' requires an argument\n",
                  argv[0], optopt);
            return -1;

         case '?':
         default:    /* invalid option */
            fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
                  argv[0], optopt);
            print_usage();
            return -1;

      }
   }

   fitter.InitializeDists();

   fitter.dists["mbl_gp"].activate = do_mbl;
   fitter.dists["mt2_220_gp"].activate = do_mt2_220;
   fitter.dists["mt2_221_gp"].activate = do_mt2_221;
   fitter.dists["maos220_gp"].activate = do_maos220;
   fitter.dists["maos210_gp"].activate = do_maos210;
   if( do_maos220 ) fitter.compute_maos220 = true;
   if( do_maos210 ) fitter.compute_maos210 = true;

   if (maoscuts220 == 1){ distcut220 = 1; }
   else if (maoscuts220 == 2){ etadisamb220 = 1; }
   else if (maoscuts220 == 3){distcut220 = 1; etadisamb220 = 1; }
   else if (maoscuts220 == 4){blmatch220 = 1; }
   else if (maoscuts220 == 5){distcut220 = 1; blmatch220 = 1; }
   else if (maoscuts220 == 6){etadisamb220 = 1; blmatch220 = 1; }
   else if (maoscuts220 == 7){distcut220 = 1; etadisamb220 = 1; blmatch220 = 1; }

   if (maoscuts210 == 1){ distcut210 = 1; } 
   else if (maoscuts210 == 2){ etadisamb210 = 1; }
   else if (maoscuts210 == 3){distcut210 = 1; etadisamb210 = 1; }
   else if (maoscuts210 == 4){blmatch210 = 1; }
   else if (maoscuts210 == 5){distcut210 = 1; blmatch210 = 1; } 
   else if (maoscuts210 == 6){etadisamb210 = 1; blmatch210 = 1; }
   else if (maoscuts210 == 7){distcut210 = 1; etadisamb210 = 1; blmatch210 = 1; }
   
   fitter.maoscuts220 = maoscuts220;
   fitter.maoscuts210 = maoscuts210;

   // Check that at least one kinematic variable's lengthscale has been entered.
   // Any additional distributions need to be added here
   if (!(fitter.dists["mbl_gp"].activate) and !(fitter.dists["mt2_220_gp"].activate) and !(fitter.dists["maos220_gp"].activate) and !(fitter.dists["maos210_gp"].activate) and (do_fit == 1 or do_templates == 1) ){
      std::cout << "At least one variable needed to do fit.  Input at least one lengthscale." << std::endl;
      print_usage();
      return -1;
   }

   fitter.LoadDatasets( datasets );

   // for event counting
   map<string, int> datacount;

   // random number seed for bootstrapping (turns on when nonzero)
   int randseed = 0;
   if( do_bootstrap ) randseed = run_number+1+10E6;

   // ********************************************************
   // events for diagnostics
   // ********************************************************
   if( do_diagnostics ){
      fitter.ReadDatasets( datasets, eventvec_datamc, "diagnostics", nsyst, fracevts, statval_numPE, statval_PE );
      fitter.GetVariables( eventvec_datamc );
      fitter.DeclareHists( hists_all_, hists2d_all_, "all" );
      fitter.FillHists( hists_all_, hists2d_all_, eventvec_datamc );
      fitter.PrintHists( hists_all_, hists2d_all_ );
      return 0;
   }

   // ********************************************************
   // events for GP training
   // ********************************************************
   
   //cout << "\nLoading datasets: training" << endl;
   //fitter.ReadDatasets( datasets, eventvec_train, "train", nsyst, fracevts, statval_numPE, statval_PE );

   // pdf systematics
   /*
   int pdfvar = -1;
   if( nsyst.find("PDFvar") != string::npos ){
         string nametemp = nsyst;
         nametemp.erase(0,6);
         pdfvar = atoi( nametemp.c_str() );
         fitter.PDFReweight( eventvec_train, pdfvar );
   }

   fitter.GetVariables( eventvec_train );
   fitter.DeclareHists( hists_train_, hists2d_train_, "train" );
   fitter.FillHists( hists_train_, hists2d_train_, eventvec_train );
   fitter.FindPTrain( hists_train_ );
   */

   hists_jvec_train_.resize(5);
   hists2d_jvec_train_.resize(5);

   string systs [] = {"JESTotalDN", "Central", "JESTotalUP"};
   for(int i=0; i < 5; i++){
      cout << "\nLoading datasets: training " << i << endl;
      vector<Event> eventvec_temp;

      //for(vector<Event>::iterator ev = eventvec_temp.begin(); ev != eventvec_temp.end(); ev++){
      //   eventvec_temp.push_back( *ev );
      //}
      //fitter.ReadDatasets( datasets, eventvec_temp, "train", systs[i], fracevts, statval_numPE, statval_PE );
      fitter.ReadDatasets( datasets, eventvec_temp, "train", nsyst, fracevts, statval_numPE, statval_PE );

      int pdfvar = -1;
      if( nsyst.find("PDFvar") != string::npos ){
         string nametemp = nsyst;
         nametemp.erase(0,6);
         pdfvar = atoi( nametemp.c_str() );
         fitter.PDFReweight( eventvec_temp, pdfvar );
      }

      ostringstream jstring;
      jstring << 100*fitter.jfactpoints[i];

      fitter.JShift( eventvec_temp, fitter.jfactpoints[i] );
      fitter.GetVariables( eventvec_temp );
      fitter.DeclareHists( hists_jvec_train_[i], hists2d_jvec_train_[i], "train"+jstring.str() );
      fitter.FillHists( hists_jvec_train_[i], hists2d_jvec_train_[i], eventvec_temp );

      // release memory in eventvec_temp
      vector<Event>().swap( eventvec_temp );

   }

   fitter.FindPTrain( hists_jvec_train_[2] );


   // release memory in eventvec_train
   //vector<Event>().swap( eventvec_train );


   // ********************************************************
   // events for fit
   // ********************************************************
   
   if( do_fit ){
      fitter.ReadDatasets( datasets, eventvec_test, "test", nsyst, fracevts, statval_numPE, statval_PE );
      fitter.GetVariables( eventvec_test );
      fitter.DeclareHists( hists_test_, hists2d_test_, "test" );
      fitter.FillHists( hists_test_, hists2d_test_, eventvec_test );
   }


   // ********************************************************
   // begin fit
   // ********************************************************
   if( do_fit ){

      vector<Event> eventvec_fit;
      map< string, map<string, TH1D*> > hists_fit_;
      map< string, map<string, TH2D*> > hists2d_fit_;

      if( use_data ){ // added 220 distribution to the data fit too, but haven't tested it

         //
         // turn this feature off for now -- will need to clean up later.
         //

         /*
         cout << "REMINDER: EVENT WEIGHTS IN MC" << endl;
         for(vector<Event>::iterator ev = eventvec_datamc.begin(); ev < eventvec_datamc.end();ev++){
            if( ev->type.find("data") != string::npos ){
               if( ev->type.find("bkgcontrol") == string::npos ){
                  eventvec_fit.push_back(*ev);
               }else{
                  eventvec_fit_bkgcontrol.push_back(*ev);
               }
            }
         }

         // flag events to be fitted
         for( vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
            ev->fit_event = true;
            for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
               if( ev->process.compare(it->first) == 0 ) datacount[it->first]+=1;
            }
         }
         for( vector<Event>::iterator ev = eventvec_fit_bkgcontrol.begin();
               ev < eventvec_fit_bkgcontrol.end(); ev++){
            ev->fit_event = true;
            for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
               if( ev->process.compare(it->first) == 0 ) datacount[it->first]+=1;
            }
         }
         cout << "Fit event count: " << endl;
         cout << setiosflags(ios::left);
         for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
            cout << "... " << setw(25) << it->first << ": " << it->second << " events" << endl;
         }

         fitter.DeclareHists( hists_fit_, "fit" );
         fitter.FillHists( hists_fit_, eventvec_fit, true );
         fitter.DeclareHists( hists_fit_bkgcontrol_, "fit_bkgcontrol" );
         fitter.FillHists( hists_fit_bkgcontrol_, eventvec_fit_bkgcontrol, true );

         double wgt=0;
         for(vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
            wgt += ev->weight;
         }
         cout << "wgt = " << wgt << endl;

         // do GP training
         cout << "Training GP." << endl;

         // GP training for mbl
         if ( lengthscale_mbl != -1){
            Shapes * fptr = new Shapes( "mbl", fitter.gplength_mbl, fitter.gplength_mt_mbl,
                  fitter.lbnd, fitter.rbnd );
            fptr->TrainGP( hists_train_ );
            fitter.aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
            fitter.aGPsig = fptr->aGPsig;
            fitter.aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
            fitter.aGPbkg = fptr->aGPbkg;
         }

         // GP training for 220
         if ( lengthscale_220 != -1){
            Shapes * fptr220 = new Shapes( "mt2_220_nomatchmbl", fitter.gplength_220, fitter.gplength_mt_220,
                  fitter.lbnd, fitter.rbnd );
            fptr220->TrainGP( hists_train_ );
            fitter.aGPsig220.ResizeTo( fptr220->aGPsig.GetNoElements() );
            fitter.aGPsig220 = fptr220->aGPsig;
            fitter.aGPbkg220.ResizeTo( fptr220->aGPbkg.GetNoElements() );
            fitter.aGPbkg220 = fptr220->aGPbkg;
         }

         // events for fitting, hists for training
         fitter.RunMinimizer( eventvec_fit );
         fitter.PlotResults( hists_fit_ ); // plot fitted events

         // fill results tree
         mcmass = 0;
         fitstatus = fitter.gMinuit->Status();
         const double *par = fitter.gMinuit->X();
         const double *par_err = fitter.gMinuit->Errors();
         mt = par[0];
         kmbl = par[1];
         k220 = par[2];
         mt_err = par_err[0];
         kmbl_err = par_err[1];
         k220_err = par_err[2];
         fitchi2 = fitter.fitchi2;
         gplength_mbl = lengthscale_mbl;
         gplength_220 = lengthscale_220;
         gplength_mt_mbl = lengthscale_mt_mbl;
         gplength_mt_220 = lengthscale_mt_220;

         // TODO
         //tree->Fill();

         eventvec_fit.clear();
         eventvec_fit_bkgcontrol.clear();
         fitter.DeleteHists( hists_fit_ );
         fitter.DeleteHists( hists_fit_bkgcontrol_ );
         */

      }else{ // loop over mc masses

         double masspnts [] = {166.5,169.5,171.5,172.5,173.5,175.5,178.5};
         for(int i=0; i < 7; i++){

            double mass = masspnts[i];
            // masspoint from command line
            if( masspnt != 0 ){
               bool check = false;
               for(int j=0; j < 7; j++){
                  if( masspnts[j] == masspnt ) check = true;
               }
               if(check){
                  mass = masspnt;
                  i = 6;
               }else{
                  cout << "masspoint " << masspnt << " not found!" << endl;
                  return -1;
               }
            }

            cout << "############# Fitting Masspoint " << mass << " ###############" << endl;

            stringstream dstr;
            dstr << floor(mass);
            string dname = "ttbar"+dstr.str();

            // load events to be fitted
            for(vector<Event>::iterator ev = eventvec_test.begin(); ev < eventvec_test.end(); ev++){
               if( ev->type.find(dname) != string::npos or ev->type.find("other") != string::npos ){
                  eventvec_fit.push_back(*ev);
               }
            }

            // release memory in eventvec_test
            vector<Event>().swap( eventvec_test );

            fitter.ReweightMC( eventvec_fit, dname );
            bool statval = statval_numPE != -1;
            if( do_bootstrap ){
               evlist = fitter.Resample( eventvec_fit, randseed, statval );
            }else{
               evlist = fitter.Resample( eventvec_fit, 1+10E6, 0 );
               cout << "RESAMPLING TO PE #1" << endl;
            }

            // flag events to be fitted
            for( vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
               ev->fit_event = true;
               for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
                  if( ev->process.compare(it->first) == 0 ) datacount[it->first]+=1;
               }
            }

            cout << "Fit event count: " << endl;
            cout << setiosflags(ios::left);
            for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
               cout << "... " << setw(25) << it->first << ": " << it->second << " events" << endl;
            }

            fitter.DeclareHists( hists_fit_, hists2d_fit_, "fit" );
            fitter.FillHists( hists_fit_, hists2d_fit_, eventvec_fit, true );

            // do GP training
            for( map<string, Distribution>::iterator it = fitter.dists.begin(); it != fitter.dists.end(); it++ ){

               string name = it->first;
               Distribution *dist = &(it->second);

               double m2llsig, m2llbkg;

               if( dist->activate ){
                  Shapes * fptr = new Shapes( name, dist->ptrain,
                        dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->lbnd, dist->rbnd );
                  cout << "aGP CENTRAL" << endl;
                  fptr->TrainGP( hists_jvec_train_, m2llsig, m2llbkg );
                  cout << "done" << endl;

                  dist->aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
                  dist->aGPsig = fptr->aGPsig;
                  dist->aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
                  dist->aGPbkg = fptr->aGPbkg;

                  delete fptr;

                  /*
                  if( fitter.fit_jfactor ){
                     Shapes * fptr2 = new Shapes( name, dist->ptrain,
                           dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->lbnd, dist->rbnd );
                     cout << "aGP UP" << endl;
                     fptr2->TrainGP( hists_trainUP_, m2llsig, m2llbkg );
                     cout << "done" << endl;

                     dist->aGPsigUP.ResizeTo( fptr2->aGPsig.GetNoElements() );
                     dist->aGPsigUP = fptr2->aGPsig;
                     dist->aGPbkgUP.ResizeTo( fptr2->aGPbkg.GetNoElements() );
                     dist->aGPbkgUP = fptr2->aGPbkg;

                     delete fptr2;

                     Shapes * fptr3 = new Shapes( name, dist->ptrain,
                           dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->lbnd, dist->rbnd );
                     cout << "aGP DN" << endl;
                     fptr3->TrainGP( hists_trainDN_, m2llsig, m2llbkg );
                     cout << "done" << endl;

                     dist->aGPsigDN.ResizeTo( fptr3->aGPsig.GetNoElements() );
                     dist->aGPsigDN = fptr3->aGPsig;
                     dist->aGPbkgDN.ResizeTo( fptr3->aGPbkg.GetNoElements() );
                     dist->aGPbkgDN = fptr3->aGPbkg;

                     delete fptr3;
                  }
                  */
               }

            }

            typedef map<string, TH1D*> tmap;
            typedef map<string, tmap> hmap;
            for( hmap::iterator h = hists_train_.begin(); h != hists_train_.end(); h++){
               for( tmap::iterator t = h->second.begin(); t != h->second.end(); t++){
                  for(int n=1; n < t->second->GetNbinsX(); n++){
                     if( t->second->GetBinContent(n) == 0 )
                        t->second->SetBinError(n, 1.0/35000);
                  }
               }
            }

            // events for fitting, hists for training
            fitter.RunMinimizer( eventvec_fit );
            fitter.PlotResults( hists_fit_, outdir+"/plotsFitResults"+run_number_str+".root" ); // plot fitted events

            cout << "Fit Chi2 = " << fitter.fitchi2 << endl;


            // fill results tree
            // any additional variables need to be added here
            mcmass = mass;
            fitstatus = fitter.gMinuit->Status();
            const double *par = fitter.gMinuit->X();
            const double *par_err = fitter.gMinuit->Errors();
            mt = par[0];
            mt_err = par_err[0];
            jesfactor = par[1];
            jesfactor_err = par_err[1];
            kmbl = par[2];
            kmbl_err = par_err[2];
            k220 = par[3];
            k220_err = par_err[3];
            kmaos220 = par[4];
            kmaos220_err = par_err[4];
            kmaos210 = par[5];
            kmaos210_err = par_err[5]; 
            k221 = par[6];
            k221_err = par_err[6];
            fitchi2 = fitter.fitchi2;

            eventvec_fit.clear();
            fitter.DeleteHists( hists_fit_, hists2d_fit_ );
         }

      }

   }

   if( do_templates ){

            // do GP training
            for( map<string, Distribution>::iterator it = fitter.dists.begin(); it != fitter.dists.end(); it++ ){

               string name = it->first;
               Distribution *dist = &(it->second);

               double m2llsig, m2llbkg;

               if( dist->activate ){
                  Shapes * fptr = new Shapes( name, dist->ptrain,
                        dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->lbnd, dist->rbnd );
                  fptr->TrainGP( hists_jvec_train_, m2llsig, m2llbkg );

                  dist->aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
                  dist->aGPsig = fptr->aGPsig;
                  dist->aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
                  dist->aGPbkg = fptr->aGPbkg;

                  dist->Ainv_sig.ResizeTo( fptr->aGPsig.GetNoElements(), fptr->aGPsig.GetNoElements() );
                  dist->Ainv_sig = fptr->Ainv_sig;
                  dist->Ainv_bkg.ResizeTo( fptr->aGPbkg.GetNoElements(), fptr->aGPbkg.GetNoElements() );
                  dist->Ainv_bkg = fptr->Ainv_bkg;

                  delete fptr;
               }

            }

            fitter.PlotTemplates( hists_jvec_train_ );

            for(int j=0; j < NMP; j++){
               tsig_mbl_chi2[j] = fitter.tsig_mbl_chi2[j];
               tbkg_mbl_chi2[j] = fitter.tbkg_mbl_chi2[j];
            }

   }

   if( do_learnparams ){
      for( map<string, Distribution>::iterator it = fitter.dists.begin(); it != fitter.dists.end(); it++ ){
         
         string name = it->first;
         Distribution *dist = &(it->second);

         if( dist->activate ){
            cout << "Learning hyperparameters for distribution " << name << "..." << endl;
            Shapes * fptrtmp = new Shapes( name, dist->ptrain,
                  dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->lbnd, dist->rbnd );
            fptrtmp->LearnGPparams( hists_jvec_train_ );

            dist->glx = fptrtmp->lx;
            dist->glmt = fptrtmp->lmass;
            dist->gnorm1 = fptrtmp->gnorm1;
            dist->gnorm2 = fptrtmp->gnorm2;

            double m2llsig, m2llbkg;
            Shapes * fptr = new Shapes( name, dist->ptrain,
                  dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2, dist->lbnd, dist->rbnd );
            fptr->TrainGP( hists_jvec_train_, m2llsig, m2llbkg );

            dist->aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
            dist->aGPsig = fptr->aGPsig;
            dist->aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
            dist->aGPbkg = fptr->aGPbkg;

            dist->Ainv_sig.ResizeTo( fptr->aGPsig.GetNoElements(), fptr->aGPsig.GetNoElements() );
            dist->Ainv_sig = fptr->Ainv_sig;
            dist->Ainv_bkg.ResizeTo( fptr->aGPbkg.GetNoElements(), fptr->aGPbkg.GetNoElements() );
            dist->Ainv_bkg = fptr->Ainv_bkg;

            cout << "begin PlotTemplates" << endl;
            fitter.PlotTemplates( hists_jvec_train_ );
            cout << "end PlotTemplates" << endl;
            delete fptr;
            delete fptrtmp;
         }
      }
   }

   if( do_fit or do_templates ){
      tree->Fill();
   }


   //
   // write fit results
   //
   if( do_fit or do_templates ){

      TFile *file = new TFile((outdir+"/fitresults"+run_number_str+".root").c_str(), "RECREATE");
      file->cd();
      tree->Write();
      file->Write();
      file->Close();
   }

   return 0;
}

