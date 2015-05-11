#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

void test_evlist(){

   TFile *file = new TFile("test.root");
   file->cd();

   TTree *tree = (TTree*)file->Get("FitResults");
   
   vector<int> *evlist = 0;
   double mcmass = 0;
   int fitStatus = 0;

   TBranch *bvpx = 0;
   tree->SetBranchAddress("evlist",&evlist,&bvpx);
   tree->SetBranchAddress("mcmass", &mcmass);
   tree->SetBranchAddress("fitStatus", &fitStatus);

   vector<int> evstatus0;
   vector<int> evstatus3;
   for(int i=0; i < tree->GetEntries(); i++){

      tree->GetEntry(i);

      if( mcmass == 166.5 ){
         if( fitStatus == 0 ){
            for(int j=0; j < evlist->size(); j++) evstatus0.push_back( (*evlist)[j] );
         }
         if( fitStatus == 3 ){
            for(int j=0; j < evlist->size(); j++) evstatus3.push_back( (*evlist)[j] );
         }
      }

   }


   vector<int> diff(1000);
   vector<int>::iterator it;
   std::sort(evstatus0.begin(),evstatus0.end());
   std::sort(evstatus3.begin(),evstatus3.end());

   it = set_difference (evstatus3.begin(), evstatus3.end(), evstatus0.begin(), evstatus0.end(), diff.begin());
   diff.resize(it-diff.begin());

   std::cout << "Events: ";
   for(int i=0; i < diff.size(); i++) std::cout << diff[i] << " ";
   std::cout << std::endl;

   return;
}

