#include "TNtuple.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEfficiency.h"
#include <iostream>
#include "K3piLib.C"
#include "FXT_K3piPlotDefs.C"
#include "FXT_data_Include.h"
#include "FXT_cuts_Include.h"

using namespace ROOT;

 //int particlesPDG[nPIDS] = {100321, 200321,-100321, -200321}; //k->3pi only, K->3pi with found K


void FXT_QA_3piVtx(){
  gROOT->Reset();
  ROOT::EnableImplicitMT(); //enambe multi threading - application must be MT safe ..your job

  //one must setup (Call) cuts here - the varibles get defined here
  K3piCut_EventCut =EventCut_2020_FXT; //just in case override 
  InitCuts();
  //these cannot be inititalized in InitCuts
  K3piCut_3piVtx_Kplus(); 
  Setup_MCvertex();


  //plotting modifiers
  int rebin=1;
  //ignires ranges set in plot definitions
  bool ignoreRange=false; // change to spot some outlayers
  //const bool normalize=false; //plot normalized
 

  // structure for results
  ResultList1D Res_Plots; 
  ResultList1D Res_eventPlots; 
   
  //loop over datasets
  for (int iFile=0;iFile<nFiles;iFile++){
   
   TString fname=files[order[iFile]].fileName;
    cout<<"opening file: "<<fname<<endl;

    auto chain_kaons = new TChain("kaons");
    auto  fcount=chain_kaons->Add(fname,nEntriefsLimit);
    cout<<" TChain Added "<<fcount<<" files from "<<fname<<endl;
    //TObjArray * ll=chain_kaons->GetListOfFiles();

    //root dataframe
    ROOT::RDF::RNode  kaons_node= ROOT::RDF::AsRNode(RDataFrame(*chain_kaons)); //raw event count
  
    //common event cuts
    auto  evCut=K3piCut_EventCut();
    if (files[order[iFile]].trigger){ //set propper trigers and apply data specific cuts
       evCut.Replace(files[order[iFile]].trigger);
    }

    //must be called after all cuts are read
    //initialized new variables defined uin FXt_cuts_include.h
    kaons_node=DefineNewVariables(kaons_node);

    //enable to include cut variations 
    //kaons_node=AddVariations(vary_EvtVz,kaons_node);

    ROOT::RDF::RNode after_evCut_node = kaons_node.Filter(evCut.Str());
    cout<<" Event cut used:  "<<endl <<evCut.Str()<<endl<<endl;
    //all 3pi candidates ... all together MC and non MC, positive and negative !!!
     
      
    //3pi with cut on  PID from KFP (100321 for K+) both real and MC
    auto Reco3piVtx_cut=K3piCut_3piVtx_Kplus();
    //for MC data add condition on MC vertex.. 
    if (files[order[iFile]].isMc) {Reco3piVtx_cut=K3piCut_3piVtx_Kplus() + Setup_MCvertex();}
    cout<<endl<<" K+ 3piVtx_ cut used:  "<<endl<< Reco3piVtx_cut.Str()<<endl<<endl;
 
   //this is cumbersome
    auto CurrentPos=Res_eventPlots.begin();

    // event plot per found 3pi+
    auto Cut=Reco3piVtx_cut+evCut;
    AddPlots4QA(Event_plots_FXT,kaons_node,Cut,Res_eventPlots,CurrentPos,"events",files[order[iFile]].lable,rebin,false);
    cout<<endl<<"Cut used:  "<<endl<<  Cut.Str()<<endl<<endl;
 

   CurrentPos=Res_Plots.begin();
    //3pi vertex 
   AddPlots4QA(RecoVtx_plots,after_evCut_node,Reco3piVtx_cut,Res_Plots,CurrentPos,"per found 3pi+",files[order[iFile]].lable,rebin,false);
      
 

   cout<<"trigger lazy evaluation"<<endl;
   auto ct= kaons_node.Count();
   ct.OnPartialResult(/*every */100000/* events*/,
                           [](auto c) { std::cout << c << '\n'; });
   cout<<*ct<<endl;
   cout<<"trigger DONE"<<endl;
 
  
   delete chain_kaons;
 
} //loop over files

//save results
TFile f("3pi_FXT_cleaned.root","recreate");

f.mkdir("event info per found 3piVtx");f.cd("event info per found 3piVtx");
DrawResults(Res_eventPlots);
f.mkdir("3piVtxs");f.cd("3piVtxs");
DrawResults(Res_Plots);

f.Write();
f.Close();

return;
}



