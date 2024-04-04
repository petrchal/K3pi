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
#include "K3piPlotDefs.C"
#include "data_Include.h"
#include "cuts_Include.h"

using namespace ROOT;

 //int particlesPDG[nPIDS] = {100321, 200321,-100321, -200321}; //k->3pi only, K->3pi with found K


void QA_matchedKaon(){
  gROOT->Reset();
  ROOT::EnableImplicitMT(); //enambe multi threading - application must be MT safe ..your job
  InitCuts();

  //2017 54 data 

   
  //K3piCut_EventCut =EventCut_2021_7p7AuAu; //just in case

  //plotting modifiers
  int rebin=1;
  bool ignoreRange=false; // change to spot some outlayers
  const bool normalize=false; //plot normalized
 

  // structure for results
  ResultList1D Res_Plots; 
   
  //loop over datasets
  for (int iFile=0;iFile<nFiles;iFile++){
   
   TString fname=files[order[iFile]].fileName;
    cout<<"opening file: "<<fname<<endl;

    auto chain_kaons = new TChain("kaons");
    auto  fcount=chain_kaons->Add(fname,nEntriefsLimit);
    cout<<" TChain Added "<<fcount<<" files from "<<fname<<endl;
    //TObjArray * ll=chain_kaons->GetListOfFiles();

    ROOT::RDF::RNode  kaons_node= ROOT::RDF::AsRNode(RDataFrame(*chain_kaons)); //raw event count
  
    //common event cuts
    auto  evCut=K3piCut_EventCut();
    if (files[order[iFile]].trigger){ //set propper trigers 
       evCut.Replace(files[order[iFile]].trigger);
    }

   //must be called after all cuts are read at least once
     kaons_node=DefineNewVariables(kaons_node);
     //kaons_node=AddVariations(vary_EvtVz,kaons_node);
     //kaons_node=AddVariations(vary_EvtVz,kaons_node);
     //kaons_node=AddVariations(vary_lastPointDiff,kaons_node);
   
    // if (files[order[iFile]].isMc) kaons_node=AddVariations(vary_DCAxy,kaons_node);

    auto d_events = kaons_node.Filter(evCut.Str());
    cout<<" Event cut used:  "<<endl <<evCut.Str()<<endl<<endl;
  
    //all 3pi candidates ... all together MC and non MC, positive and negative !!!
    //auto d_3piCandidates = d_events.Filter(K3piCut_VertexCanditate().Str());
     
    auto CurrentPos=Res_Plots.begin();
    
    //matched kaons with cut on mother PID from KFP (100321 for K+) both real and MC
    auto kaon_cut=K3piCut_Matched_Kplus();
    //for MC filed I only add cut on MC vertex !!!!not on the KAON - allows to study mis matches
    if (files[order[iFile]].isMc) {kaon_cut=K3piCut_Matched_Kplus() + Setup_MCvertex();}
    cout<<endl<<" K+ Mother cut used:  "<<endl<< kaon_cut.Str()<<endl<<endl;
 
     
    // event plot per matched K
    auto ev_cut=kaon_cut+evCut;
    //AddPlots4QA(Event_plots,kaons_node,ev_cut,Res_Plots,CurrentPos,"per matched K+",files[order[iFile]].lable,rebin,false);
   
    //Mother(3pi vertex) per matched K
   // AddPlots4QA(Mother_plots,d_events,kaon_cut,Res_Plots,CurrentPos,"per matched K+",files[order[iFile]].lable,rebin,false);
      
    //Properties of matched kaon
    AddPlots4QA(MatchedKaon_plots,d_events,kaon_cut,Res_Plots,CurrentPos,"of matched K+",files[order[iFile]].lable,rebin,false);
    

 // .. it is a problem, since histogram may come from different trees (nodes)
  //AddProgressBar(event_node);

   cout<<"trigger lazy evaluation"<<endl;
 
   auto ct= kaons_node.Count();
   
   cout<<*ct<<endl;

   cout<<"trigger DONE"<<endl;
 
  
   delete chain_kaons;
 
} //loop over files

DrawResults(Res_Plots);

return;
}



