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
  //here must be called cuts that cannot be called in InitCuts in roder to define variables
  Setup_3piVertexQA();
  Setup_KaonMatching();

 
  //plotting modifiers
  int rebin=1;
  bool ignoreRange=false; // change to spot some outlayers
  const bool normalize=false; //plot normalized
  gIgnoreCutMods=false;

   // structure for results
  ResultList1D Res_EventPlots;
  ResultList2D Res_EventPlots_2D; 
  ResultList1D Res_3piPlots;
  ResultList2D Res_3piPlots_2D; 
  ResultList1D Res_KaonPlots; 
  ResultList2D Res_KaonPlots_2D; 
   
  //loop over datasets
  for (int iFile=0;iFile<nFiles;iFile++){
   
   TString fname=files[order[iFile]].fileName;
    cout<<"opening file: "<<fname<<endl;

    auto chain_kaons = new TChain("kaons");
    auto  fcount=chain_kaons->Add(fname,nEntriesLimit);
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

    auto after_evCut_node = kaons_node.Filter(evCut.Str());
    cout<<" Event cut used:  "<<endl <<evCut.Str()<<endl<<endl;
  
     
    //matched kaons with cut on mother PID from KFP (100321 for K+) both real and MC
    auto MatchedKaon_cut=K3piCut_Matched_Kplus();
    //for MC filed I only add cut on MC vertex !!!!not on the KAON - allows to study mis matches
    if (files[order[iFile]].isMc) {MatchedKaon_cut=K3piCut_Matched_Kplus() + Setup_MCvertex();}
    cout<<endl<<"matched K+ 3piVtx_ cut used:  "<<endl<< MatchedKaon_cut.Str()<<endl<<endl;
 
     
    // event plot per matched K
    auto Cut=MatchedKaon_cut+evCut;
    Res_EventPlots.resetPosition();
    Res_EventPlots_2D.resetPosition(); 
    //this takes long
    AddPlots4QA(Event_plots,kaons_node,Cut,Res_EventPlots,"per matched K+",files[order[iFile]].lable,rebin,false);
    AddPlots4QA(Event_plots_2D,kaons_node,Cut,Res_EventPlots_2D,"per matched K+",files[order[iFile]].lable,rebin,false);
   
    //properties of 3pi vertex with matched K
    Res_3piPlots.resetPosition();
    Res_3piPlots_2D.resetPosition(); 
    //3pi vertex 
    AddPlots4QA(RecoVtx_plots,after_evCut_node,MatchedKaon_cut,Res_3piPlots,"per matched K+",files[order[iFile]].lable,rebin,false);
    AddPlots4QA(RecoVtx_plots_2D,after_evCut_node,MatchedKaon_cut,Res_3piPlots_2D,"per atched K+",files[order[iFile]].lable,rebin,false);


    //Properties of matched kaons
    Res_KaonPlots.resetPosition(); 
    Res_KaonPlots_2D.resetPosition(); 
    AddPlots4QA(Kaon_plots,after_evCut_node,MatchedKaon_cut,Res_KaonPlots,"of matched K+",files[order[iFile]].lable,rebin,false);
    AddPlots4QA(Kaon_plots_2D,after_evCut_node,MatchedKaon_cut,Res_KaonPlots_2D,"of matched K+",files[order[iFile]].lable,rebin,false);
   


   cout<<"trigger lazy evaluation"<<endl;
 
   auto ct= kaons_node.Count();
   ct.OnPartialResult(/*every */100000/* events*/,
                           [](auto c) { std::cout << c << '\n'; });
   cout<<*ct<<endl;

   cout<<"trigger DONE"<<endl;
 
  
   delete chain_kaons;
 
} //loop over files

  TFile *f=new TFile("matchedComp_2019_SL24_vs_SL23_noEmbed_noCuts.root","recreate");
  f->mkdir("events");f->cd("events"); 
  DrawResults(Res_EventPlots);
  DrawResults(Res_EventPlots_2D); 
  f->mkdir("3pi");f->cd("3pi"); 
  DrawResults(Res_3piPlots);
  DrawResults(Res_3piPlots_2D); 
  f->mkdir("kaons");f->cd("kaons"); 
  DrawResults(Res_KaonPlots); 
  DrawResults(Res_KaonPlots_2D);
f->Write();
//f.Close(); //dono tclose to see resutls


return;
}



