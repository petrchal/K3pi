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


using namespace ROOT;

 //int particlesPDG[nPIDS] = {100321, 200321,-100321, -200321}; //k->3pi only, K->3pi with found K


void QA_3pi_vs_matched(){
  gROOT->Reset();
  ROOT::EnableImplicitMT(); //enambe multi threading - application must be MT safe ..your job
  InitCuts();

  //2017 54 data 
  TFileDescription data54GeV[]={
    {"../../ntup/2017_54GeV_full/*11.root","data",0},
   };


 //27 data vs embeddings
  TFileDescription data27GeV[]={
    {"../../ntup/2018_27GeV_Sept22/*.root","data",0},
    {"../../ntup/2018_27GeV_2023Jan/*.root","data rerun",0},
    //{"../../ntup/2018_27GeV_2023Jan/*11.root","data rerun",0,"n3pi:(Evt.nK3piP==1)"},
    {"../../ntup/2018_27Embed_tuned/kaon*.root","tuned embed",1,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"},
   
    {"../../ntup/2018_27GeV_embed_ptW_1M/kaon*.root","embed",0},
    {"../../ntup/2018_27GeV_embed_ptW_1M/kaon*.root","embed isMC",1},
    {"../../ntup/2018_embed_flatPt_forcedDec/kaons.root","MC flat Pt",1},
    {"/media/petrchal/XSD-Linux/tpcAna/27GeV_Apr22/*1.root","old data",0},
    }; 


    //compare runs
    TFileDescription compare[]={
     {"../../ntup/2018_27GeV_Sept22/*.root","27GeV data",0,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"},
    //{"/media/petrchal/XSD-Linux/tpcAna/27GeV_Apr22/*1.root","old data",0},
    {"../../ntup/2018_27Embed_tuned/kaon*.root","tuned embed",1,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"},
     // {"../../ntup/2017_54GeV_full/*1.root","54Gev data",0,"trigger:K.Evt.isTrigger(trigList_2017_54AuAu)"},
     {"../../ntup/2017_54GeV_full/*11.root","54Gev data",0,"trigger: 1"},
     {"/media/petrchal/XSD-Linux/tpcAna/19GeV_May5/*1.root","19GeV data",0,"trigger: 1"}
   };

   TFileDescription data19GeV[]={
    {"../../ntup/2019_19GeV_Mar2023/*1.root","19GeV data",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //new data
    {"../../ntup/2019_19GeV_embed/kaon*.root","embed",1,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"},
     };

  //TFileDescription* files=data54GeV;
  TFileDescription* files=data19GeV;
  //TFileDescription* files=compare;
  const int nFiles=2; 
  const int order[]={0,1,2,1,2};
  const Long64_t nEntriefsLimit=TTree::kMaxEntries;//# or your number of entries

  //plotting modifiers
  int rebin=1;
  bool ignoreRange=false; // change to spot some outlayers
  const bool normalize=false; //plot normalized
 

    
  //loop over datasets
  for (int iFile=0;iFile<nFiles;iFile++){
 
    // structure for results
    ResultList1D Res_Plots; 


    TString fname=files[order[iFile]].fileName;
    cout<<"opening file: "<<fname<<endl;

    auto chain_kaons = new TChain("kaons");
    auto  fcount=chain_kaons->Add(fname,nEntriefsLimit);
    cout<<" TChain Added "<<fcount<<" files from "<<fname<<endl;
    //TObjArray * ll=chain_kaons->GetListOfFiles();
    //rott datafrme
        RDataFrame kaons_node(*chain_kaons); //raw event count
  
    //common event cuts
    auto  evCut=K3piCut_EventCut();
    if (files[order[iFile]].trigger){ //set propper trigers 
       evCut.Replace(files[order[iFile]].trigger);
    }
    auto d_events = kaons_node.Filter(evCut.Str());
    cout<<" Event cut used:  "<<endl <<evCut.Str()<<endl<<endl;
  
    //all 3pi candidates ... all together MC and non MC, positive and negative !!!
    //auto d_3piCandidates = d_events.Filter(K3piCut_VertexCanditate().Str());
     
    auto CurrentPos=Res_Plots.begin();
    

    //3pi with cut on mother PID from KFP (100321 for K+) both real and MC
    auto mother_cut=K3piCut_Mother_Kplus();
    if (files[order[iFile]].isMc) {mother_cut=K3piCut_Mother_Kplus() + Setup_MCvertex();}

    //matched kaons with cut on mother PID from KFP (100321 for K+) both real and MC
    auto kaon_cut=K3piCut_Matched_Kplus();
    if (files[order[iFile]].isMc) {kaon_cut=K3piCut_Matched_Kplus() + Setup_MCvertex();}
    
    
    // event plot per matched K
     auto tmppos=std::distance(Res_Plots.begin(),CurrentPos);
     AddPlots4QA(Event_plots,kaons_node,kaon_cut+evCut,Res_Plots,CurrentPos,"per matched K+",files[order[iFile]].lable,rebin,false);
     CurrentPos=Res_Plots.begin()+tmppos;
     AddPlots4QA(Event_plots,kaons_node,mother_cut+evCut,Res_Plots,CurrentPos,"per matched K+",files[order[iFile]].lable,rebin,false);
   
   
    //Mother(3pi vertex) per matched K
    tmppos=std::distance(Res_Plots.begin(),CurrentPos);
    AddPlots4QA(Mother_plots,d_events,kaon_cut+evCut,Res_Plots,CurrentPos,"per matched K+",files[order[iFile]].lable,rebin,false);
    CurrentPos=Res_Plots.begin()+tmppos;
    AddPlots4QA(Mother_plots,d_events,mother_cut+evCut,Res_Plots,CurrentPos,"per matched K+",files[order[iFile]].lable,rebin,false);
      
    //Properties of matched kaon
    tmppos=std::distance(Res_Plots.begin(),CurrentPos);
    AddPlots4QA(MatchedKaon_plots,d_events,kaon_cut+evCut,Res_Plots,CurrentPos,"of matched K+",files[order[iFile]].lable,rebin,false);
    CurrentPos=Res_Plots.begin()+tmppos;
    AddPlots4QA(MatchedKaon_plots,d_events,mother_cut+evCut,Res_Plots,CurrentPos,"of matched K+",files[order[iFile]].lable,rebin,false);
   

 // .. it is a problem, since histogram may come from different trees (nodes)
  //AddProgressBar(event_node);

   cout<<"trigger lazy evaluation"<<endl;
 
   auto ct= kaons_node.Count();
   
   cout<<*ct<<endl;

   cout<<"trigger DONE"<<endl;
 
  
   delete chain_kaons;

   DrawResults(Res_Plots);
 
} //loop over files


return;
}



