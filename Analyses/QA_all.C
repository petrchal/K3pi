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


void QA_all(){
  gROOT->Reset();
  ROOT::EnableImplicitMT(); //enambe multi threading - application must be MT safe ..your job
  InitCuts();
  //here must be called cuts that cannot be called in InitCuts in roder to define variables
  Setup_3piVertexQA();

  //plotting modifiers
  int rebin=1;
  bool ignoreRange=false; // change to spot some outlayers, that ratios then do not make sense...
  //const bool normalize=false; //plot normalized
 

  // structure for results
  ResultList1D Res_events; 
  ResultList1D Res_Vtx; 
  ResultList1D Res_kaons; 
  //matched kaons
  ResultList1D Res_events_K; 
  ResultList1D Res_Vtx_K; 
  ResultList1D Res_kaons_K; 

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
    if (files[order[iFile]].trigger){ //set propper trigers and modifications
       evCut.Replace(files[order[iFile]].trigger);
    }

    //must be called after all cuts are read
    //defines varibles that can be later varied for systematics
    kaons_node=DefineNewVariables(kaons_node);
     //kaons_node=AddVariations(vary_EvtVz,kaons_node);
     //kaons_node=AddVariations(vary_EvtVz,kaons_node);
     //kaons_node=AddVariations(vary_lastPointDiff,kaons_node);
   


    ROOT::RDF::RNode d_events = kaons_node.Filter(evCut.Str());
    cout<<" Event cut used:  "<<endl <<evCut.Str()<<endl<<endl;
     //all 3pi candidates ... all together MC and non MC, positive and negative !!!
    //auto d_3piCandidates = d_events.Filter(K3piCut_VertexCanditate().Str());
    
  
    //3pi with cut on  PID from KFP (100321 for K+) both real and MC
    auto Reco3piVtx_cut=K3piCut_3piVtx_Kplus();
    //for MC data add condition on MC vertex.. 
    if (files[order[iFile]].isMc) {Reco3piVtx_cut=K3piCut_3piVtx_Kplus() + Setup_MCvertex();}
    //cout<<endl<<" K+ 3piVtx_ cut used:  "<<endl<<  3piVtx__cut.Str()<<endl<<endl;
 
    // event plot per found 3pi+  ... this is different cthen QA_events
    auto Cut=Reco3piVtx_cut+evCut;
    auto CurrentPos=Res_events.begin();
    AddPlots4QA(Event_plots,kaons_node,Cut,Res_events,CurrentPos,"per found 3pi+",files[order[iFile]].lable,rebin,false);
   
  
     //d_events=AddVariations(vary_EvtVz,d_events);
     /*
     if (files[order[iFile]].trigger){ //set propper trigers 
       Cut.Replace(files[order[iFile]].trigger);
    }
    */
    cout<<endl<<"Cut used:  "<<endl<<  Cut.Str()<<endl<<endl;
 

     
    
      CurrentPos=Res_Vtx.begin();
      AddPlots4QA(RecoVtx_plots,d_events,Reco3piVtx_cut,Res_Vtx,CurrentPos,"per found 3pi+",files[order[iFile]].lable,rebin,false);
      
      //!!! this is without kaon matching cuts
      //ROOT::RDF::RNode d_Vtx = kaons_node.Filter(Reco3piVtx_cut.Str());
      CurrentPos=Res_kaons.begin();
      AddPlots4QA(MatchedKaon_plots,d_events,Reco3piVtx_cut,Res_kaons,CurrentPos,"per found 3pi+",files[order[iFile]].lable,rebin,false);
     
    //---------------now matched kaons-------------------------
   
  
    //matched kaons with cut on mother PID from KFP (100321 for K+) both real and MC
    auto MatchedKaon_cut=K3piCut_Matched_Kplus();
    //for MC filed I only add cut on MC vertex !!!!not on the KAON - allows to study mis matches
    if (files[order[iFile]].isMc) {MatchedKaon_cut=K3piCut_Matched_Kplus() + Setup_MCvertex();}
    cout<<endl<<" K+ Mother cut used:  "<<endl<< MatchedKaon_cut.Str()<<endl<<endl;
   
    Cut=evCut+MatchedKaon_cut;
    CurrentPos=Res_events_K.begin();
    AddPlots4QA(Event_plots,kaons_node,Cut,Res_events_K,CurrentPos,"per matched K+",files[order[iFile]].lable,rebin,false);

    CurrentPos=Res_Vtx_K.begin();
    AddPlots4QA(RecoVtx_plots,d_events,MatchedKaon_cut,Res_Vtx_K,CurrentPos,"per matched K+",files[order[iFile]].lable,rebin,false);
    
    CurrentPos=Res_kaons_K.begin();
    AddPlots4QA(MatchedKaon_plots,d_events,MatchedKaon_cut,Res_kaons_K,CurrentPos,"per matched K+",files[order[iFile]].lable,rebin,false);


   cout<<"trigger lazy evaluation"<<endl;
 
   auto ct= kaons_node.Count();
   ct.OnPartialResult(/*every */100000/* events*/,
                           [](auto c) { std::cout << c << '\n'; });
   cout<<*ct<<endl;

   cout<<"trigger DONE"<<endl;
 
  
   delete chain_kaons;
 
} //loop over files

TFile f("QA_dataComp_ZDC:200:600.root","recreate");

f.mkdir("3pi");
f.mkdir("3pi/event info per found 3piVtx");f.cd("3pi/event info per found 3piVtx");
DrawResults(Res_events);
f.mkdir("3pi/3piVtxs");f.cd("3pi/3piVtxs");
DrawResults(Res_Vtx);
f.mkdir("3pi/matched kaons without cuts");f.cd("3pi/matched kaons without cuts");
DrawResults(Res_kaons);

f.mkdir("matched Kaons");
f.mkdir("matched Kaons/event per matched K");f.cd("matched Kaons/event per matched K");
DrawResults(Res_events_K);
f.mkdir("matched Kaons/3piVtxs");f.cd("matched Kaons/3piVtxs");
DrawResults(Res_Vtx_K);
f.mkdir("matched Kaons/kaons");f.cd("matched Kaons/kaons");
DrawResults(Res_kaons_K);


f.Write();
f.Close();

return;
}



