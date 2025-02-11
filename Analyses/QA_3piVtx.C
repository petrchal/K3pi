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


void QA_3piVtx(){
  gROOT->Reset();
  ROOT::EnableImplicitMT(); //enambe multi threading - application must be MT safe ..your job
  InitCuts();
  //here must be called cuts that cannot be called in InitCuts in roder to define variables
  Setup_3piVertexQA();
  
  //plotting modifiers
  int rebin=1;
  bool ignoreRange=false; // change to spot some outlayers
  //const bool normalize=false; //plot normalized
 

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
    auto  fcount=chain_kaons->Add(fname,nEntriefsLimit);
    cout<<" TChain Added "<<fcount<<" files from "<<fname<<endl;
    //TObjArray * ll=chain_kaons->GetListOfFiles();

    //root dataframe
    ROOT::RDF::RNode  kaons_node= ROOT::RDF::AsRNode(RDataFrame(*chain_kaons)); //raw event count
  
    //common event cuts
    auto  evCut=K3piCut_EventCut();
    if (files[order[iFile]].trigger){ //set propper trigers 
       evCut.Replace(files[order[iFile]].trigger);
    }

    //must be called after all cuts are read
    kaons_node=DefineNewVariables(kaons_node);

    ROOT::RDF::RNode d_events = kaons_node.Filter(evCut.Str());
    cout<<" Event cut used:  "<<endl <<evCut.Str()<<endl<<endl;
     //all 3pi candidates ... all together MC and non MC, positive and negative !!!
    //auto d_3piCandidates = d_events.Filter(K3piCut_VertexCanditate().Str());
    
      
    //3pi with cut on  PID from KFP (100321 for K+) both real and MC
    auto Reco3piVtx_cut=K3piCut_3piVtx_Kplus();
    //for MC data add condition on MC vertex.. 
    if (files[order[iFile]].isMc) {Reco3piVtx_cut=K3piCut_3piVtx_Kplus() + Setup_MCvertex();}
    cout<<endl<<" K+ 3piVtx_ cut used:  "<<endl<<  Reco3piVtx_cut.Str()<<endl<<endl;
 
    Res_EventPlots.resetPosition();
    Res_EventPlots_2D.resetPosition(); 
    Res_3piPlots.resetPosition();
    Res_3piPlots_2D.resetPosition(); 
    Res_KaonPlots.resetPosition(); 
    Res_KaonPlots_2D.resetPosition(); 
   

    // event plot per 3pi+
    auto Cut=Reco3piVtx_cut+evCut;
    AddPlots4QA(Event_plots,kaons_node,Cut,Res_EventPlots,"per found 3pi+",files[order[iFile]].lable,rebin,false);
    AddPlots4QA(Event_plots_2D,kaons_node,Cut,Res_EventPlots_2D,"per found 3pi+",files[order[iFile]].lable,rebin,false);
   
  
     //d_events=AddVariations(vary_EvtVz,d_events);
     /*
     if (files[order[iFile]].trigger){ //set propper trigers 
       Cut.Replace(files[order[iFile]].trigger);
    }
    */
    cout<<endl<<"Cut used:  "<<endl<<  Cut.Str()<<endl<<endl;
 

    Cut=Reco3piVtx_cut;
    //Mother(3pi vertex) per found 3pi vertex 
    AddPlots4QA(RecoVtx_plots,d_events,Cut,Res_3piPlots,"per found 3pi+",files[order[iFile]].lable,rebin,false);
    AddPlots4QA(RecoVtx_plots_2D,d_events,Cut,Res_3piPlots_2D,"per found 3pi+",files[order[iFile]].lable,rebin,false);
    //!!! this is without kaon matching cuts
    AddPlots4QA(Kaon_plots,d_events,Cut,Res_KaonPlots,"per found 3pi+",files[order[iFile]].lable,rebin,false);
    AddPlots4QA(Kaon_plots_2D,d_events,Cut,Res_KaonPlots_2D,"per found 3pi+",files[order[iFile]].lable,rebin,false);
     
   
   cout<<"trigger lazy evaluation"<<endl;
   auto ct= kaons_node.Count();
   ct.OnPartialResult(/*every */100000/* events*/,
                           [](auto c) { std::cout << c << '\n'; });
   cout<<*ct<<endl;
   cout<<"trigger DONE"<<endl;
 
  
   delete chain_kaons;
 
} //loop over files

TFile *f=new TFile("3piComp_2019_SL24_noCuts.root","recreate");
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



