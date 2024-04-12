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
  ResultList1D Res_Plots; 
   
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
    cout<<endl<<" K+ 3piVtx_ cut used:  "<<endl<<  3piVtx__cut.Str()<<endl<<endl;
 
    // event plot per 3pi+
    auto Cut=Reco3piVtx_cut+evCut;
    //AddPlots4QA(Event_plots,kaons_node,ev_cut,Res_Plots,CurrentPos,"per found 3pi+",files[order[iFile]].lable,rebin,false);
   
  
     //d_events=AddVariations(vary_EvtVz,d_events);
     /*
     if (files[order[iFile]].trigger){ //set propper trigers 
       Cut.Replace(files[order[iFile]].trigger);
    }
    */
    cout<<endl<<"Cut used:  "<<endl<<  Cut.Str()<<endl<<endl;
 

    auto CurrentPos=Res_Plots.begin();
    
    //Mother(3pi vertex) per found 3pi vertex 
    //AddPlots4QA(Mother_plots,d_events,mother_cut,Res_Plots,CurrentPos,"per found 3pi+",files[order[iFile]].lable,rebin,false);
      AddPlots4QA(RecoVtx_plots,d_events,Cut,Res_Plots,CurrentPos,"per found 3pi+",files[order[iFile]].lable,rebin,false);
      //!!! this is without kaon matching cuts
      //AddPlots4QA(MatchedKaon_plots,d_events,Cut,Res_Plots,CurrentPos,"per found 3pi+",files[order[iFile]].lable,rebin,false);
      
   
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



