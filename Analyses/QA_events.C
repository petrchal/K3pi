#include "TNtuple.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEfficiency.h"
#include <iostream>
#include "K3piLib.C"
#include "K3piPlotDefs.C"
#include "cuts_Include.h"
#include "data_Include.h"


using namespace ROOT;

 //int particlesPDG[nPIDS] = {100321, 200321,-100321, -200321}; //k->3pi only, K->3pi with found K

 //the used cut  - all include event cuts
  //cut_found3pi - any 3pi vertex - real or MC or mixed;
  //cut_found3pi_MC;//cut_K_helix;//cut_K_KF;//cut_found3pi;//cut_K_helix;
 


 

void QA_events(){
  gROOT->Reset();
  ROOT::EnableImplicitMT(); //enambe multi threading - application must be MT safe ..your job
 
  //initilaize variables
  InitCuts(false);

 //K3piCut_EventCut =EventCut_2019_19AuAu; //just in case override 

  //plotting modifiers
  int rebin=1;
  bool ignoreRange=false; // change to spot some outlayers
  gIgnoreCutMods=false;
  //const bool normalize=false; 
 
 
  // strucuture for results
  ResultList1D Res_eventPlots; 
  ResultList2D Res_eventPlots_2D; 
   
  //loop over datasets
  for (int iFile=0;iFile<nFiles;iFile++){
   
    TString fname=files[order[iFile]].fileName;
    cout<<"opening file: "<<fname<<endl;     

    //event count trees
   // auto chain_events = new TChain("events");
    auto chain_events = new TChain("events");
    auto  fcount=chain_events->Add(fname,nEntriesLimit);
    cout<<" TChain Added "<<fcount<<" files from "<<fname<<endl;
   //TObjArray * ll=chain_kaons->GetListOfFiles();

   //use tchain as input to data frame
   // RDataFrame event_node(*chain_events); //raw event count
   ROOT::RDF::RNode event_node = ROOT::RDF::AsRNode(RDataFrame(*chain_events)); //raw event count

   Res_eventPlots.resetPosition();
   Res_eventPlots_2D.resetPosition();
   
   
   
    auto evCut=K3piCut_EventCut(); 

      if (files[order[iFile]].trigger){ //set propper trigers 
        evCut.Replace(files[order[iFile]].trigger);
     }
   cout<<" Event cut used:  "<<endl <<evCut.Str()<<endl<<endl;
   
    event_node=DefineNewVariables(event_node);
   //  event_node=AddVariations(vary_EvtVz,event_node);
    AddPlots4QA(Event_plots,event_node,evCut,Res_eventPlots,"events",files[order[iFile]].lable,rebin,false);
    AddPlots4QA(Event_plots_2D,event_node,evCut,Res_eventPlots_2D,"events",files[order[iFile]].lable,rebin,false);
   
    
    cout<<"trigger lazy evaluation"<<endl;
    auto ct= event_node.Count();
    ct.OnPartialResult(/*every */100000/* events*/,
                           [](auto c) { std::cout << c << '\n'; });
    cout<<*ct<<endl;
    cout<<"trigger DONE"<<endl;
 
  
   delete chain_events; 
 
} //loop over files

TFile *f=new TFile("eventsComp_2019_19GeV_SL24vsSL24.root","recreate");
DrawResults(Res_eventPlots);
DrawResults(Res_eventPlots_2D);
f->Write();
f->Flush();
//f->Close(); //don't close to see resutls

return;
}



