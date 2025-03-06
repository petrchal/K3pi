#include "TNtuple.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TEfficiency.h"
#include <iostream>
#include "K3piLib.C"
#include "K3piPlotDefs.C"
#include "data_Include.h"
#include "cuts_Include.h"


//int particlesPDG[nPIDS] = {100321, 200321,-100321, -200321}; //k->3pi only, K->3pi with found K


using namespace ROOT;
using RNode = ROOT::RDF::RNode;


  //


  const int   color[]={kBlack,kBlue,kRed,kGreen};
  //!!!!the _MC is important
/*
  //!!!! the first index is file
   RNode** (nodes[][2])={  //denumerator , numerator
                           {&node_cut_found3pi,&node_cut_K_helix},
                           {&node_cut_found3pi_MC,&node_cut_K_helix_MC},
                           {&node_cut_found3pi,&node_cut_K_helix},
                           {&node_cut_found3pi_MC,&node_cut_K_helix_MC},
                           {&node_cut_found3pi,&node_cut_K_helix},
                           {&node_cut_found3pi_MC,&node_cut_K_helix_MC},
                           {&node_cut_found3pi_MC,&node_cut_K_helix_MC},
                           };
                          
*/
    
  const Long64_t nEntriesLimit=TTree::kMaxEntries;;//1000000;//TTree::kMaxEntries;//100000;// -1;
  const int rebin=2; 
  bool ignoreRange=false; // change to spot some outlayers
  
  

  //std::vector<ROOT::RDF::Experimental::RResultMap<TH1D>> num;
  //std::vector<ROOT::RDF::Experimental::RResultMap<TH1D>> den;
    

  
  TEffList Res_Plots;
  

/*
  //TNtuple *nt;
  TChain *nt;
  TCanvas *c;

  TCanvas* cn[nVars]; 
  TCanvas* rcn[nVars]; 
 */
TCanvas *c;


int tmpVarCount;
void AddEffPlots(TPlotDefinitions& plotDefs, RNode numerator_node, RNode denom_node, std::vector<TEffFromSingleDef<TH1D>>& Res,const char * prefix);
void DrawEffs(std::vector<TEffFromSingleDef<TH1D>> &Res);
THStack *plotStackRatios(THStack * stack, bool plotNormalized=false, bool plotNormalizedRatio=false);
void plotEffStack(TEffFromSingleDef<TH1D> &res);

TH1D*  MakeRatioPlot(TH1D *num,TH1D *den);

void plotEfficiencies(){
  gROOT->Reset();
  c= new TCanvas();
  //const Float_t pi_m=0.494;
  ROOT::EnableImplicitMT(); //enambe multi threading - application must be MT safe ..your job
  TH1::SetDefaultSumw2(); //!!!

  //intilaize cuts
  InitCuts();
  Setup_3piVertexQA(); //must be called after InitCuts() !
  Setup_KaonMatching();

  TPlotDefinitions PlotDefs= Efficiency_plots;
  

  //loop over files
  for (int iFile=0;iFile<nFiles;iFile++){
    tmpVarCount=0; //reset the counter ..not realy necessary
    
   
    //TString TheCut=*cut[order[iFile]];
    TString fname=files[order[iFile]].fileName;
     cout<<"opening file: "<<fname<<endl;

    auto chain = new TChain("kaons");
    int nfiles=chain->Add(fname,nEntriefsLimit);
    cout<<" TChain Added "<<nFiles<<" nfiles.";
  

   // cout<<"Getentries="<<nt->GetEntries()<<endl;
   // delete nt;
   // continue;

    ROOT::RDF::RNode  data_node= ROOT::RDF::AsRNode(RDataFrame(*chain)); //K data
    
     //must be called after all cuts are read
     data_node=DefineNewVariables(data_node);
     //data_node=AddVariations(vary_EvtVz,data_node);
   
     
     data_node=AddVariations(vary_lastPointDiff,data_node);
     data_node=AddVariations(vary_3piVtx_chi2ndf,data_node);
     data_node=AddVariations(vary_dpDecay,data_node);
     data_node=AddVariations(vary_Minv,data_node);
     data_node=AddVariations(vary_daughter_Nhits,data_node);
     

    //start filtering event
    auto  evCut=K3piCut_EventCut();
    if (files[order[iFile]].trigger){ //set propper trigers 
       evCut.Replace(files[order[iFile]].trigger);
    } 
    auto node_Events=data_node.Filter(evCut.Str());
    cout<<" cut used for event selection pion:"<<endl<<evCut.Str()<<endl<<endl;


    //common denominator: any 3pi vertex ..this should speed it up
    //now 3pi candidates ... all together MC and non MC, positive and negative !!!
    //note - the previous cuts are actually repeated, likely does not matter much
    auto d_3piFound = node_Events.Filter(K3piCut_VertexCanditate().Str());
 
    
    //3pi+  - adding PID  
    auto _3pi_cut=K3piCut_3piVtx_Kplus();
    RNode node_cut_found3pi= d_3piFound.Filter(_3pi_cut.Str()); 
    cout<<" cut used for 3pi selection:"<<endl<<_3pi_cut.Str()<<endl<<endl;


    //3pi+ from MC
    auto _3pi_MC_cut= _3pi_cut+ Setup_MCvertex();
    RNode node_cut_found3pi_MC= node_cut_found3pi.Filter(_3pi_MC_cut.Str()); 

    //matched kaon
    auto matched_cut=K3piCut_Matched_Kplus();
    RNode node_cut_K_helix= node_cut_found3pi.Filter(matched_cut.Str()); 
    cout<<" cut used for matched kaon:"<<endl<<matched_cut.Str()<<endl<<endl;

    //matched from MC   
    //The qa truth is important!!! similar cut as in data
    auto matched_MC_cut=matched_cut+Setup_MCvertex();matched_MC_cut["qaT"]="(d.qaTruth[4]>50)";
    RNode node_cut_K_helix_MC=  node_cut_K_helix.Filter(matched_MC_cut.Str());
   


    
    if (files[order[iFile]].isMc)
      AddEffPlots(PlotDefs,node_cut_K_helix_MC,node_cut_found3pi_MC,Res_Plots,"");
    else
      AddEffPlots(PlotDefs,node_cut_K_helix,node_cut_found3pi,Res_Plots,"");
   

    /* does not work... not sure why
    cout<<"triggering evaluation"<<endl;
    if (Res_Plots.begin()!=Res_Plots.end()){
      //(Res_Plots.begin()->res.end()-1)->NumRMap()[0].DrawClone();
      //Res_Plots.begin()->res.back().NumRMap()[0].DrawClone();
      Res_Plots.begin()->num[0]["nominal"].DrawClone();
      }
   delete nt;
   */

    cout<<"trigger lazy evaluation"<<endl;
   auto ct= node_Events.Count();
   ct.OnPartialResult(/*every */100000/* events*/,
                           [](auto c) { std::cout << c << '\n'; });
   cout<<*ct<<endl;

   cout<<"trigger DONE"<<endl;

  } //file loop

#ifdef __EXTERNAl_OVERRIDE__

  TString side; 
   #ifdef _rightSIDE
    side="_rightVz";
  #endif
  #ifdef _leftSIDE
     side="_leftVz";
  #endif
  TString nm="eff_2019_AuAu19GeV_SL23";nm+=side;nm+="_DCA";nm+=c_DCA;nm+="_nhits";nm+=c_nhits;nm+=".root";
  #else
  
  TString nm="eff_2019_AuAu19GeV_SL24_rightHalf_10_DCA1_Nhits20.root";
 #endif

  cout<<"saving to"<<nm<<endl;
  TFile *f=new TFile(nm,"recreate");
  DrawEffs(Res_Plots);
  f->Write();
  f->Flush();
  //f->Close(); //to see after closing
  return;
}



//-----------------------------------------------------------------------
TH1D*  MakeRatioPlot(TH1D *num,TH1D *den){
   //cout<<"MakeRatioPlot"<<endl;
   TH1D* res= (TH1D*)num->Clone();
   TH1D* tmpden= (TH1D*)den->Clone();
 
 //for (int i=0;i<=tmpden->GetNbinsX();i++)
//    tmpden->SetBinError(i,0);  
 //      res->Divide(tmpden,"B");
   res->Divide(num,den,1.,1.,"B");
   delete tmpden;
   return res;
}


//-----------------------------------------------------------------------
void DrawEffs(std::vector<TEffFromSingleDef<TH1D>> &Res){
//draw 1D histograms in stack
gStyle->SetPadTickY(1);
gStyle->SetTickLength(0.02,"Y");
//gStyle->SetOptStat(0);
//

  /*
  TH1D* h_low;
  TH1D* h_hi;

  c = new TCanvas("c_tmp"); //to be able to call THNtuple->Draw(>>hist);
  THStack *s[nVars]; //efficiencies
  THStack *d[nVars]; //denominaotor
  THStack *n[nVars]; //numerator
  TLegend *l[nVars]; 
  TRatioPlot* r[nVars];
   */

 for ( auto &plot : Res){ //over vector of TEffFromSingleDef (based on single plot definiton)
    cout<<" plotting "<<plot.def.expr<<endl;
  
    TString nm="cn_";nm+=plot.def.title;
    c=new TCanvas(nm,nm,1200,600); 
    c->Divide(2,1);
    c->cd(1);
    auto l = new TLegend(0.1,0.7,0.48,0.9);
    l->SetHeader(plot.def.title,"C");
  
    nm="eff1_";nm+=plot.def.title;
    TCanvas* rc=new TCanvas(nm,nm,800,600); 

    //auto numMap=plot.NumRMap();
   // auto denMap=plot.NumRMap();
    
   // auto keys=plot.num.begin()->GetKeys(); //assuming all keys are the same!!!
    auto keys=plot.numMap.begin()->GetKeys(); //assuming all keys are the same!!!
    //that's actually not true sice matching cuts are only applied to numerator
    int nRes=plot.numMap.size(); //is is actually thenumber of input files at this points
    cout<<" I have "<<nRes<<" files for this plot"<<endl;
    int  nVariations =keys.size();
    cout<<" nVariations="<<nVariations<<endl;
    for (int iv=0;iv<nVariations;iv++){ //loop over variations
       
        int iii=0;

        auto currentVariation=new TEffResult<TH1D>(); 
        plot.res.push_back(currentVariation);
      
        //for ( auto &h : plot.res){ //over files (real data, mc)h is of type TEffResult
        for (int hi=0;hi<nRes;hi++){ 

        auto currentKey=keys[iv];
      
        cout<<"  var key= "<<currentKey<<endl;
 
         cout<<" files loop hi="<<hi<<endl;
         
         /*
         cout<<"chp1 "<< plot.num.size()<<endl;
         auto nummap=plot.num[hi];
         cout<<"chp1a "<< endl;
         cout<<"keysize="<< nummap.GetKeys().size()<<endl;
         nummap[currentKey].Draw();
         cout<<" aaaaa"<<endl;
         cout<<"entries="<<nummap[currentKey].GetEntries()<<endl;
          */
         /*
         auto keys2=plot.num[hi].GetKeys();
         int nk2=keys2.size();
         cout<<"test nkeys="<<nk2<<endl;
        for (int k2=0;k2<nk2;k2++) cout<<keys2[k2]<<endl;
        */
        
         TH1D* hNum=(TH1D*)plot.numMap[hi][currentKey].Clone();
         /*
         cout<<"chp1b"<<endl;
         auto keys2=plot.den[hi].GetKeys();
         int nk2=keys2.size();
         cout<<"test nkeys="<<nk2<<endl;
         for (int k2=0;k2<nk2;k2++) cout<<keys2[k2]<<endl;
          */

         //some variations are aplicable only numerator
         auto denKey=currentKey;
         auto keys2=plot.denMap[hi].GetKeys();
        if (std::find(keys2.begin(), keys2.end(), currentKey) == keys2.end()){
         cout<< " NO curret Key="<<currentKey<<" for denominator - using nominal"<<endl;
         denKey="nominal";
        }
        TH1D* hDen=(TH1D*)plot.denMap[hi][denKey].Clone();//(TH1D*) h.DenRMap()[currentKey].Clone();
         
         TString nm="denom of ";nm+=hDen->GetTitle();
         hDen->SetName(nm);
         hDen->SetTitle(nm);
         hDen->SetDirectory(0); //must be after the MakeRatioPlot
         hDen->SetLineColor(color[iii]);
         hDen->SetLineColor(kBlack);
         hDen->GetXaxis()->SetTitle(plot.def.axisTitle);
         nm="num of ";nm+=hNum->GetTitle();nm+=" "; nm+= currentKey;
         hNum->SetName(nm);
         hNum->SetTitle(nm);
         hNum->SetDirectory(0); //must be after the MakeRatioPlot
         hNum->SetLineColor(color[iii]);
         hNum->SetLineColor(color[iii]);
         hNum->GetXaxis()->SetTitle(plot.def.axisTitle);
        
        //make the efficiencies
        hNum->Draw();hDen->Draw();
       
        currentVariation->hNum.push_back(hNum);currentVariation->hDen.push_back(hDen);

        //  plot.res.push_back(h);
         auto ef=new TEfficiency(*hNum,*hDen);
         currentVariation->eff.push_back(ef);
         ef->SetStatisticOption(TEfficiency::kBBayesian);
         ef->SetConfidenceLevel(0.68);
      /* does not work
      auto c2 = new TCanvas();
      h.eff->Draw("EY");
      c2->Update(); 
     */
         ef->SetTitle(plot.def.title);
         //ef->SetLineColor(color[iii]);
         ef->SetLineColor(color[hi]);

         auto hr=MakeRatioPlot(hNum,hDen);
         TString r="eff:";nm+=hi;nm+=" of ";
         nm=hNum->GetTitle(); nm.ReplaceAll("num of ",r);
         hr->SetName(nm);
         hr->SetTitle(nm);
        
         hr->SetLineColor(color[hi]);
         currentVariation->ratio.push_back(hr);
         //assert(h.ratio);
      //h.ratio=(TH1D*)(h.eff->GetPaintedHistogram()->Clone());


     iii++;
     /*
    if (plot.stack_num==NULL){ //first file
            //create output stacks
            cout<<"  creating THStacks"<<endl;
            plot.stack_num=new THStack();//r.def.var,r.def.title);  ...name,title ..will do
            plot.stack_den=new THStack(); //"denumerator - found 3pi"
            plot.stack_ratio=new THStack(); 
        } //first file
*/
    //now add to stacks
    plot.stack_num.Add(hNum);
    plot.stack_den.Add(hDen);
    //plot.stack_ratio.Add(hr);
    //plot.stack_ratio->Add(h.eff->GetPaintedHistogram());
    
    //eff ratios
    //if (!first)
    //new TRatioPlot((TH1*)plot.res.begin()->eff,(TH1*) h.eff);
    //,opt);  --asi dulezite
   
  } //files loop
 
  }//variations loop

  cout<<" Variations loops done for "<< plot.def.expr<<endl;
  //r[v] = new TRatioPlot(s[v],(TH1*)hists->At(0),"gauss");
  //r[v]->Draw();
  //r[v]->GetLowerRefYaxis()->SetTitle("ratio");
  //r[v]->GetUpperRefYaxis()->SetTitle("efficiency");
  c->cd(1);
  //plot.stack_num->Draw("ehistnostack");
  plotStackRatios(&(plot.stack_num),true,true);
 
  c->cd(2);
  plotStackRatios(&(plot.stack_den),true,true);
  //l[v]->Draw();

            //hK->Draw("ehist");
            //l[v] = new TLegend(0.1,0.7,0.48,0.9);
            //l[v]->SetHeader(vars[v],"C");
  //c->Update();

  

  rc->cd();
  /*
  //must be done in loop since Tefficiency produces graph - cannot put to stack
  first=true; iii=0;
  for ( auto &h : plot.res) {
     TEfficiency *gEff=h.eff;
     if (first) h.ratio->Draw("hist"); else h.ratio->Draw("samehist");
     //gEff->Draw("A*"); else
     gEff->Draw("same0*");

     first=false;
     }  
  //assert(gEff);
  //cout<<gEff->GetTitle()<<endl;
  //for ( auto &h : plot.res) h.eff->Draw("A same");
 //plot stacks
 */

  plotEffStack(plot);
  rc->Write();


 // rc->Update();

 c->cd();
 c->Update();
 c->Write();
} //ploted variables

} //DrawEfs


//-----------------------------------------------------------
//plot hisograms with ratios ...i fth eratios do nt exist, create them
THStack *plotStackRatios(THStack * stack, bool plotNormalized, bool plotNormalizedRatio){
   cout<<" plotting stack ratio .."<< stack->GetNhists()<<endl;;
   Double_t pdiv = 0.3;
   auto topPad=gPad;
   TPad *pad1 = new TPad("p1", "p1", 0., pdiv, 1., 1.); // upper
   TPad *pad2 = new TPad("p2", "p2", 0., 0., 1., pdiv); // lower
   pad1->Draw();
   pad2->Draw();

   pad1->cd();

  //global constants with range of the plot
  //const Double_t gxmin = 17., gymin = 15., gxmax = 1990., gymax = 1000.;

 // TH1F* frame1 = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
  
//  frame1->GetXaxis()->SetMoreLogLabels();
  gPad->SetLeftMargin(0.09);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.11);//0.025
  gPad->SetBottomMargin(0.);

  pad2->SetLeftMargin(0.09);
  pad2->SetRightMargin(0.05);
  pad2->SetBottomMargin(0.2);//0.025
  pad2->SetTopMargin(0.);

  pad1->cd();
  float nlines=stack->GetNhists();
  auto l = new TLegend(0.55,0.75-nlines*0.05,0.9,0.85);
  l->SetTextSize(0.04);

  //actually this is a mem leak ---dont care for plotting
  //create stack of normalized histograms
   //normalize distributions for plotting
  THStack *norm=NULL;
  
  //normalized histograms
  if (plotNormalized || plotNormalizedRatio){
  TList *hlist=  stack->GetHists();
  norm=new THStack();
  for (int i=0;i<stack->GetNhists();i++){
        TH1 *h=(TH1*)hlist->At(i)->Clone();
        h->Scale(1./h->Integral());
        norm->Add(h);
        //l->AddEntry(h,label[order[i]],"l");
      }
  }

 //create ratios
 THStack *ratios=new THStack();
 //ratios from normalized histograms
 TList *nl;
 //ratios from normalized histograms
  if (plotNormalizedRatio)  nl= norm->GetHists();
    else nl= stack->GetHists();

  TH1 *hDen=(TH1*)nl->At(0);
  l->SetHeader(hDen->GetTitle(),"C");
  for (int i=1;i<nl->GetEntries();i++){
        TH1 *hNum=(TH1*)nl->At(i)->Clone();
        TString nm= hNum->GetName();
        nm="normrat: "+nm;
        hNum->SetName(nm);
        hNum->Divide(hDen);
        ratios->Add(hNum);
      }
  
  //the one to plot
  THStack *st;
  if( plotNormalized) st=norm; else st=stack;
  
  st->Draw("ehistnostack");
  st->SetMinimum(0.5);
  
  return NULL;

  auto frame1 = st->GetHistogram();
  Float_t siz = 0.045;
    
  frame1->SetTitleSize(siz);       frame1->SetLabelSize(siz);
  frame1->SetTitleSize(siz, "Y");  frame1->SetLabelSize(0.8*siz, "Y");
  frame1->SetLabelSize(0.00001, "X");
  frame1->GetXaxis()->SetTitleOffset(1.4);
  if (plotNormalized)frame1->GetYaxis()->SetTitle("normalized distribution      ");
  else frame1->GetYaxis()->SetTitle("counts      ");
  frame1->GetYaxis()->SetTitleOffset(1.); //0.64
  //frame1->GetYaxis()->SetTitleSize(0.05); 
  l->Draw("same");

  //plot normalized ratios
  pad2->cd();
  

  ratios->Draw("enostack");
  ratios->GetXaxis()->SetTitle(hDen->GetXaxis()->GetTitle());
  ratios->SetMinimum(0.5);
  ratios->SetMaximum(2);

  auto frame2 = ratios->GetHistogram();
  siz = siz*(1.+(1.-pdiv));
  frame2->SetTitleSize(siz);       frame2->SetLabelSize(siz);
  frame2->SetTitleSize(siz, "Y");  frame2->SetLabelSize(siz, "Y");
  frame2->GetXaxis()->SetTitleOffset(1.);
  frame2->GetXaxis()->SetTitleSize(0.1);
  if (plotNormalizedRatio) frame2->GetYaxis()->SetTitle("normalized ratio    "); else frame2->GetYaxis()->SetTitle("ratio    ");
  frame2->GetYaxis()->SetTitleOffset(0.4); //0.64
  frame2->GetYaxis()->SetTitleSize(0.1); 
  frame2->SetLabelSize(0.1, "X");
  frame2->SetLabelSize(0.08, "Y");
  
  
  topPad->Modified();

 return ratios;
}    

//---------------------------------------------------------------
void AddEffPlots(TPlotDefinitions& plotDefs, RNode numerator_node, RNode denom_node, TEffList& Res,const char * prefix)
 //std::vector<TEffFromSingleDef<TH1D>>::iterator &iter){
 {
    auto iter=Res.begin();

     for (auto def:plotDefs){

       //for fist file this means to create the structure for others only add histogram
        if (iter==Res.end()){ 
            // THstack *s=new THStack(plotDefs[v].var,nm); ...coul be done here
            TEffFromSingleDef<TH1D> r;
            r.def=def; 
            cout<<"pushing "<<def.expr<<endl;
            Res.push_back(r);
            iter=Res.end();--iter;
        }


        //if collumn name existed I woudl not need to do the Define.. but how to simply find out?
        TString var="tmpVar";var+=tmpVarCount++;
        // I could also use
        //it(iterator) - vec.begin()
        TString title=prefix; title+=def.title;title+=";"; title+=def.axisTitle;
        try{
         auto h_num=numerator_node.Histo1D(ROOT::RDF::TH1DModel(def.expr,title, 100/rebin, (ignoreRange)?0:def.lo, (ignoreRange)?0:def.hi),def.expr);
         iter->numMap.push_back(ROOT::RDF::Experimental::VariationsFor(h_num));
        }
        catch (const std::runtime_error &e) { //must Define the variable
         auto h_num=numerator_node.Define(var.Data(),def.expr) //this must be done for calculated variables
                    .Histo1D(ROOT::RDF::TH1DModel(def.expr,title, 100/rebin, (ignoreRange)?0:def.lo, (ignoreRange)?0:def.hi),var.Data()); 
         iter->numMap.push_back(ROOT::RDF::Experimental::VariationsFor(h_num));
       }

       try{
        auto h_den=denom_node.Histo1D(ROOT::RDF::TH1DModel(def.expr,title, 100/rebin, (ignoreRange)?0:def.lo, (ignoreRange)?0:def.hi),def.expr); 
        iter->denMap.push_back(ROOT::RDF::Experimental::VariationsFor(h_den));
        }
        catch (const std::runtime_error &e) { //must Define the variable
         auto h_den=denom_node.Define(var.Data(),def.expr) 
        .Histo1D(ROOT::RDF::TH1DModel(def.expr,title, 100/rebin, (ignoreRange)?0:def.lo, (ignoreRange)?0:def.hi),var.Data()); 
        iter->denMap.push_back(ROOT::RDF::Experimental::VariationsFor(h_den));
       }

       
        //cout<<" currently at "<<it1D->def.var<<" "<<it1D-Res1D.begin()<<endl;
        iter++;
  }
}


//-----------------------------------------------------------
//plot hisograms with ratios ...i fth eratios do nt exist, create them
void plotEffStack(TEffFromSingleDef<TH1D> &results){
   cout<<"plotting efficiences"<<endl;
   Double_t pdiv = 0.3;
   auto topPad=gPad;
   TPad *pad1 = new TPad("p1", "p1", 0., pdiv, 1., 1.); // upper
   TPad *pad2 = new TPad("p2", "p2", 0., 0., 1., pdiv); // lower
   pad1->Draw();
   pad2->Draw();

   pad1->cd();

  //global constants with range of the plot
  //const Double_t gxmin = 17., gymin = 15., gxmax = 1990., gymax = 1000.;

 // TH1F* frame1 = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
  
//  frame1->GetXaxis()->SetMoreLogLabels();
  gPad->SetLeftMargin(0.09);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.11);//0.025
  gPad->SetBottomMargin(0.);

  pad2->SetLeftMargin(0.09);
  pad2->SetRightMargin(0.05);
  pad2->SetBottomMargin(0.2);//0.025
  pad2->SetTopMargin(0.);
  pad2->SetGridy();

  pad1->cd();
  
 /* 
  THStack *stack=&(res.stack_ratio);
  float nlines=stack->GetNhists();
  auto l = new TLegend(0.55,0.75-nlines*0.05,0.9,0.85);
  l->SetHeader(res.def.title,"C");
  l->SetTextSize(0.04);

  //actually this is a mem leak ---dont care for plotting
  //create stack of normalized histograms
   //normalize distributions for plotting
  THStack *norm=NULL;
  
  TList *hlist= stack->GetHists();

  norm=new THStack();
  for (int i=0;i<stack->GetNhists();i++){
          TH1 *h=(TH1*)hlist->At(i);
       //   l->AddEntry(h,label[order[i]],"l");
      }
  /* 
  //create ratios
  //ratios of efficiencies
  THStack *ratios=new THStack();
  TH1 *hDen=(TH1*)hlist->At(0);
  for (int i=1;i<stack->GetNhists();i++){
        TH1 *hNum=(TH1*)hlist->At(i)->Clone();
        hNum->Divide(hDen);
        ratios->Add(hNum);
      }
  */

 THStack *AllEffs=new THStack();

  

  //create ratios of efficiencies for each variation
 cout<<" results.res.size()="<<results.res.size()<<endl;
  for (auto effres : results.res) { //loop over variations
    int nSets=effres->ratio.size();
    cout<<"nSets ="<<nSets<<endl;
    if (nSets<=1) continue;
    TH1 *hDen=effres->ratio[0];
    AllEffs->Add(hDen);
    for (int i=1;i<nSets;i++){ //loop over datasets(file) - usually MC and data
        TH1 *hNum=(TH1*)effres->ratio[i]->Clone();
        AllEffs->Add((TH1*)hNum->Clone());
        TString nm= hNum->GetName();
        nm="rat: "+nm;
        hNum->SetName(nm);
        hNum->Divide(hDen);
        results.stack_ratio.Add(hNum);
      }
    }
  
 
  AllEffs->Draw("ehistnostack");
  AllEffs->SetMinimum(0.00001);
  
  auto frame1 = AllEffs->GetHistogram();
  Float_t siz = 0.045;
    
  frame1->SetTitleSize(siz);       frame1->SetLabelSize(siz);
  frame1->SetTitleSize(siz, "Y");  frame1->SetLabelSize(0.8*siz, "Y");
  frame1->SetLabelSize(0.00001, "X");
  frame1->GetXaxis()->SetTitleOffset(1.4);
  frame1->GetYaxis()->SetTitle("efficiency      ");
 frame1->GetYaxis()->SetTitleOffset(1.); //0.64
  //frame1->GetYaxis()->SetTitleSize(0.05); 

 
  auto l = new TLegend(0.65,0.75,0.9,0.9);
  l->SetHeader(results.def.title,"C");
  auto &effres=results.res[0]; //first of variation
  for (int i=0;i<nFiles;i++)l->AddEntry(effres->ratio[i],files[order[i]].lable,"l");
    l->Draw("same");


  //plot normalized ratios
  pad2->cd();

  results.stack_ratio.Draw("enostack");
  results.stack_ratio.GetXaxis()->SetTitle(results.res[0]->ratio[0]->GetXaxis()->GetTitle());
  results.stack_ratio.SetMinimum(0.5);
  results.stack_ratio.SetMaximum(2);

  auto frame2 = results.stack_ratio.GetHistogram();
  siz = siz*(1.+(1.-pdiv));
  frame2->SetTitleSize(siz);       frame2->SetLabelSize(siz);
  frame2->SetTitleSize(siz, "Y");  frame2->SetLabelSize(siz, "Y");
  frame2->GetXaxis()->SetTitleOffset(1.);
  frame2->GetXaxis()->SetTitleSize(0.1);
  frame2->GetYaxis()->SetTitle("ratio    ");
  frame2->GetYaxis()->SetTitleOffset(0.4); //0.64
  frame2->GetYaxis()->SetTitleSize(0.1); 
  frame2->SetLabelSize(0.1, "X");
  frame2->SetLabelSize(0.08, "Y");
  

  
  topPad->Modified();


// return ratios;
}    
     /*
      else{
        cn[v]->cd();
        //hK->Draw("samehiste");
      }
      l[v]->AddEntry(hK,label[plot[i]],"l");
      s[v]->Add(hK);
      hDen->Scale(1./hDen->Integral());
      hNum->Scale(1./hNum->Integral());
      d[v]->Add(hDen);
      n[v]->Add(hNum);
      //n[v]->Add(res);
    }
    */

    /*
    ef[i]=new TEfficiency(*hK,*h3pi);
    ef[i]->SetDirectory(0);
    TString cn="cn";cn+=i;
    new TCanvas(cn,cn);
    ef[i]->Draw();  
    */
        /*
    TH1D* h=(TH1D*)hK->Clone();
    TString nm="hk";nm+=i;
    h->SetName(nm);
    h->Divide(h3pi);
    h->SetDirectory(0);
    h->SetLineColor(color[i]);
    c1->cd();
    if (i==0){
      h->SetTitle("\"K+3pi\" / \"3pi only\"");
      h->Draw("");
    }
    else h->Draw("same");
    h->SetLineColor(color[plot[i]]);
    l1->AddEntry(h,label[plot[i]],"l");
   */
  

/*
for (int v=0;v<nVars;v++){
  cn[v]->Divide(2,1);
  cn[v]->cd(1);
  d[v]->Draw("ehistnostack");
  cn[v]->cd(2);
  n[v]->Draw("ehistnostack");
  l[v]->Draw();
  TList *hists =s[v]->GetHists();
  TString opt="pois";
  //TString opt="diffsig";
  
  rcn[v]=new TCanvas(); 
 
  r[v] = new TRatioPlot((TH1*)hists->At(0),(TH1*) hists->At(1),opt);
  //r[v] = new TRatioPlot(s[v],(TH1*)hists->At(0),"gauss");
  r[v]->Draw();
  r[v]->GetLowerRefYaxis()->SetTitle("ratio");
  r[v]->GetUpperRefYaxis()->SetTitle("efficiency");
  //*/

  //hack
  //s[v]->Draw("nostack,e1p"); 




/*
void SaveResults(const char* dir){
  TString nDir="pics/";nDir+=dir;
  gSystem->mkdir(nDir,1);
  
 for (int v=0;v<nVars;v++){
  TString nm=nDir;
  nm+="/";nm+=vars[v];
  TString nm_=nm+".eps";
  rcn[v]->SaveAs(nm_);
  nm_=nm+".C";
  rcn[v]->SaveAs(nm_);
  nm=nDir;
  nm+="/";nm+=vars[v];
  nm_=nm+"_dist.eps";
  cn[v]->SaveAs(nm_);
  nm_=nm+"_dist.C";
  cn[v]->SaveAs(nm_);
 }
}
*/