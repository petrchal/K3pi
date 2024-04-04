#include "K3piLib.C"

//this is mean as global setup to make sure that I use alway the same files and cut
//one can override in the script

  //2017 54 data 
  TFileDescription data54GeV[]={
    {"../../ntup/2017_54GeV_full/*11.root","data",0},
   };


 //27 data vs embeddings
  TFileDescription data27GeV[]={
    {"../../ntup/2018_27GeV_Sept22/*.root","27GeV data",0,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"}, //new data
    {"../../ntup/2018_27GeV_2023Jan/*.root","27GeV data rerun",0,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"},
    {"../../ntup/2018_27GeV_embed_ptW_1M/kaon*.root","embed",0},
    {"../../ntup/2018_27GeV_embed_ptW_1M/kaon*.root","embed isMC",1},
    {"../../ntup/2018_embed_flatPt_forcedDec/kaons.root","MC flat Pt",1}
    }; 


    //compare runs
    TFileDescription compare[]={
     {"../../ntup/2018_27GeV_Sept22/*.root","27GeV data",0,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"}, //new data
     {"../../ntup/2018_27GeV_Sept22/*1.root","27GeV data no M_adds",0,"trigger:Evt.isTrigger(trigList_2018_27AuAu),MotherAdds:1"}, //new data
     {"../../ntup/2018_27GeV_2023Jan/*1.root","27GeV data rerun",0,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"},
    //{"/media/petrchal/XSD-Linux/tpcAna/27GeV_Apr22/*1.root","old data",0},
      {"../../ntup/2018_27Embed_tuned/kaon*.root","embedding ",1,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"},
     {"../../ntup/2018_27Embed_tuned/kaon*.root","tuned embed no M_adds",1,"trigger:Evt.isTrigger(trigList_2018_27AuAu),MotherAdds:1"},
      //{"../../ntup/2018_27Embed_tuned/kaon*.root","tuned embed without hit rat`",1,"trigger:Evt.isTrigger(trigList_2018_27AuAu),nhits_posrat:1"},
      
      //   {"../../ntup/2018_27GeV_embed_ptW_1M/kaon*.root","old embedding",1,"trigger:K.Evt.isTrigger(trigList_2018_27AuAu)"},
      //  {"../../ntup/2018_27GeV_embed_tuned/kaon*.root","old embedding",1,"trigger:K.Evt.isTrigger(trigList_2018_27AuAu)"},
        // {"../../ntup/2017_54GeV_full/*1.root","54Gev data",0,"trigger:K.Evt.isTrigger(trigList_2017_54AuAu)"},
        //{"../../ntup/2017_54GeV_full/*.root","54Gev data",0,"trigger: 1"},
        //{"/media/petrchal/XSD-Linux/tpcAna/19GeV_May5/*1.root","19GeV data",0,"trigger: 1"}
   };

  TFileDescription data19GeV[]={
   //SL21 data - older run with all events
   // {"../../ntup/2019_19GeV_Mar2023/*.root","19GeV data SL21 nM=4",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu),nK3piP:Evt.nK3piP==4"}, //March - new data withouut XY DCA
    {"../../ntup/2019_19GeV_Mar2023/*.root","data SL21",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //March - new data withouut XY DCA
   //{"../../ntup/2019_19GeV_Mar2023/*.root","19GeV data SL21 lowZDC",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu),ZDC:(Evt.ZDCx>100) && (Evt.ZDCx<300)"}, //March - new data withouut XY DCA
   //{"../../ntup/2019_19GeV_Mar2023/*.root","19GeV data SL21  highZDC",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu),ZDC:(Evt.ZDCx>400) && (Evt.ZDCx<600)"}, //March - new data withouut XY DCA
 
   
    //SL21 data - new but the data have lower stat 
    {"../../ntup/2019_19GeV_SL21/*.root","data SL21 ",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"},
    //{"../../ntup/2019_19GeV_SL21/*.root","19GeV data SL21 nK=4",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu),nK3piP:Evt.nK3piP==4"},
   
    //SL21 embed
    {"../../ntup/2019_19GeV_SL21_emb/*.root","embedding SL21 ",1,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"},
    //{"../../ntup/2019_19GeV_SL21_emb/*.root","embedding SL21 nK=4",1,"trigger:Evt.isTrigger(trigList_2019_19AuAu),nK3piP:Evt.nK3piP==4"},

    
   //SL23 data 
   {"../../ntup/2019_19GeV_SL23/*.root","data SL23 rerun",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //new data
   // {"../../ntup/2019_19GeV_SL23/*.root","19GeV data SL23 rerun nK=4",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu),nK3piP:Evt.nK3piP==4"}, //new data
   //{"../../ntup/2019_19GeV_SL23/*.root","19GeV data SL23 lowZDC",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu),ZDC:(Evt.ZDCx>100) && (Evt.ZDCx<300)"}, //new data
   //{"../../ntup/2019_19GeV_SL23/*.root","19GeV data SL23 hiZDC",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu),ZDC:(Evt.ZDCx>400) && (Evt.ZDCx<600)"}, //new data
    
    //SL23 embed   
    {"../../ntup/2019_19GeV_SL23_embedding/*.root","embedding SL23",1,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //new data
    {"../../ntup/2019_19GeV_SL23_embedding/*.root","embedding SL23 nK=3",1,"nK3piP:Evt.nK3piP==3"}, //new data
   };

   TFileDescription data_2020_FXT[]={
    {"../../ntup/2020_FXT_31p2_P23id/*.root","31.2GeV FXT SL23d",0,NULL}, 
    {"../../ntup/2020_FXT_13p5_P23ie/*.root","13.5GeV FXT SL23e",0,NULL}, 
    {"../../ntup/2020_FXT_5p75_P23ie/*.root","5.75GeV FXT Sl23e",0,NULL}, 
    }; 


   TFileDescription data_2021_7p7[]={
    {"../../ntup/2021_7p7_official/*.root","7.7GeV official data",0,"trigger:Evt.isTrigger(trigList_2021_7p7AuAu)"}, 
    {"../../ntup/kaons_TFGsim_hiMult.root","TFG sim",1,"trigger:1,VPDdif:1"},
    {"../../ntup/2021_7p7_TFG/*.root","7.7GeV TFG data",0,"trigger:Evt.isTrigger(trigList_2021_7p7AuAu)"}, 
    };    

  //TFileDescription* files=data54GeV;
  //TFileDescription* files=data27GeV;
  //TFileDescription* files=data_2021_7p7;
 TFileDescription* files=data19GeV;
  //TFileDescription* files=data_2020_FXT;
// K3piCut_EventCut =EventCut_2021_7p7AuAu; //just in case override 
   //to compare to TFG
 
  const int nFiles=2; 
  const int order[]={1,3};//{1,3};//{3,5};//{1,9};//{2,0,4};//{2,4,3,3};
  const Long64_t nEntriefsLimit=TTree::kMaxEntries;//# or your number of entrie


// I want calculate fourier transform




