#include "K3piLib.C"

//this is mean as global setup to make sure that I use alway the same files and cut
//one can override in the script

  //2017 54 data 
  TFileDescription data54GeV[]={
    {"../../ntup/2017_54GeV_full/*11.root","data",0},
   };


 //27 data vs embeddings
  TFileDescription data27GeV[]={
    {"/media/petrchal/XSD-Linux/tpcAna/2018_27GeV_Sept22/*.root","27GeV data",0,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"}, //new data
    //{"../../ntup/2018_27GeV_Sept22/*.root","27GeV data",0,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"}, //new data
    //{"../../ntup/2018_27GeV_2023Jan/*.root","27GeV data rerun",0,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"},
    {"/media/petrchal/XSD-Linux/tpcAna/2018_27GeV_embed_ptW_1M/kaon*.root","embed",1,"trigger:Evt.isTrigger(trigList_2018_27AuAu)"},
    //{"../../ntup/2018_27GeV_embed_ptW_1M/kaon*.root","embed isMC",1},
    //{"../../ntup/2018_embed_flatPt_forcedDec/kaons.root","MC flat Pt",1}
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
    //{"../../ntup/2019_19GeV_Mar2023/*.root","data SL21",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //March - new data withouut XY DCA
  
   
    //SL21 data - new but the data have lower stat 
    //{"../../ntup/2019_19GeV_SL21/*.root","data SL21 ",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"},
   
    //SL21 embed
  //  {"../../ntup/2019_19GeV_SL21_emb/*.root","embedding SL21 ",1,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"},

    
   //SL23 data 
    {"../../ntup/2019_19GeV_SL23/*.root","data SL23 ",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //new data
    //{"../../ntup/2019_19GeV_SL23/*.root","data SL23  with runId selection",0,"trigger:Evt.isTrigger(trigList_2019_19AuAu),runId:goodRunId(Evt.runId)"}, 
    
    //SL23 embed  - final full stat
    {"../../ntup/2019_19GeV_SL23_embed_final/*.root","embedding SL23",1,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //new data
    {"../../ntup/2019_19GeV_SL23_embed_final_devKFP/*.root","embedding SL23 under devKFP",1,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //new data, KFP under dev
    

    //SL24y
    {"../../ntup/2019_19GeV_SL24y_FullRerunJuly2025_devKFP/*.root","data SL24y devKFP TOFmatch_low",0,  "trigger:Evt.isTrigger(trigList_2019_19AuAu),TOF:Evt.nBTOFMatch<100"}, //new data
    //{"../../ntup/2019_19GeV_SL24y_FullRerunJuly2025_devKFP/*.root","data SL24y devKFP TOFmatch_low",0,  "trigger:Evt.isTrigger(trigList_2019_19AuAu),TOF:(Evt.nBTOFMatch>100)&&(Evt.nBTOFMatch<200)"}, //new data
    //{"../../ntup/2019_19GeV_SL24y_FullRerunJuly2025_devKFP/*.root","data SL24y devKFP TOFmatch_hi",0,  "trigger:Evt.isTrigger(trigList_2019_19AuAu),TOF:(Evt.nBTOFMatch>300)&&(Evt.nBTOFMatch<400)"}, //new data
    //{"../../ntup/2019_19GeV_SL24y_FullRerunJuly2025_devKFP/*.root","data SL24y devKFP nK3piP>0",0,  "nK3piP:Evt.nK3piP>0,trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //new data
    // {"../../ntup/2019_19GeV_SL24y_FullRerunJuly2025_devKFP/*.root","data SL24y devKFP nK3piP==1",0,"nK3piP:Evt.nK3piP==1,trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //new data
    //embed with SL24y  
    {"/home/petrchal/tpcAna/ntup/2019_19GeV_SL24embed_devKFP/*.root","embedding SL24y devKFP",1,    ",trigger:Evt.isTrigger(trigList_2019_19AuAu),TOF:Evt.nBTOFMatch>100"}, //emb test
    //{"/home/petrchal/tpcAna/ntup/2019_19GeV_SL24embed_devKFP/*.root","embedding SL24y devKFP devKFP nK3piP>0",1,    "nK3piP:Evt.nK3piP>0,trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //emb test
    //{"/home/petrchal/tpcAna/ntup/2019_19GeV_SL24embed_devKFP/*.root","embedding SL24y devKFP nK3piP==1",1,          "nK3piP:Evt.nK3piP==1,trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //emb test
    //obsolete
    {"/home/petrchal/tpcAna/ntup/2019_19GeV_SL24embed/*.root","embedding SL24y",1,"trigger:Evt.isTrigger(trigList_2019_19AuAu)"}, //emb test

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
  const int order[]={3,4,6};//{5,3,6,4};//{5,6};//{5,3,6,4}; //{5,6,4};
  //const int order[]={3,5,1};//{3,1,6,2};//{1,3};//{3,5};//{1,9};//{2,0,4};//{2,4,3,3};
  const Long64_t nEntriesLimit=TTree::kMaxEntries;//# or your number of entrie







