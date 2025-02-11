#include "K3piLib.C"

//this is mean as global setup to make sure that I use alway the same files and cut
//one can override in the script
 

 /*  just an example from collider data
 - inludes embedding
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
*/

  /* TFileDescription data_2020_FXT[]={
    {"../../ntup/2020_FXT_31p2_SL23d_noLcut/*.root","31.2GeV FXT SL23d",0,NULL}, 
    {"../../ntup/2020_FXT_13p5_SL23e_noLcut/*.root","13.5GeV FXT SL23e",0,NULL}, 
    {"../../ntup/2020_FXT_5p75_SL23e_noLcut/*.root","5.75GeV FXT SL23e",0,NULL}, 
    {"../../ntup/2020_FXTTFG_5p75_withMap/*.root","5.75GeV FXT TFG topoMap",0,NULL}, 
      }; */

 TFileDescription data_2019_FXT_4p59[]={
    {"../../ntup/2019FXT_4p59_SL23d_picoDst/*.root","SL23d_pico",0,NULL}, 
    {"../../ntup/2019FXT_4p59_TFG24c_picoDst/*.root","TFG24c picoDst",0,NULL}, 
    {"../../ntup/2019FXT_4p59_TFG24c_MuDst/*.root","TFG24c MuDst",0,NULL}, 
    {"../../ntup/2019FXT_4p59_TFG24d_picoDst//*.root","TFG24d picoDst",0,NULL}, 
    {"../../ntup/2019FXT_4p59_TFG24d_MuDst/*.root","TFG24d MuDst",0,NULL}, 
};



//=======================
//This select what data will be used globally    

  //TFileDescription* files=data19GeV;
  //TFileDescription* files=data_2020_FXT;
  TFileDescription* files=data_2019_FXT_4p59;
 
  const int nFiles=5; //for efficiency plots only nFiles<=2 possible
  const int order[]={0,1,2,3,4};//{0,1,2};
  const Long64_t nEntriefsLimit=TTree::kMaxEntries;//# or your number of entrie





