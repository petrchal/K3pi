#ifndef CUTS_INCLUDE
#define CUTS_INCLUDE


#include "TString.h"
#include <iostream>
#include <string>
#include <map>
#include "tpcPadPlanes.C"
#include "K3piLib.C"


void AddNewVar(const char* var, const char *val) {
    cout<<"Adding variable: "<<var<<"="<<val<<endl;
    NewVars[var]=val;
}



const Float_t K_m=0.494;
//int particlesPDG[nPIDS] = {100321, 200321,-100321, -200321}; //k->3pi only, K->3pi with found K

//----------------- cut variations ----------------------------------------
TCutVariation vary_3piVtx_chi2ndf={"cut_3pi_ch2ndf",{"cut_3pi_ch2ndf+0.2","cut_3pi_ch2ndf+0.4"}};
TCutVariation vary_lastPointDiff={"cut_lastPointDiff",{"cut_lastPointDiff+3.","cut_lastPointDiff-3."}};
TCutVariation vary_EvtVz={"cut_Vz",{"cut_Vz+15.","cut_Vz-15."}};
TCutVariation vary_dpDecay={"cut_dpDecay",{"cut_dpDecay+0.05","cut_dpDecay-0.02"}};
TCutVariation vary_Minv={"cut_mom_Minv",{"cut_mom_Minv+0.02", "cut_mom_Minv-0.05"}};//0.015 
TCutVariation vary_daughter_Nhits={"cut_daugh_nhits",{"12.", "13."}};//11
TCutVariation vary_DCAxy={"PvtxDcaXY_corrected",{"float(PvtxDcaXY_corrected+0.1)"}};//11

//runId list
std::vector<unsigned int>  runIds_2019;

//------ trigger lists   ----------------------------------------
std::vector<unsigned int> trigList_2019_19AuAu;  //2019 19GeV AuAu same as in embedding
std::vector<unsigned int> trigList_2018_27AuAu;  //2018 27GeV AuAu same as in embedding
std::vector<unsigned int> trigList_2017_54AuAu;
std::vector<unsigned int> trigList_2021_7p7AuAu;
/* does not work so far
bool IsGoodTrigger( std::vector<unsigned int> &trigs,  std::vector<unsigned int> &goodList){
    //should be already in data, but...
    std::sort(trigs.begin(), trigs.end());
  
    std::vector<int> v_intersection;
    std::set_intersection(trigs.begin(), trigs.end(),
                          goodList.begin(), goodList.end(),
                          std::back_inserter(v_intersection));

    return (v_intersection.size()>0);
}

bool IsGoodTrigger_2018( std::vector<unsigned int> &trigs){
  return IsGoodTrigger(trigs, trigList_2018);
}
*/

void InitTriggerLists(){
    //vector<int> vect1{ 10, 20, 30 };
 
   //vector<int> vect2(vect1.begin(), vect1.end());
  std::vector<unsigned int> tr27{610001,610011,610021,610031,610041,610051}; //2018 27GeV AuAu same as in embedding
  
  std::vector<unsigned int> tr54_all{580001,580011,580021,580051,580054,580054}; //2017 54GeV AuAu 
  std::vector<unsigned int> tr54_minb{580001,580011,580021}; //2017 54GeV AuAu 
  //these are not mutually exclusive
  //5800[0-2]1 minbias trigger, mainly 580021 - 1.1B evts 
  //580051 - minbias-50-hlt   700M evts
  //580054 - minbias-vpd-hlt  660M evts
  //580055 - minbias-zdc-hlt  640M evst 

 //2019 19GeV AuAu
 std::vector<unsigned int> tr19_minb{640011,640021,640031,640041,640051}; //2019 19GeV AuAu minbias ...there are 1.3B according to web
 std::vector<unsigned int> tr19_hlt{640012,640022,640032}; //2019 19GeV AuAu hltgood ..shoudl be subset of minb

//2021 7.7GeV AuAu 
 std::vector<unsigned int> tr7p7_minb{ 810010, 810020, 810030,810040}; 


  trigList_2018_27AuAu=tr27;
  trigList_2017_54AuAu=tr54_minb;
  trigList_2019_19AuAu=tr19_minb;
  trigList_2021_7p7AuAu=tr7p7_minb;
  //for (unsigned int x : trigList_2018) cout <<"Trig="<< x << " "<<endl;
  
  std::sort(trigList_2018_27AuAu.begin(), trigList_2018_27AuAu.end());
  std::sort(trigList_2017_54AuAu.begin(), trigList_2017_54AuAu.end());
  std::sort(trigList_2019_19AuAu.begin(), trigList_2019_19AuAu.end());
  std::sort(trigList_2021_7p7AuAu.begin(), trigList_2021_7p7AuAu.end());
}

//--------------------runlists
std::vector<unsigned int> runList_2019_19AuAu;  //2019 19GeV AuAu same as in embedding
void InitRunLists(){
  //read runlist from file
  std::ifstream is("/home/petrchal/tpcAna/ntup/2019_19GeV_SL23_embed_final/runIds");
  if (!is.is_open()) {  std::cerr << "Error opening file!" << std::endl;}
  unsigned int num;
  while (is >> num) {
        runList_2019_19AuAu.push_back(num);
    }
  is.close();
 
 std::cout << "Read " << runList_2019_19AuAu.size() << " numbers" << std::endl;
 
 //print the numbers to stdout
  std::cout << "numbers read in:\n";
  std::copy(runList_2019_19AuAu.begin(), runList_2019_19AuAu.end(), 
            std::ostream_iterator<unsigned int>(std::cout, " "));
  std::cout << std::endl;

  std::sort(runList_2019_19AuAu.begin(), runList_2019_19AuAu.end());
}

bool goodRunId(unsigned int id) {
  //cout<<" petr runId="<<id<<endl;
  bool res= std::find(runList_2019_19AuAu.begin(), runList_2019_19AuAu.end(), id) != runList_2019_19AuAu.end();
  //TString s=" petr runId=";s+=id;s+=" found=";s+=res;
  //cout<<s.Data()<<endl;
  return res;
}


//------EVENT CUTS here----------------------------------------
//------ this a to ensure that all theprocedures use the same event cut
K3PiCut EventCut_2019_19AuAu(){

  K3PiCut Event_cut;

  Event_cut["trigger"]="Evt.isTrigger(trigList_2019_19AuAu)";
  
  //2019 has rather wide Vz distribution, I'm loosing about 60%
  AddNewVar("cut_Vz","80."); //50 alows cut variation
 //cumulative loss of 70% for 2019
  Event_cut["Vz"]="(Evt.Vz>-cut_Vz)&&(Evt.Vz<cut_Vz)"; //embedding width
  //Event_cut["Vz"]="(Evt.Vz>20)&&(Evt.Vz<80)"; //right half
  //Event_cut["Vz"]="(Evt.Vz<-20)&&(Evt.Vz>-80)"; //left half
  
  Event_cut["VPDdif"]="fabs(Evt.vzVpd-Evt.Vz)<5"; //VPD cut ..maybe to tight, but ok .. takes off another 30%

  //Event_cut["nK3piP"]="Evt.nK3piP<=1"; //"there is a significant change from SL21 to SL23"
 
  //extracted from emebedding - using exatly same data
  Event_cut["runId"]="goodRunId(Evt.runId)";
 
 /*
  //optimization for 27GeV data
  Event_cut["BBC"]="(Evt.BBCx<500000)";
  Event_cut["ZDC"]="(Evt.ZDCx>800) && (Evt.ZDCx<1800)";
  Event_cut["gRefMult"]="(Evt.gRefMult<320)";
  Event_cut["TOFmatch"]="(Evt.nBTOFMatch<340)";
  */

 //optimization for 19GeV data
 //Event_cut["gRefMult"]="(Evt.gRefMult<420)";
 //Event_cut["TOFmatch"]="(Evt.nBTOFMatch<460)&&(Evt.nBTOFMatch>100)";

  //Event_cut["ZDC"]="(Evt.ZDCx>50) && (Evt.ZDCx<800)";
  //Event_cut["ZDC"]="(Evt.ZDCx>200) && (Evt.ZDCx<600)";
  //Event_cut["eventId"]="(Evt.eventId<1200000)";


  return Event_cut;
} 

//------ this a to ensure that all theprocedures use the same event cut
K3PiCut EventCut_2018_27AuAu(){

  K3PiCut Event_cut;

  Event_cut["trigger"]="Evt.isTrigger(trigList_2018_27AuAu)";
  
   AddNewVar("cut_Vz","50.");

  //Event_cut["Vz"]="(Evt.Vz>-cut_Vz)&&(Evt.Vz<cut_Vz)"; //embedding width
  //Event_cut["Vz"]="(Evt.Vz>30)&&(Evt.Vz<80)"; //right half
  //Event_cut["Vz"]="(Evt.Vz>10)&&(Evt.Vz<50)"; //right half
   Event_cut["Vz"]="(Evt.Vz<-10)&&(Evt.Vz>-50)"; //left half
  
  Event_cut["VPDdif"]="fabs(Evt.vzVpd-Evt.Vz)<5"; //VPD cut ..maybe to tight, but ok .. takes off another 30%

   
  //extracted from emebedding - using exatly same data
  //Event_cut["runId"]="goodRunId(Evt.runId)";
 
 
  //optimization for 27GeV data
  Event_cut["BBC"]="(Evt.BBCx<500000)";
  Event_cut["ZDC"]="(Evt.ZDCx>800) && (Evt.ZDCx<1800)";
  Event_cut["gRefMult"]="(Evt.gRefMult<320)";
  Event_cut["TOFmatch"]="(Evt.nBTOFMatch<340)";
  
  return Event_cut;
} 

//-------------
K3PiCut EventCut_2021_7p7AuAu(){

  K3PiCut Event_cut;

  Event_cut["trigger"]="Evt.isTrigger(trigList_2021_7p7AuAu)";
 
  AddNewVar("cut_Vz","50.");
  Event_cut["Vz"]="(Evt.Vz>-cut_Vz)&&(Evt.Vz<cut_Vz)"; //embedding width
  

  Event_cut["VPDdif"]="fabs(Evt.vzVpd-Evt.Vz)<5"; //VPD cut ..maybe to tight, but ok

  //when comparing to TFG production
  Event_cut["runId"]="(Evt.runId<22042000)";

  //good run && luminosity
  //bad run list should be part of the candidate extraction
 
  return Event_cut;
}

//-------------
K3PiCut EventCut_2020_FXT(){ //so far empty cut

  K3PiCut Event_cut;
  Event_cut["passAll"]="1";

  //Event_cut["trigger"]="Evt.isTrigger(trigList_2021_7p7AuAu)";
 

  return Event_cut;
}    
//!!!!!!!assign which cut is globaly used!!!!
//std::function<K3PiCut()> K3piCut_EventCut =EventCut_2021_7p7AuAu;
std::function<K3PiCut()> K3piCut_EventCut =EventCut_2019_19AuAu;
//std::function<K3PiCut()> K3piCut_EventCut =EventCut_2018_27AuAu;
//std::function<K3PiCut()> K3piCut_EventCut =EventCut_2020_FXT;

//==========END of Event cuts =============================



// =============== 3pi vertex reconstruction quality =====
//Cuts used for reconstruction of the 3pi vertex

//return maximim number of daughter hist based on Vr
//not working for iTPC geometry!!!
//values are from simulations 
/*
float NhitsOK(float nhits, float Vr){
    float y1=-0.2*Vr+57;
    float y2=-0.5*Vr+94.5;
    float maxN=y1;
    if (y1>y2) maxN=y2;

    return maxN;
    if (nhits>maxN) return 0;
 
    if (nhits/maxN < 0.5) return 0;

    return 1;   
}
*/

//get number of hits from radial position for old TPC
//NOTE: Must Call InitPadRadii First
int MaxHits(double r){
   int n=PadFromR(r);
   n++;
   if (n<0) return 0;
   return n;
}

//for old TPC
int MaxHitsDaughter(double r){
  return nPadRows-MaxHits(r);
}

//quality of reconstruction of the 3pi vertex + PID
//can be varied for systematic checks
//on purpose does not contain PID cut
K3PiCut Setup_3piVertexQA(){
  
     K3PiCut VertexQA_cut;

     AddNewVar("MaxHitsDaughter","MaxHitsDaughter(decay_Vr)");

     AddNewVar("cut_3pi_ch2ndf",".2");  //0.1 - strict cut, most of 3pi are below 0.2
     VertexQA_cut["3piVtx_chi"]="(mother_chi2ndf<cut_3pi_ch2ndf)"; //30 -cut off in extraction, the DNF shoudl be 5?
     
     //TODO - possible to include
     //3piVtx_PV_l;
     //3piVtx_PV_dl;
     //3piVtx_PV_chi2;
     
    AddNewVar("cut_mom_Minv","0.015");
    VertexQA_cut["minv"]="(fabs(mother_m-0.494)<cut_mom_Minv)"; //sigma=0.005 (even little less)
          
    AddNewVar("cut_daugh_nhits","11.");
    VertexQA_cut["nhits_daughters"]="(d.nhits[0]>=cut_daugh_nhits && d.nhits[1]>=cut_daugh_nhits && d.nhits[2]>=cut_daugh_nhits)"; 
     
    //this seems to remove all short tracks - BAD seem the nhits_pos is not filled correctly for secondaries
    // 3piVtxQA_cut["nhits_posrat"]="((float)d.nhits[0]/(float)d.nhits_pos[0]>0.51)&&((float)d.nhits[1]/(float)d.nhits_pos[1]>0.51)&&((float)d.nhits[3]/(float)d.nhits_pos[3]>0.51)";
     
 
     //additional cleanup cuts - used for 2019 19GeV
     //noteL 
     K3PiCut tmp;
     //cleans 3pi signal - more sensitive in data
    
     //similar effects - take only secondaries
     tmp["daughter_DCA_PV"]="(d.PvtxDca_official[0]>20) && (d.PvtxDca_official[1]>20)&& (d.PvtxDca_official[2]>20)";
     tmp["daughter_chi2"]="(d.match_chi2[0]>500) && (d.match_chi2[1]>500) && (d.match_chi2[2]>500) ";
     //tmp["daughter_lastHit"]="(d.lastPointR[0]>155) && (d.lastPointR[1]>155)&& (d.lastPointR[2]>155)";
     //tmp["nhits_posrat"]="( (NhitsOK(d.nhits[0],decay_Vr)) &&  (NhitsOK(d.nhits[1],decay_Vr)) &&  (NhitsOK(d.nhits[2],decay_Vr)) )"; not working for iTPC
     
      
     VertexQA_cut["3piVtxAdds"]=tmp.Str();
    


     return VertexQA_cut;
}
/// ===END of  3pi vertex reconstruction quality =====


//================= Kinematics cuts========================
//this sets kinematics of the found 3pi vertex - that select the kinemaics of the mother kaon
//IMPORTANT: hence all kinematics cuts are done at after back propagation to primary vertex!!!
K3PiCut Setup_3piVtxKinematics(){

  K3PiCut kin_cut;

  //base cut
  kin_cut["pt"]="(mother_pt_PVX>0.25)&&(mother_pt_PVX<0.9)"; //bellow 200MeV it is bad ,above 1GeV  also , but not sure why
  
  //kin_cut["pt"]="(mother_pt_PVX>0.25)&&(mother_pt_PVX<0.4)";
  //pt="&&(mother_pt_PVX>0.2)&&(mother_pt_PVX<1.)"; //basic pt cut
    
   //kin_cut["eta"]="(fabs(mother_eta_PVX)<0.8)"; //standard
    kin_cut["eta"]="(fabs(mother_eta_PVX)<1.2)";
   
   // 110cm -base cut for old TPC - nor 3pi bellow 110cm
   //kin_cut["decay_Vr"]="(decay_Vr>110)"; //simulation cutoff for flat pt

   //base cut for iTPC - zero efficiency below 
   //!!!MUST be applied!!!
   kin_cut["decay_Vr"]="(decay_Vr>80)"; 

   // long track without iTPC: 2018
   //kin_cut["decay_Vr"]="(decay_Vr>140)&&(decay_Vr<160)"; //170
 
   // long track with iTPC: from 2019 up
   kin_cut["decay_Vr"]="(decay_Vr>130)&&(decay_Vr<160)"; //160 may be safer, could go to 120
   //kin_cut["decay_Vr"]="(decay_Vr<120)&&(decay_Vr>80)"; //short tracks
 
   // inner/outer divide
   //kin_cut["decay_Vr"]="(decay_Vr>115)&&(decay_Vr<122)"; 
   

   //kin_cut["decay_Vr"]="(decay_Vr>130)&&(decay_Vr<140)";
   //kin_cut["decay_Vr"]="(decay_Vr>150)"; //160 may be safer, could go to 120
  
   //short tracks
   //kin_cut["decay_Vr"]="(decay_Vr<110)"; //160 may be safer, could go to 120
   
   //--Vz---
   //skip central membrane
   //tmp["decay_Z"]="(fabs(decay_Vz)<130)&&(fabs(decay_Vz)>30)";
   // kin_cut["decay_Z"]="decay_Vz>10";
   //this does not work ...not sure why ..mixing two branches?
    //kin_cut["sameSide"]=" (Evt.Vz<-10 && decay_Vz<-10)|| (Evt.Vz>10 && decay_Vz>10)";

   //kin_cut["decay_Vz"]="((decay_Vz<-5)||(decay_Vz>5))";

     //kin_cut["membrane"]="((decay_Vz<-0 && Evt.Vz>0)||(decay_Vz>0 && Evt.Vz<0))";
  //kin_cut["nonmembrane"]=="((decay_Vz<0 && Evt.Vz<0)||(decay_Vz>0 && Evt.Vz>0))";
  return kin_cut;
}

//--------------FXT---------------------
//IMPORTANT: hence all kinematics cuts are done at after back propagation to primary vertex!!!
K3PiCut Setup_FXT_3piVtxKinematics(){

  K3PiCut kin_cut;

  //base cut
  kin_cut["pt"]="(mother_pt_PVX>0.2)&&(mother_pt_PVX<0.9)"; 
  kin_cut["eta"]="mother_eta_PVX>-2.)&&(mother_eta_PVX)<0.)";
  
  //nor sure what the Vr cut should be for FXT
   // long track without iTPC: 2018
   //kin_cut["decay_Vr"]="(decay_Vr>140)&&(decay_Vr<165)"; //170
 
   // long track with iTPC: from 2019 up
   //kin_cut["decay_Vr"]="(decay_Vr>130)&&(decay_Vr<160)"; //160 may be safer, could go to 120
   //kin_cut["decay_Vr"]="(decay_Vr>120)&&(decay_Vr<160)";
   //kin_cut["decay_Vr"]="(decay_Vr>150)"; //160 may be safer, could go to 120
  
   
  return kin_cut;
}


//-----------------------------
//this is THE analysis like cut on track propertie
//Watch out - this must NOT include the kinematics cuts (pt,eta, phi)
K3PiCut K3piCut_KaonTrackCut(){
   K3PiCut res;

   //ANALYSIS CUTS
   //res["kaon_DCA"]="d.PvtxDca_official[K_match]<1";
   res["kaon_DCA"]="PvtxDca_corrected<2";
   res["kaon_nhits"]=" d.nhits[K_match]>20";
   //res["kaon_hits_ratio"]="(d.nhits[K_match]/d.nhits_pos[K_match])>0.5"; 

   //note: nhits/npos .... not tested yet the npos is not calculated correctly for track not reaching outer edge of TPC
   //res["nhits_dEdx"]="";

   return res;
}

//-----------------------------------------------------
// selecting 3pi vertex coming from MC
K3PiCut Setup_MCvertex(){
   K3PiCut res;
   //standard
   //should remove original (non-embeded 3pi vertexes) 
   // at least one track is MC - allows mixing to backround
   res["qaTruth"]="((d.qaTruth[0]+d.qaTruth[1]+d.qaTruth[2])>0)";
 
   //this is basicaly equivalent
   //when is MC is on, one has to compare to clean MC tracks .. .(d.qaTruth[K_match]>50)
   res["isMC"]="(mother_isMc==1)";
   //res["qaTruth"]="(d.qaTruth[0]>50)&&(d.qaTruth[1]>50)&&(d.qaTruth[2]>50)";
  
    return res;
}

//-----------------------------------------------------
// describing matching of primary kaon to 3pi vertex
// this is matching via helix - dp
K3PiCut Setup_KaonMatching(){
 
    K3PiCut KaonMatching_cut;
    //matchedGeom is global variable which select matching by dp
    KaonMatching_cut["K_match"]="( matchedGeom && (d.nhits[K_match]>10)) ";//&& d.nhits[K_match]>20";

     //DCA and matching
     //isBest this is useful only for the parent track, otherwise -1
     // -1 ..was not touched
     //0 .. was touched but, eventualy better candidate was found 
     // 1 .. best matched track by geometry
     // 2 .. if matching was done via dca at decay vertex, not dp, this would not be the best candidate
    KaonMatching_cut["K_best"]="(d.isBest[K_match]>=0)"; //this should not be necessary..but is .there is some bug...there are isBest=-1 in the data
    //should be Eequivalent to selecting matchedGeom==1
  
    AddNewVar("cut_lastPointDiff","5.");//5
    KaonMatching_cut["K_lastR"]="((d.lastPointR[K_match]-decay_Vr)<cut_lastPointDiff)";// previously 15,10

    AddNewVar("cut_dpDecay","0.1");
    KaonMatching_cut["K_dp"]="(d.dp_Decay[K_match]<cut_dpDecay)";  //standard cut deduceed from pure simulations is 200MeV


    return KaonMatching_cut;
}

//------------------------------------
//3piVtx selection without PID (both charges, matched and unmatched)
   K3PiCut K3piCut_VertexCanditate(){
   K3PiCut res= Setup_3piVertexQA()+Setup_3piVtxKinematics();
   return res;
}


//---------------------------------------
K3PiCut K3piCut_3piVtx_Kplus(){
   K3PiCut res= K3piCut_VertexCanditate();
   res["PID"]="(mother_PID==100321)";
   return res;
}

//-----------------------------
K3PiCut K3piCut_3piVtx_Kminus(){
   K3PiCut res= K3piCut_VertexCanditate();
   res["PID"]="(mother_PID==-100321)";
   return res;
}


//-----------------------------
K3PiCut K3piCut_Matched_Kplus(){
   K3PiCut res=K3piCut_3piVtx_Kplus()+Setup_KaonMatching()+K3piCut_KaonTrackCut();
   return res;
}


void InitCuts(){

  gROOT->ProcessLine(".L K3pi.cxx+"); //this allows to use isStrigger

  NewVars.Clear();
  InitTriggerLists(); //done once per job
  InitRunLists();
  InitPadRadii_iTPC();


  //for DCA adjustments - move and smear embedding
  //this should be done in a smarter way - separately for each data set
  AddNewVar("PvtxDcaXY_corrected","d.PvtxDcaXY_official[K_match]");
  AddNewVar("PvtxDca_corrected","d.PvtxDca_official[K_match]");
  
  //AddNewVar("PvtxDcaXY_corrected","float((d.qaTruth[K_match]>0.0)?(((d.PvtxDcaXY_official[K_match]-0.04965)*1.1389+0.1458)):d.PvtxDcaXY_official[K_match])");
  //AddNewVar("PvtxDca_corrected","float(sqrt(PvtxDcaXY_corrected*PvtxDcaXY_corrected+d.PvtxDcaZ_official[K_match]*d.PvtxDcaZ_official[K_match]))");


  //print all cuts for information: also to define variable
  cout<<endl<<"EventCut:"<<endl;
  cout<<"  "<<K3piCut_EventCut().Str()<<endl;

  cout<<"3piVtxKinematics:"<<endl;
  cout<<"  "<<Setup_3piVtxKinematics().Str()<<endl;
  cout<<"Setup_3piVtxQA:"<<endl;
  //cout<<"  "<<Setup_3piVertexQA().Str()<<endl;
  cout<<"Setup_KaonMatching:"<<endl;
  //cout<<"  "<<Setup_KaonMatching().Str()<<endl;
  cout<<"Setup_MCvertex:"<<endl;
  //cout<<"  "<<Setup_MCvertex().Str()<<endl;
  
  cout<<endl<<endl;
 
  
}

#endif