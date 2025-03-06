#ifndef CUTS_INCLUDE
#define CUTS_INCLUDE


#include "TString.h"
#include <iostream>
#include <string>
#include <map>
#include "tpcPadPlanes.C"
#include "K3piLib.C"


void AddNewVar(const char* var, const char *val) {NewVars[var]=val;}


const Float_t K_m=0.494;
//int particlesPDG[nPIDS] = {100321, 200321,-100321, -200321}; //k->3pi only, K->3pi with found K

//----------------- cut variations ----------------------------------------
TCutVariation vary_3piVtx_chi2ndf={"cut_mom_ch2ndf",{"cut_mom_ch2ndf+0.2","cut_mom_ch2ndf+0.4"}};
TCutVariation vary_lastPointDiff={"cut_lastPointDiff",{"cut_lastPointDiff+3.","cut_lastPointDiff-3."}};
TCutVariation vary_EvtVz={"cut_Vz",{"cut_Vz+15.","cut_Vz-15."}};
TCutVariation vary_dpDecay={"cut_dpDecay",{"cut_dpDecay+0.05","cut_dpDecay-0.02"}};
TCutVariation vary_Minv={"cut_mom_Minv",{"cut_mom_Minv+0.02", "cut_mom_Minv-0.05"}};//0.015 
TCutVariation vary_daughter_Nhits={"cut_daugh_nhits",{"12.", "13."}};//11
TCutVariation vary_DCAxy={"PvtxDcaXY_corrected",{"float(PvtxDcaXY_corrected+0.1)"}};//11


//------ trigger lists   ----------------------------------------
std::vector<unsigned int> trigList_2020_FXT_5p75AuAu;
std::vector<unsigned int> trigList_2019_FXT_4p59AuAu;

void InitTriggerLists(){
  // no triggers for FXT sofar!!!

 //2020_FXT_5p75AuAu
 //720000 (epde-or-bbce-or-vpde-tof1) [It is probably this one because that has a bigger data set]
 //720007 (epde-or-bbce-or-vpde-tof1-etof) //not used so far - small statistics
 trigList_2020_FXT_5p75AuAu= std::vector<unsigned int> {720000}; 
 trigList_2019_FXT_4p59AuAu= std::vector<unsigned int> {}; //no triggers so far
  
 std::sort(trigList_2020_FXT_5p75AuAu.begin(), trigList_2020_FXT_5p75AuAu.end());
 std::sort(trigList_2019_FXT_4p59AuAu.begin(), trigList_2019_FXT_4p59AuAu.end());
}

//------EVENT CUTS here----------------------------------------
//------ this a to ensure that all theprocedures use the same event cut
/* example -------------
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
*/

//-------------
K3PiCut EventCut_2020_FXT(){ //so far empty cut

  K3PiCut Event_cut;
  Event_cut["trigger"]="1"; //pass all

  Event_cut["trigger"]="Evt.isTrigger(trigList_2020_FXT_5p75AuAu)"; 
  Event_cut["Vz"]="(Evt.Vz>150)&&(Evt.Vz<250)"; //embedding width
 // Event_cut["VPDdif"]="fabs(Evt.vzVpd-Evt.Vz)<5"; //VPD cut ..not usable for FXT


  return Event_cut;
}  


//-------------
K3PiCut EventCut_2019_FXT(){ //so far empty cut

  K3PiCut Event_cut;
  Event_cut["trigger"]="1"; //pass all

  //no triggers sofar
  Event_cut["Vz"]="(Evt.Vz>150)&&(Evt.Vz<250)"; //embedding width
 // Event_cut["VPDdif"]="fabs(Evt.vzVpd-Evt.Vz)<5"; //VPD cut ..not usable for FXT


  return Event_cut;
}  
//!!!!!!!assign which cut is globaly used!!!!
//std::function<K3PiCut()> K3piCut_EventCut =EventCut_2021_7p7AuAu;
//std::function<K3PiCut()> K3piCut_EventCut =EventCut_2020_FXT;
std::function<K3PiCut()> K3piCut_EventCut =EventCut_2019_FXT;

//==========END of Event cuts =============================



// =============== 3pi vertex reconstruction quality =====
//Cuts used for reconstruction of the 3pi vertex

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

     cout<<"Petr -CHP1"<<endl;

     AddNewVar("cut_3pi_ch2ndf","1.");
     VertexQA_cut["3piVtx_chi"]="(mother_chi2ndf<cut_3pi_ch2ndf)"; //30 -cut off in extraction, the DNF shoudl be 5?
     //3piVtxQA_cut["3piVtx_chi"]="(mother_chi2ndf<1.5)&&(mother_chi2ndf>0.8)"; //30 -cut off in extraction, the DNF shoudl be 5?
     
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
     //cleans of 3pi signal - more sensitive in data -
     //improves signal/bg but also removes statistics - needs to be tested for effects on final results
     tmp["daughter_DCA_PV"]="(d.PvtxDca_official[0]>20) && (d.PvtxDca_official[1]>20)&& (d.PvtxDca_official[2]>20)";
     tmp["daughter_chi2"]="(d.match_chi2[0]>500) && (d.match_chi2[1]>500) && (d.match_chi2[2]>500) ";
     
     //tmp["daughter_lastHit"]="(d.lastPointR[0]>155) && (d.lastPointR[1]>155)&& (d.lastPointR[2]>155)";
     //tmp["nhits_posrat"]="( (NhitsOK(d.nhits[0],decay_Vr)) &&  (NhitsOK(d.nhits[1],decay_Vr)) &&  (NhitsOK(d.nhits[2],decay_Vr)) )"; not working for iTPC
     
     //skip central membrane
     //tmp["Z_decay"]="(fabs(decay_Vz)<130)&&(fabs(decay_Vz)>30)";
     
     //VertexQA_cut["3piVtxAdds"]=tmp.Str();
    


     return VertexQA_cut;
}
/// ===END of  3pi vertex reconstruction quality =====


//================= Kinematics cuts========================

//--------------FXT---------------------
//IMPORTANT: hence all kinematics cuts are done at after back propagation to primary vertex!!!
K3PiCut Setup_FXT_3piVtxKinematics(){

  K3PiCut kin_cut;

   //base cut
    kin_cut["pt"]="(mother_pt_PVX>0.2)&&(mother_pt_PVX<0.9)"; 
    kin_cut["eta"]="(mother_eta_PVX>-2.)&&(mother_eta_PVX<0.)";
  
   //should always be on  - nothing is macthed below 80 
   //is this true for FXT?
   kin_cut["decay_Vr"]="(decay_Vr>80)";

   //not sure what the Vr cut should be for FXT
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
   //res["kaon_DCA"]="PvtxDca_corrected<1";
   //res["kaon_nhits"]=" d.nhits[K_match]>20";
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
    //KaonMatching_cut["K_dp"]="(d.dp_Decay[K_match]<cut_dpDecay)";  //standard cut deduceed from pure simulations is 200MeV


    return KaonMatching_cut;
}

//------------------------------------
//3piVtx selection without PID (both charges, matched and unmatched)
   K3PiCut K3piCut_VertexCanditate(){
   K3PiCut res= Setup_3piVertexQA()+Setup_FXT_3piVtxKinematics();
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
  InitPadRadii_iTPC();


  //for DCA adjustments - move and smear embedding
  //AddNewVar("PvtxDcaXY_corrected","float((d.qaTruth[K_match]>0.0)?((d.PvtxDcaXY_official[K_match]*1.1+0.0456)):d.PvtxDcaXY_official[K_match])");
  //AddNewVar("PvtxDca_corrected","float(sqrt(PvtxDcaXY_corrected*PvtxDcaXY_corrected+d.PvtxDcaZ_official[K_match]*d.PvtxDcaZ_official[K_match]))");


  //print all cuts for information: also to define variable
  cout<<endl<<"EventCut:"<<endl;
  cout<<"  "<<K3piCut_EventCut().Str()<<endl;

  cout<<"3piVtxKinematics:"<<endl;
  cout<<"  "<<Setup_FXT_3piVtxKinematics().Str()<<endl;
  cout<<"Setup_3piVtxQA:"<<endl;
  //cout<<"  "<<Setup_3piVertexQA().Str()<<endl;
  cout<<"Setup_KaonMatching:"<<endl;
  cout<<"  "<<Setup_KaonMatching().Str()<<endl;
  cout<<"Setup_MCvertex:"<<endl;
  cout<<"  "<<Setup_MCvertex().Str()<<endl;
  
  cout<<endl<<endl;
 
  
}

#endif