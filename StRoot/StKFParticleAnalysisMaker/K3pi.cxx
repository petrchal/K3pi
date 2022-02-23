#include "K3pi.h"

#include <algorithm>
#include <list>
#include <iterator>

//--- Mu classes ---
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#//--- pico classes ---
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"

ClassImp(TDaughter)
ClassImp(TEvInfo)
ClassImp(TK3pi)


void TDaughter::Clear(){
    id=-1,index=-1,charge=0,
    nhits=-1, nhits_dEdx=-1,nhits_pos=-1,dEdx=0,lastPointR=-1,
    p=0,pt=0,eta=0,phi=0, px=0,py=0,pz=0,
    decay_p=0,decay_pt=0,decay_eta=0,decay_phi=0, decay_px=0,decay_py=0,decay_pz=0, phi_wrt_Vr=-2;
    DecayDca_KF=10000,DecayDca_mu=10000,PvtxDca_KF=10000,PvtxDca_official=10000,
    PvtxDca_mu=10000,isBest=-1,dp_Decay=-10000,dp_decay_KF=-10000,dp_PVX=-10000,
    helix_R=-1, helix_Cr=-1, helix_lowR=0,helix_hiR=0,
    decay_dl=-10000,
    pdg=0,idTruth=-1,qaTruth=-1;
}

//======================================
TK3pi::TK3pi():TObject(),Evt(),d("TDaughter", 5){
      for (int i=0;i<5;i++){
       new (d[i]) TDaughter;
       daughter(i).Clear();
     }
    }

void TK3pi::Clear(){
     mother_PID=-1; mother_isMc=-1;
     // decay position 
     decay_Vr=0; decay_Vx=0; decay_Vy=0; decay_Vz=0;
    //momentum at the decay vertex from KFP
     mother_pt=0;     mother_px=0;     mother_py=0;     mother_pz=0;    mother_eta=0;      
     mother_phi=0; 
     //recalculated at PVTX
     mother_pt_PVX=0; mother_px_PVX=0; mother_py_PVX=0; mother_pz_PVX=0; mother_eta_PVX=0; mother_phi_PVX=0;
     mother_m=-1;  mother_PV_chi2=-1, mother_PV_l=-1; mother_PV_dl=-1;
   for (int i=0;i<5;i++)daughter(i).Clear();
}


 //=================================================


void TEvInfo::Clear(){
    if( !triggerIds.empty() ) {
         triggerIds.clear();
     }
     runId=-1;eventId=-1;
     Vx=0;Vy=0;Vz=0;vzVpd=0;
     ZDCx=-1,BBCx=-1;

     int refMult=-100; 
     int gRefMult=-100;
}


bool TEvInfo::isTrigger(unsigned int id) const {
  return std::find(triggerIds.begin(), triggerIds.end(), id) != triggerIds.end();
}

bool TEvInfo::isTrigger(std::vector<unsigned int> &trigs){
    std::sort(trigs.begin(), trigs.end());
    std::sort(triggerIds.begin(), triggerIds.end());
 
    std::vector<int> v_intersection;
 
    std::set_intersection(trigs.begin(), trigs.end(),
                          triggerIds.begin(), triggerIds.end(),
                          std::back_inserter(v_intersection));

    return (v_intersection.size()>0);
}

bool TEvInfo::Fill(StMuDst *dst, KFParticle &primVtx){
    Clear();
    runId   = dst->event()->runId();
    eventId = dst->event()->eventId(); 
    Vx=primVtx.GetX();
    Vy=primVtx.GetY();
    Vz=primVtx.GetZ();
    assert(0);
    /*
    Int_t runId,eventId;
   
     Float_t Vx,Vy,Vz,vzVpd,
     ZDCx,BBCx; // coincidence rates  
     //same as in picoDST
     int refMult; //via tracks (-0.5<eta<0.5)
     int gRefMult;//global tracks in |eta|<0.5
  private:
     std::vector<unsigned int> triggerIds;
     */
}

bool TEvInfo::Fill(StPicoDst *dst, KFParticle &primVtx){
    Clear();
    Vx=primVtx.GetX();
    Vy=primVtx.GetY();
    Vz=primVtx.GetZ();
    runId   = dst->event()->runId();
    eventId = dst->event()->eventId();
    assert(0);
}
