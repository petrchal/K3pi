#include "K3pi.h"

#include <algorithm>
#include <list>
#include <iterator>

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

void TK3pi::Clear(){
     runId=-1;eventId=-1;
     Vx=0;Vy=0;Vz=0;
     mother_PID=-1; mother_isMc=-1;
     // decay position 
     decay_Vr=0; decay_Vx=0; decay_Vy=0; decay_Vz=0;
    //momentum at the decay vertex from KFP
     mother_pt=0;     mother_px=0;     mother_py=0;     mother_pz=0;    mother_eta=0;      
     mother_phi=0; 
     //recalculated at PVTX
     mother_pt_PVX=0; mother_px_PVX=0; mother_py_PVX=0; mother_pz_PVX=0; mother_eta_PVX=0; mother_phi_PVX=0;
     mother_m=0;  mother_PV_l=0; mother_PV_dl=0;
   for (int i=0;i<5;i++)daughter(i).Clear();
}

void TEvInfo::Clear(){
    if( !triggerIds.empty() ) {
         triggerIds.clear();
     }
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
