                                                             
#ifndef STAR_K3pi
#define STAR_K3pi
//#define __DEVT__


#include "TObject.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "KFParticle/KFParticle.h"
#include "StarClassLibrary/StThreeVectorD.hh"

class StKFParticleInterface;
class StKFParticlePerformanceInterface;
// class KFParticle;
class StPicoDst;
class StMuDst;
class TNtuple;
class TFile;
class TChain;
class TTree;
class StRefMultCorr;

class TDaughter : public TObject {
   public:
    TDaughter(){Clear();}
    void Clear();
    //TODO pridat an naboj a pt projekci...
    //track properties
    Int_t 
      id=0,index=0,charge=0, 
      nhits=0, nhits_dEdx=0,nhits_pos=0;
    Float_t dEdx=0,lastPointR=0;    //lastPointR - radial position of last hit from muDst or from picoDst hitmap
    
    
    //reconstructed (at beam DCA) as save in mu(pico)DST .. i.e. global momentum
    //note:th p,pt,eta are saved for convenience only..could save some space
    Float_t p=0,pt=0,eta=0,phi=0, px=0,py=0,pz=0;
    //recalculated at decay point
    Float_t decay_p=0,decay_pt=0,decay_eta=0,decay_phi=0, decay_px=0,decay_py=0,decay_pz=0,
     phi_wrt_Vr; //decay angle in respect to decay position vector
     
     //DCA and matching
     //isBest this is useful only for the parent track, otherwise -1
     //
    Bool_t  isBest; 
    Float_t 
      DecayDca_KF=0,DecayDca_mu=0,     //DCA to 3pi   vertex from KF and from muDst(picoDst)
      PvtxDca_KF=0,PvtxDca_official=0, 
      PvtxDca_mu=0,
      dp_Decay=0,dp_decay_KF=0,dp_PVX=0,
     //information from helix 
     helix_R=-1, helix_Cr=-1, helix_lowR=0,helix_hiR=0, //radius and distance of centre from beam in transverse plane
     //from KFParticle
     decay_dl=0,
     //others
     pdg=0,idTruth=-5,qaTruth=-5;

   ClassDef(TDaughter,1) 
   };

//basically copied from picoDst class
class TEvInfo: public TObject{
   public:
     TEvInfo(){Clear();}
     void Clear();
     std::vector<unsigned int> getTriggerIds() const { return triggerIds; }
     bool isTrigger(unsigned int) const;
     bool isTrigger(std::vector<unsigned int> &trigs);  
     //TODO setters ..trigers, from event

     //to make sure that the righ PV is used  the values of PV position are filled from KFP 
     bool Fill(StMuDst *dst, KFParticle &primVtx);
     bool Fill(StPicoDst *dst, KFParticle &primVtx);

     StThreeVectorD primVtx_StVec(){return StThreeVectorD(Vx,Vy,Vz);}
     TVector3 primVtx_TVec(){return TVector3(Vx,Vy,Vz);}

    // event properties   
     Int_t runId,eventId;
   
     Float_t Vx,Vy,Vz,vzVpd,
     ZDCx,BBCx; // coincidence rates  
     //same as in picoDST
     int refMult; //via tracks (-0.5<eta<0.5)
     int gRefMult;//global tracks in |eta|<0.5
  private:
     std::vector<unsigned int> triggerIds;
     ClassDef(TEvInfo,1) 
};

class TK3pi : public TObject {
   public:
    TK3pi();
    void Clear();
    TEvInfo &EvInfo(){return Evt;}
    TDaughter& daughter(int i){return *((TDaughter*)(d[i]));}
    TDaughter* pdaughter(int i){return (TDaughter*)d[i];}
       
    Int_t mother_PID;
    Char_t mother_isMc;
    
    Int_t matchedKF=0;
    Int_t matchedGeom=0;

     // decay position 
     Float_t decay_Vr, decay_Vx, decay_Vy, decay_Vz,
    //momentum at the decay vertex from KFP
     mother_pt, mother_px,     mother_py,     mother_pz,    mother_eta,      mother_phi, 
     //recalculated after transport to PVTX
     mother_pt_PVX, mother_px_PVX, mother_py_PVX, mother_pz_PVX, mother_eta_PVX, mother_phi_PVX,
     mother_m,  mother_PV_chi2, mother_PV_l, mother_PV_dl;

    //first three are decay product
    //[3] - matched by KFP
    //[4] - matched by geometrical cuts
    TClonesArray d;
    TEvInfo Evt;
    //TDaughter d[5]; 

    ClassDef(TK3pi,1) 
   };


#endif

