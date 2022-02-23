//*-- Author : Yuri Fisyak 02/02/2016
#include "StKFParticleAnalysisMaker.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TSystem.h"
//--- KF particle classes ---
#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "K3pi.h"
#include "KFParticleTopoReconstructor.h"
#include "KFPartEfficiencies.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"
//--- Pico classes ---
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
//--- Mu classes ---
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"

//hit topology and position
#include "StEvent/StTrackTopologyMap.h"
#include "StDetectorDbMaker/St_tpcPadConfigC.h"
#include "StDetectorDbMaker/St_tpcPadPlanesC.h"
#include "StDetectorDbMaker/St_itpcPadPlanesC.h"


ClassImp(StKFParticleAnalysisMaker);

//________________________________________________________________________________
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name) : StMaker(name),  fIsPicoAnalysis(true), 
  fProcessSignal(kAllTracks), fCollectTrackHistograms(false), fCollectPIDHistograms(false),
  fDecays(0), fIsProduce3DEfficiencyFile(false), f3DEfficiencyFile(""), 
  fStoreCandidates(false), fPartcileCandidate(), fIsStoreCandidate(KFPartEfficiencies::nParticles, false), fCandidateFile(nullptr), fCandidatesTree(nullptr),
  fKaonAnalysis(false), fKaonFileName("kaons.root"){
  memset(mBeg,0,mEnd-mBeg+1);
  
}
//________________________________________________________________________________
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker() 
{
  SafeDelete(fStKFParticleInterface);
  SafeDelete(fStKFParticlePerformanceInterface);
}

//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Init()
{
  TFile *f = GetTFile();
  if(f) 
  {
    f->cd();
    BookVertexPlots();
    if(fCollectTrackHistograms)
      fStKFParticleInterface->CollectTrackHistograms();
    if(fCollectPIDHistograms)
      fStKFParticleInterface->CollectPIDHistograms();
  }

    
  if(fStoreCandidates)
  {
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    
    fCandidateFile = new TFile("candidates.root", "RECREATE");
    fCandidatesTree = new TTree("Candidates", "Candidates");
    
    fCandidatesTree->Branch("Candidates", &fPartcileCandidate, 32000, 0);

    gFile = curFile;
    gDirectory = curDirectory;
  }

  
  

  if(fKaonAnalysis)
  {
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    

    fKaonFile = new TFile(fKaonFileName, "RECREATE");
     
    fKaonTree = new TTree("kaons","tree of K to 3 pi");
    fKaonTree->Branch("K","TK3pi",&fK,32000,2);
      
    gFile = curFile;
    gDirectory = curDirectory;
  }
  return kStOK;
}
//________________________________________________________________________________
Int_t StKFParticleAnalysisMaker::InitRun(Int_t runumber) 
{
//   assert(StPicoDstMaker::instance());
//   if (StPicoDstMaker::instance()->IOMode() == StPicoDstMaker::ioRead) {
    //TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO Ask Yuri
//     StPicoDstMaker::instance()->SetStatus("*":0);
//     const Char_t *ActiveBranches[] = {
//       "MuEvent"
//       ,"PrimaryVertices"
//       ,"PrimaryTracks"
//       ,"GlobalTracks"
//       ,"StStMuMcVertex"
//       ,"StStMuMcTrack"
//       ,"CovPrimTrack"
//       ,"CovGlobTrack"
//       ,"StStMuMcVertex"
//       ,"StStMuMcTrack"
//       ,"KFTracks"
//       ,"KFVertices"
//       ,"StBTofHit"
//       ,"StBTofHeader"
//     }; 
//     Int_t Nb = sizeof(ActiveBranches)/sizeof(Char_t *);
//     for (Int_t i = 0; i < Nb; i++) StPicoDstMaker::instance()->SetStatus(ActiveBranches[i],1); // Set Active braches
//   }
  return StMaker::InitRun(runumber);
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::PrintMem(const Char_t *opt)
{
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  cout << opt 
       << "\tMemory : Total = " << info.fMemTotal 
       << "\tUsed = " << info.fMemUsed
       << "\tFree = " << info.fMemFree
       << "\tSwap Total = " << info.fSwapTotal
       << "\tUsed = " << info.fSwapUsed
       << "\tFree = " << info.fSwapFree << endl;
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::BookVertexPlots()
{
  TDirectory *dirs[2] = {0};
  dirs[0] = TDirectory::CurrentDirectory(); assert(dirs[0]);
  dirs[0]->cd();
  if (! dirs[0]->GetDirectory("Particles")) {
    dirs[0]->mkdir("Particles");
  }
  dirs[1] = dirs[0]->GetDirectory("Particles"); assert(dirs[1]);
  dirs[1]->cd();
  PrintMem(dirs[1]->GetPath());
  
  fStKFParticleInterface = new StKFParticleInterface;
  for(unsigned int iDecay=0; iDecay<fDecays.size(); iDecay++)
    fStKFParticleInterface->AddDecayToReconstructionList( fDecays[iDecay] );
  bool storeMCHistograms = false;
  if(!fIsPicoAnalysis && fProcessSignal==kMcTracksOnly) storeMCHistograms = true;
  fStKFParticlePerformanceInterface = new StKFParticlePerformanceInterface(fStKFParticleInterface->GetTopoReconstructor(), storeMCHistograms, fIsProduce3DEfficiencyFile);
  if(!f3DEfficiencyFile.IsNull()) {
    fStKFParticlePerformanceInterface->Set3DEfficiency(f3DEfficiencyFile);
  }
  dirs[0]->cd();
  PrintMem(dirs[1]->GetPath());
}
//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Make()
{  
  if(fIsPicoAnalysis)
  {
    fPicoDst = StPicoDst::instance();
    if(!fPicoDst) return kStOK;
  }
  else
  {  
    fMuDst = StMuDst::instance();
    if(!fMuDst) return kStOK;
    else { if(StMuDst::instance()->numberOfPrimaryVertices() == 0 ) return kStOK; }
  }
  
  //find max global track index
  int maxGBTrackIndex = -1;
  if(fIsPicoAnalysis)
  {
    for(unsigned int iTrack = 0; iTrack < fPicoDst->numberOfTracks(); iTrack++){
      StPicoTrack *gTrack = fPicoDst->track(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  else //muDst 
  {
    for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfGlobalTracks(); iTrack++){
      StMuTrack *gTrack = fMuDst->globalTracks(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  vector<KFMCTrack> mcTracks(0);
  vector<int> mcIndices(maxGBTrackIndex+1);
  for(unsigned int iIndex=0; iIndex<mcIndices.size(); iIndex++)
    mcIndices[iIndex] = -1;
  
  
  vector<int> triggeredTracks;
  bool isGoodEvent = false;
  
  //Process the event - fill event and call KFP
  if(maxGBTrackIndex > 0)
    fStKFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex+1);
  if(fIsPicoAnalysis)
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fPicoDst, triggeredTracks);
  else
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fMuDst, mcTracks, mcIndices, fProcessSignal);

//   bool openCharmTrigger = false;
//   if(isGoodEvent) openCharmTrigger =  fStKFParticleInterface->OpenCharmTrigger();
//   fStKFParticleInterface->OpenCharmTriggerCompression(triggeredTracks.size(), fPicoDst->numberOfTracks(), openCharmTrigger);
  //collect histograms
  

  if(isGoodEvent) return kStOk;
  
 
#if 1
    //clean clusters
    int nTracks = 0;

    std::vector<bool> isValidTrack(maxGBTrackIndex+1);
    for(int iTrack=0; iTrack<maxGBTrackIndex+1; iTrack++) isValidTrack[iTrack] = true;

    for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++) {
      const KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];

      if(particle.NDaughters() == 1)
        nTracks = iParticle + 1;

      if( (abs(particle.GetPDG()) > 3001) && (abs(particle.GetPDG()) <= 3029) ) {

        KFParticle cluster = particle;

        std::vector<int> trackIds;
        bool isValidParticle = true;
        for(int iD=0; iD<particle.NDaughters(); iD++) {
          const int daughterId = particle.DaughterIds()[iD];
          const KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterId];
          const int daughterTrackId = daughter.DaughterIds()[0];
          trackIds.push_back(daughterTrackId);
          isValidParticle &= isValidTrack[daughterTrackId];
        }
//         if(!isValidParticle) continue;

        for(int iTrack=0; iTrack<nTracks; iTrack++) {
          KFParticle track = fStKFParticleInterface->GetParticles()[iTrack];
          
          //check that not the same track
          bool isSameTrack = false;
          for(unsigned int iD=0; iD<trackIds.size(); iD++)
            if(trackIds[iD] == track.DaughterIds()[0])
              isSameTrack = true;
          if(isSameTrack) continue;  
          
          // use only secondary tracks
          if(particle.NDaughters() == 2) {
            const float chiprimCut = 8;
            if(track.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex()) < chiprimCut) continue;
          }
          
          // should be close to the initial position
          const float chiSec2 = track.GetDeviationFromVertex(particle);
          if(chiSec2 > 18) continue;

          const float dev = track.GetDeviationFromVertex(cluster);
          if( dev > 10. ) continue;
//           if(particle.NDaughters() == 2) {
//             if( dev > 3. ) continue;
//           }
//           else {
//             if( dev > 10. ) continue;
//           }

          //add track to cluster
          KFParticle clusterTmp = cluster;
          clusterTmp += track;
          if(clusterTmp.GetChi2()/float(clusterTmp.GetNDF()) < 3) {
            cluster = clusterTmp;
            trackIds.push_back(track.DaughterIds()[0]);
          }
        }  //iTrack
        
        if(cluster.NDaughters() > 5) {
          
//           float l, dl;
//           particle.GetDistanceToVertexLine(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex(), l, dl);
//           std::cout << "event " << iEvent << "   pdg: " << particle.GetPDG() << " l " << l <<" l/dl " << (l/dl) << "  killed by cluster" << std::endl;
        
          fStKFParticleInterface->RemoveParticle(iParticle);
          for(unsigned int iTrackId=0; iTrackId<trackIds.size(); iTrackId++)
            isValidTrack[trackIds[iTrackId]] = false;
        }
      } //PDG cut 3001<pDG<=3029
    } //clean clusters iParticle loop over reconstructed particles

#if 1 //FIXME
    for(int iParticle=nTracks; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++) {
      const KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
      
      if( (abs(particle.GetPDG()) > 3001) && (abs(particle.GetPDG()) <= 3103) ) {

        const float dmCut = (particle.NDaughters() == 2) ? 2.5e-3f : 2.0e-3f;
        if(particle.GetErrMass() > dmCut) {
          fStKFParticleInterface->RemoveParticle(iParticle);
        }
#if 0
        bool isValidParticle = true;
        for(int iD=0; iD<particle.NDaughters(); iD++) {
          const int daughterId = particle.DaughterIds()[iD];
          const KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterId];
          const int daughterTrackId = daughter.DaughterIds()[0];
          isValidParticle &= isValidTrack[daughterTrackId];
        }
        if(!isValidParticle) {
          fStKFParticleInterface->RemoveParticle(iParticle);
          continue;
        }
#endif
        float l, dl;
        particle.GetDistanceToVertexLine(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex(), l, dl);
        if(dl > 3.f) {
          fStKFParticleInterface->RemoveParticle(iParticle);
          continue;
        }
      }
    }
#endif
#endif




#if 1 //FIXME
//clean primary lambdas
    const int nParticles0 = fStKFParticleInterface->GetParticles().size();
    for(int iParticle=0; iParticle<nParticles0; iParticle++) {
      KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];

      if(particle.GetPDG()==3006 || 
         particle.GetPDG()==3007 || 
         particle.GetPDG()==3012 || 
         particle.GetPDG()==3013 ||
         particle.GetPDG()==3028 ||
         particle.GetPDG()==3029)
      {
        KFParticle pion = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
        KFParticle fragment = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
        KFParticle proton = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[2]];
        
        KFParticle ppi;
        ppi += pion;
        ppi += proton;
        ppi.SetNonlinearMassConstraint(1.115683);
        const float chiPrimPPi = ppi.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
        
        const float chiPrimCutPPi = ((particle.GetPDG()==3012) || (particle.GetPDG()==3013)) ? 18.f : 8.f;
        
        if(chiPrimPPi < chiPrimCutPPi) {
          fStKFParticleInterface->RemoveParticle(iParticle);
          continue;
        }
        
        KFParticle lambdaFragment;
        lambdaFragment += fragment;
        lambdaFragment += ppi;
        
        float l, dl;
        lambdaFragment.GetDistanceToVertexLine(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex(), l, dl);
//         if((l/dl < 3) && (ppi.GetChi2()/float(ppi.GetNDF() < 10))) {
//           fStKFParticleInterface->RemoveParticle(iParticle);
//           continue;
//         }
        if(l/dl < 3) {
          fStKFParticleInterface->RemoveParticle(iParticle);
          continue;
        }

//         const float chiPrimCutFPi = ((particle.GetPDG()==3006) || (particle.GetPDG()==3007)) ? 3.f : 8.f;
//         KFParticle fpi;
//         fpi += pion;
//         fpi += fragment;
//         const float chiPrimFPi = fpi.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
// 
//         if(chiPrimFPi < chiPrimCutFPi) {
//           fStKFParticleInterface->RemoveParticle(iParticle);
//           continue;
//         }
//         const float chiF = fragment.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
//         if(chiF < 6 && chiPrimPPi < 18) {
//           fStKFParticleInterface->RemoveParticle(iParticle);
//           continue;
//         }
//         if( (chiPrimPPi < 18.f) && (chiF < 18.f) ) {
//           fStKFParticleInterface->RemoveParticle(iParticle);
//           continue;
//         }
      }
    }  //clean primary lambdas
#endif

    
    if(fStoreCandidates) {
      KFPartEfficiencies parteff;
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++) {
        const KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
        if(particle.GetPDG() == -1) continue;
        const int particleIndex = parteff.GetParticleIndex(particle.GetPDG());

        if(!fIsStoreCandidate[particleIndex]) continue;

        fPartcileCandidate = fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex();
        fCandidatesTree->Fill();

        fPartcileCandidate = particle;
        fCandidatesTree->Fill();

        if(particle.NDaughters() == 1) continue;
        
        for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++) {
          const KFParticle daughter = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[iDaughter]];
          fPartcileCandidate = daughter;
          fCandidatesTree->Fill();
        }
      }
    }

    int eventId = -1;
    int runId = -1;
    
   
    
    int centralityWeight = 1;
    int centralityBin= 1;
    
    fStKFParticlePerformanceInterface->SetMCTracks(mcTracks);
    fStKFParticlePerformanceInterface->SetMCIndexes(mcIndices);    
    fStKFParticlePerformanceInterface->SetCentralityBin(centralityBin);
    fStKFParticlePerformanceInterface->SetCentralityWeight(centralityWeight);
    Int_t nevent = 100000;
    fStKFParticlePerformanceInterface->SetPrintEffFrequency(nevent);
    fStKFParticlePerformanceInterface->PerformanceAnalysis();
    
   
   if(fKaonAnalysis) Fill_KaonNtuples();
  
  
  
  return kStOK;
}

//------------------------------------------------------------
bool StKFParticleAnalysisMaker::FillKFDaughters(KFParticle& particle){
  const KFParticleTopoReconstructor *topoRec=fStKFParticleInterface->GetTopoReconstructor();
 
  StThreeVectorD mother_p_decay(fK.mother_px,fK.mother_py,fK.mother_pz);
  StThreeVectorD mother_p_PVX(fK.mother_px_PVX,fK.mother_py_PVX,fK.mother_pz_PVX);
 
      //now fill the daughter tracks  

      for(int iD=0; iD<particle.NDaughters(); iD++) {
              TDaughter &daughter=fK.daughter(iD);
              daughter.index=iD;
              int id_kf=particle.DaughterIds()[iD];
              KFParticle daugh=(fStKFParticleInterface->GetParticles())[id_kf];
              daughter.id=daugh.DaughterIds()[0]; //id of parent MuDST tracks

             
              daughter.pt = daugh.GetPt();
              daughter.p  = daugh.GetP();
              daughter.pdg = daugh.GetPDG();
              daughter.eta = daugh.GetEta();
              daughter.phi = daugh.GetPhi();
              daughter.px = daugh.GetPx();
              daughter.py = daugh.GetPy();
              daughter.pz = daugh.GetPz();
             
              StMuTrack *mutrack=NULL;
              StPicoTrack *picotrack=NULL;

              const int iDataTrack = fStKFParticleInterface->TrackIdToI()[daughter.id]; //also index forpico track


              if(fIsPicoAnalysis){
               picotrack= fPicoDst->track(iDataTrack);
               daughter.charge= picotrack->charge();
               daughter.nhits=picotrack->nHitsFit(); 
               daughter.nhits_dEdx=picotrack->nHitsDedx(); 
               daughter.nhits_pos=picotrack->nHitsPoss(); 
               daughter.dEdx=picotrack->dEdx();
               daughter.idTruth =picotrack->idTruth();
               daughter.qaTruth =picotrack->qaTruth();

                //last point radius - workaround from topology map
               /// Toplogy Map data0 and data1. See StEvent/StTrackTopologyMap.cxx
             // daughter.lastPointR=picotrack->lastPoint().perp(); //from PV to last hist
               daughter.lastPointR=0;
            
               #if !defined (__TFG__VERSION__)
               StTrackTopologyMap map(picotrack->topologyMap(0),picotrack->topologyMap(1),picotrack->iTpcTopologyMap());
               #else
               StTrackTopologyMap map(picotrack->topologyMap(0),picotrack->topologyMap(1),picotrack->topologyMap(2));
               #endif  
               if (map.hasHitInDetector(kTpcId)){
                int lastTpcRow=45;
                while (lastTpcRow>0 && !map.hasHitInRow(kTpcId, lastTpcRow)){lastTpcRow--;}
                 if (lastTpcRow>0) {
                  if (lastTpcRow<=13) daughter.lastPointR=St_tpcPadPlanesC::instance()->innerRowRadii(0)[lastTpcRow-1];
                  else daughter.lastPointR=St_tpcPadPlanesC::instance()->outerRowRadii(0)[lastTpcRow-14];
                }
               }//hits in TPC
               //now iTPC
               if (map.hasHitInDetector(kiTpcId)){
                float Ri=0;
                int lastiTpcRow=40;
                  while ( lastiTpcRow>0 && !map.hasHitInRow(kiTpcId, lastiTpcRow)) {lastiTpcRow--;}
                  if (lastiTpcRow>0) Ri=St_itpcPadPlanesC::instance()->innerRowRadii(0)[lastiTpcRow-1];
                  if (Ri>daughter.lastPointR) daughter.lastPointR=Ri;
               }//iTPC hits
              

               //St_tpcPadConfigC::instance()->
             }//picoDst
              else {
                mutrack = (StMuTrack *) fMuDst->globalTracks(iDataTrack); 
                daughter.charge= mutrack->charge();
                daughter.nhits=mutrack->nHitsFit(); 
                daughter.nhits_dEdx=mutrack->nHitsDedx(); 
                daughter.nhits_pos=mutrack->nHitsPoss(); 
                daughter.lastPointR=mutrack->lastPoint().perp(); //from PV to last hist
                daughter.dEdx=mutrack->dEdx();
                daughter.idTruth =mutrack->idTruth();
                daughter.qaTruth =mutrack->qaTruth();

            }

              cout<<" fK.mother_isMc="<<fK.mother_isMc<<endl;
              cout<<" iD="<<iD<<" .. daughter.idTruth="<<daughter.idTruth<<endl;
              
              fK.mother_isMc=fK.mother_isMc && (daughter.idTruth>0)&&(daughter.idTruth<10000); //above 10000 it comes from real data
              cout<<" fK.mother_isMc="<<fK.mother_isMc<<endl;
             
              //DCA's - mainly useful for mother track - iD==3
              //a) to decay Vtx: from KFP, from MuTrack
                 //from KFP
              daughter.DecayDca_KF = daugh.GetDistanceFromVertex(particle);
              KFParticle tmp=daugh;
              tmp.TransportToParticle(particle);
              StThreeVectorD daugher_p_decay_KF(tmp.Px(),tmp.Py(),tmp.Pz());
              daughter.dp_decay_KF=(daugher_p_decay_KF-mother_p_decay).mag();

              StThreeVectorD decayVtx(fK.decay_Vx, fK.decay_Vy, fK.decay_Vz); 
               
              //muDST
              StThreeVectorD daughter_p_decay;
              StThreeVectorD daughter_p_PVX;
              
              if(fIsPicoAnalysis){
                 StPicoPhysicalHelix helix=picotrack->helix(fPicoDst->event()->bField()); //no kilogauss here!!

                double cx=helix.xcenter();double cy=helix.ycenter();
                 daughter.helix_R=1./helix.curvature(); daughter.helix_Cr=sqrt(cx*cx+cy*cy);
                 daughter.helix_lowR=daughter.helix_Cr-daughter.helix_R;
                 daughter.helix_hiR=daughter.helix_Cr+daughter.helix_R;
                 TVector3 decayVtx_(fK.decay_Vx, fK.decay_Vy, fK.decay_Vz); //for picoDst analyses
                 double pathlength = helix.pathLength(decayVtx_, true ); // false- do not scan periods
                 daughter.DecayDca_mu=helix.distance(decayVtx_);//(helix.at(pathlength)-decayVtx).mag();
                //momentum at decay vtx - from global track!
                TVector3 tmp= helix.momentumAt(pathlength,fPicoDst->event()->bField()*kilogauss);
                daughter_p_decay.set(tmp.X(),tmp.Y(),tmp.Z());
                daughter.dp_Decay=(daughter_p_decay-mother_p_decay).mag();

                //B)to prim Vtx:
                //from KFP
                 daughter.PvtxDca_KF = daugh.GetDistanceFromVertex(topoRec->GetPrimVertex());
                //from picoTrack                
                 TVector3 pVtx_=fK.EvInfo().primVtx_TVec();
                 pathlength = helix.pathLength(pVtx_, true );
                 daughter.PvtxDca_mu=helix.distance(pVtx_);//(helix.at(pathlength)-pVtx).mag();
                 tmp= helix.momentumAt(pathlength,fPicoDst->event()->bField()*kilogauss);
                 daughter_p_PVX.set(tmp.X(),tmp.Y(),tmp.Z());
                 daughter.dp_PVX=(daughter_p_PVX-mother_p_PVX).mag();

                //DCA from MuDST
                daughter.PvtxDca_official=picotrack->gDCA(fPicoDst->event()->primaryVertex()).Mag(); 

               

              }
               else { //MuDst
                StPhysicalHelixD helix= mutrack->helix();

                double cx=helix.xcenter();double cy=helix.ycenter();
                daughter.helix_R=1./helix.curvature(); daughter.helix_Cr=sqrt(cx*cx+cy*cy);
                daughter.helix_lowR=daughter.helix_Cr-daughter.helix_R;
                daughter.helix_hiR=daughter.helix_Cr+daughter.helix_R;


                StThreeVectorD decayVtx(fK.decay_Vx, fK.decay_Vy, fK.decay_Vz); 
                double pathlength = helix.pathLength(decayVtx, true ); // false- do not scan periods
                 daughter.DecayDca_mu=helix.distance(decayVtx);//(helix.at(pathlength)-decayVtx).mag();
                //momentum at decay vtx - from global track!
                 daughter_p_decay= helix.momentumAt(pathlength,fMuDst->event()->runInfo().magneticField()*kilogauss);
                 daughter.dp_Decay=(daughter_p_decay-mother_p_decay).mag();

                //B)to prim Vtx:
                //from KFP
                 daughter.PvtxDca_KF = daugh.GetDistanceFromVertex(topoRec->GetPrimVertex());
                //from MuTrack 
                 StThreeVectorD pVtx=fK.EvInfo().primVtx_StVec();
                 pathlength = helix.pathLength(pVtx, true );
                 daughter.PvtxDca_mu=(helix.at(pathlength)-pVtx).mag();
              //daughter.PvtxDca_mu=helix.distance(pVtx); o
                 daughter_p_PVX= helix.momentumAt(pathlength,fMuDst->event()->runInfo().magneticField()*kilogauss);
                 daughter.dp_PVX=(daughter_p_PVX-mother_p_PVX).mag();

                 //DCA from MuDST
                 daughter.PvtxDca_official=mutrack->dcaGlobal().mag();   
 
              }
                    

              daughter.decay_dl=daugh.GetDeviationFromVertex(particle);

              daughter.decay_p=daughter_p_decay.mag();daughter.decay_pt=daughter_p_decay.perp(); daughter.decay_eta=daughter_p_decay.pseudoRapidity();
              daughter.decay_phi=daughter_p_decay.phi();
              daughter.decay_px=daughter_p_decay.x();daughter.decay_py=daughter_p_decay.y();daughter.decay_pz=daughter_p_decay.z();
              
              /*
              double x=daughter.decay_px*fK.decay_Vx + daughter.decay_py*fK.decay_Vy;
              double y=daughter.decay_px*fK.decay_Vy - daughter.decay_py*fK.decay_Vx;
             daughter.phi_wrt_mother=TMath::ATan2(y,x); 
              */
              //better this way
              daughter.phi_wrt_Vr=decayVtx.angle(daughter_p_decay);

             } //NDaughters loop
  return true;
}

void StKFParticleAnalysisMaker::Fill_KaonNtuples() {
  cout<<"StKFParticleAnalysisMaker::Fill_KaonNtuples() fIsPicoAnalysis="<<fIsPicoAnalysis<<endl;
 
  const int nPIDS=4;
  int particlesPDG[nPIDS] = {100321, 200321,-100321, -200321}; //K->3pi with found K (200321), K->3pi only (100321)
  int nDaughters[nPIDS] = {3, 4, 3, 4};

//Note There can be mutiple primary vertexes. I should take only the  one that was intarfaced to KFP
// the 3pi should actually point "somewhere" around this vertex - I do not want potential kaons from other vertexes
// and then match the primary kaon

 const KFParticleTopoReconstructor *topoRec=fStKFParticleInterface->GetTopoReconstructor();
 // cout<<" Primary vertices:"<<endl<<"in KF: "<<topoRec->NPrimaryVertices()<< " , ID="<<topoRec->GetPrimVertex().Id()<<endl;
  //cout<<"in MuDst: "<<fMuDst->numberOfPrimaryVertices()<<endl;
  int npt=fStKFParticleInterface->GetParticles().size();
  cout<<"fStKFParticleInterface->GetParticles().size()="<<npt<<endl;
  /*test 
  for(unsigned int iParticle=0; iParticle<npt; iParticle++) {
    cout<<"iParticle="<<iParticle<<" ID="<<fStKFParticleInterface->GetParticles()[iParticle].Id()<<endl;
  }  
*/

  //primary vertex
  KFParticle primVtx=topoRec->GetPrimVertex();
 
  for(unsigned int iParticle=0; iParticle<npt; iParticle++) {
    //cout<<"i="<<iParticle<<" KFP_iD="<<fStKFParticleInterface->GetParticles()[iParticle].Id()
    //<<" PDG="<<fStKFParticleInterface->GetParticles()[iParticle].GetPDG()<<" pt="<<
    //fStKFParticleInterface->GetParticles()[iParticle].GetPt()<<endl;
    //cout<< fStKFParticleInterface->GetParticles()[iParticle]<<endl;
  
    for(int iPDG = 0; iPDG < nPIDS; iPDG++){
  
      if(fStKFParticleInterface->GetParticles()[iParticle].GetPDG() != particlesPDG[iPDG]) continue;
      cout<<"found 3pi from KFParticle PDG="<<fStKFParticleInterface->GetParticles()[iParticle].GetPDG()<<endl;

      KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];

      if(particle.NDaughters() != nDaughters[iPDG]) {
        cout << "Wrong number of daughters! for"<< particlesPDG[iPDG] << " expected "<< nDaughters[iPDG]<<" but found "<< particle.NDaughters()<<endl;
        return false;
      }

      //start fillig
      fK.Clear();

      //fill event info
      if(fIsPicoAnalysis) fK.EvInfo().Fill(fPicoDst,primVtx); else fK.EvInfo().Fill(fMuDst,primVtx);   

      //fill mother (3pi vertex)  info
      // skip decays out of TPC --- save disk space
      if (particle.GetR()<50){ cout<<" ..skipping K.decay_Vr="<<fK.decay_Vr<<endl; continue;} 


      fK.mother_PID=particle.GetPDG();
      fK.mother_m = particle.GetMass();

        //decay point
      fK.decay_Vx=particle.GetX();
      fK.decay_Vy=particle.GetY();
      fK.decay_Vz=particle.GetZ();
      fK.decay_Vr=particle.GetR();
      cout<<"fK.decay_Vr="<<fK.decay_Vr<<endl;
     
     
      //momentum at decay point
      fK.mother_pt = particle.GetPt();
      fK.mother_px = particle.GetPx();
      fK.mother_py = particle.GetPy();
      fK.mother_pz = particle.GetPz();
      fK.mother_eta = particle.GetEta();
      fK.mother_phi = particle.GetPhi();          

      // move to primary vertex
      particle.TransportToPoint(primVtx.Parameters());
      //momentum at primary vertex
      fK.mother_px_PVX = particle.GetPx();
      fK.mother_py_PVX = particle.GetPy();
      fK.mother_pz_PVX = particle.GetPz();
      fK.mother_pt_PVX = particle.GetPt();
      fK.mother_eta_PVX = particle.GetEta();
      fK.mother_phi_PVX = particle.GetPhi();
      

       
         //distance to PV - usefull to see if it comes from the primary vertex...I coudl have seom decays of secondary kaons
        //float_v lCandidate, dlCandidate;
         //KFParticleSIMD part(particle);
        //KFParticleSIMD pvt(topoRec->GetPrimVertex());
        //daugh.GetDistanceToVertexLine(pvt, lCandidate, dlCandidate);
/*
  const float chiPrimPPi = ppi.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
        
        const float chiPrimCutPPi = ((particle.GetPDG()==3012) || (particle.GetPDG()==3013)) ? 18.f : 8.f;
        
        if(chiPrimPPi < chiPrimCutPPi) {
          fStKFParticleInterface->RemoveParticle(iParticle);
          continue;
        }
        
        KFParticle lambdaFragment;
        lambdaFragment += fragment;
        lambdaFragment += ppi;
        
        float l, dl;
        lambdaFragment.GetDistanceToVertexLine(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex(), l, dl);
//         if((l/dl < 3) && (ppi.GetChi2()/float(ppi.GetNDF() < 10))) {
//           fStKFParti
*/
     // fK.mother_PV_chi2 =particle.GetDeviationFromVertex(primVtx);
      float l, dl;
      particle.GetDistanceToVertexLine(primVtx,l,dl);
      fK.mother_PV_l=l; fK.mother_PV_dl=dl; 
      //wrong... fK.mother_PV_l=particle.GetDistanceFromVertex(primVtx);
     
     //fill daughters
     fK.mother_isMc=1; //will changed by FillDaughters
     if (!FillKFDaughters(particle)) continue;

     
     // matching of tracks from MuDST to the found decay vertex
    
     if(fabs(fStKFParticleInterface->GetParticles()[iParticle].GetPDG())!=100321) goto FILL_TREE;// run only decay vertices
 
    
     MatchMotherKaon(particle);

     
  FILL_TREE:
 
     fKaonTree->Fill();
     
    }//nPIDs loop
       
  }//FKparticles loop
}

void StKFParticleAnalysisMaker::MatchMotherKaon(KFParticle& particle){
  StThreeVectorD mother_p_decay(fK.mother_px,fK.mother_py,fK.mother_pz);
  StThreeVectorD mother_p_PVX(fK.mother_px_PVX,fK.mother_py_PVX,fK.mother_pz_PVX);


  if (fIsPicoAnalysis){ //for picoDstAnalysis
      cout<<"matching.. from pico"<<endl;
      StPicoTrack *picotrack;
      TDaughter best_dt;
      best_dt.Clear();
      float bestDca=10e20; //ket track of best dca value
      StThreeVectorD best_dt_p_PVX;
      best_dt.dp_Decay=10e20;
      best_dt.DecayDca_mu=10e20;
      //loop over global

      TDaughter &daughter =fK.daughter(4);
        
      for (UInt_t k = 0; k < fPicoDst->numberOfTracks(); k++) {
           picotrack = fPicoDst->track(k);
           if (! picotrack) continue;
          daughter.Clear();
          //DCA from MuDST
          daughter.PvtxDca_official=picotrack->gDCA(fPicoDst->event()->primaryVertex()).Mag();
          if (daughter.PvtxDca_official>10 ) continue;//junk

          daughter.index=4;
          daughter.id=picotrack->id();

          daughter.charge=picotrack->charge();
          daughter.pt = picotrack->gPt();
          daughter.p  = picotrack->gPtot();
          daughter.eta = picotrack->gMom().Eta();
          daughter.phi = picotrack->gMom().Phi();
          daughter.px = picotrack->gMom().X();
          daughter.py = picotrack->gMom().Y();
          daughter.pz = picotrack->gMom().z();
          daughter.nhits=picotrack->nHitsFit(); 
          daughter.nhits_dEdx=picotrack->nHitsDedx();
          daughter.nhits_pos=picotrack->nHitsPoss(); 
          daughter.dEdx=picotrack->dEdx(); 
          daughter.idTruth =picotrack->idTruth();
          daughter.qaTruth =picotrack->qaTruth();

        
           StPicoPhysicalHelix helix = picotrack->helix(fPicoDst->event()->bField());
           TVector3 decayVtx_(fK.decay_Vx, fK.decay_Vy, fK.decay_Vz); //for picoDst analyses
           StThreeVectorD decayVtx(fK.decay_Vx, fK.decay_Vy, fK.decay_Vz); 
           double pathlength = helix.pathLength(decayVtx_, true ); // false- do not scan periods
           //daughter.pdg=pathlength;// ok, that's dirty...
           daughter.DecayDca_mu=helix.distance(decayVtx_);//(helix.at(pathlength)-decayVtx).mag();
        
           //momentum at decay vtx - from global track!
           TVector3 tmp= helix.momentumAt(pathlength,fPicoDst->event()->bField()*kilogauss);
           StThreeVectorD mother_p_decay(fK.mother_px,fK.mother_py,fK.mother_pz);
           StThreeVectorD daughter_p_decay; daughter_p_decay.set(tmp.X(),tmp.Y(),tmp.Z());
           daughter.dp_Decay=(daughter_p_decay-mother_p_decay).mag();



          //pathlength = helix.pathLength(pVtx, true );
          //daughter.PvtxDca_mu=(helix.at(pathlength)-pVtx).mag();
          TVector3 pVtx_(fK.EvInfo().Vx, fK.EvInfo().Vy, fK.EvInfo().Vz);
          daughter.PvtxDca_mu=helix.distance(pVtx_); 
             //daughter.DecayDca_mu=fabs(helix.geometricSignedDistance(v));
          tmp= helix.momentumAt(pathlength,fPicoDst->event()->bField()*kilogauss);
          StThreeVectorD daughter_p_PVX;
          daughter_p_PVX.set(tmp.X(),tmp.Y(),tmp.Z());
          daughter.dp_PVX=(daughter_p_PVX-mother_p_PVX).mag();
          
          if (daughter.dp_Decay>2.) continue; //junk                
          if (daughter.DecayDca_mu>2.) continue; //junk
 
          cout<<" match ok"<<endl;
          
            //daughter.lastPointR=picotrack->lastPoint().perp(); //from PV to last hist
          //this is special for picoDST since the last point is not saved - we need to go via hitmap
          daughter.lastPointR=0; //from PV to last hist
          #if !defined (__TFG__VERSION__)
          StTrackTopologyMap map(picotrack->topologyMap(0),picotrack->topologyMap(1),picotrack->iTpcTopologyMap());
          cout<<"topo NO-TFG data: "<<picotrack->topologyMap(0)<<" "<<picotrack->topologyMap(1)<<" "<<picotrack->iTpcTopologyMap()<<endl;
          #else
          StTrackTopologyMap map(picotrack->topologyMap(0),picotrack->topologyMap(1),picotrack->topologyMap(2));
         cout<<"topo TFG data: "<<picotrack->topologyMap(0)<<" "<<picotrack->topologyMap(1)<<" "<<picotrack->topologyMap(2)<<endl;
          #endif  
          if (map.hasHitInDetector(kTpcId)){
            int lastTpcRow=45;
            while ( lastTpcRow>0 && !map.hasHitInRow(kTpcId, lastTpcRow)) {lastTpcRow--;}
            cout<<"last tpc row="<<lastTpcRow;
            if (lastTpcRow>0) {
              if (lastTpcRow<=13) daughter.lastPointR=St_tpcPadPlanesC::instance()->innerRowRadii(0)[lastTpcRow-1];
              else daughter.lastPointR=St_tpcPadPlanesC::instance()->outerRowRadii(0)[lastTpcRow-14];
            }
            cout<<" R="<<daughter.lastPointR<<endl;
           }//hits in TPC
           if (map.hasHitInDetector(kiTpcId)){
            float Ri=0;
            int lastiTpcRow=40;
            while ( lastiTpcRow>0 && !map.hasHitInRow(kiTpcId, lastiTpcRow)) {lastiTpcRow--;}
            cout<<"last iTpc row="<<lastiTpcRow;
            if (lastiTpcRow>0) Ri=St_itpcPadPlanesC::instance()->innerRowRadii(0)[lastiTpcRow-1];
            if (Ri>daughter.lastPointR) daughter.lastPointR=Ri;
            cout<<" Ri="<<Ri<<endl;
          }//iTPC hits
         
         
           //if (daughter.dp_Decay>1000.) continue; //junk from initialization best_dt.dp_Decay=10e20;
          
           
          //daughter.decay_dl=daugh.GetDeviationFromVertex(particle);
          
          daughter.decay_p=daughter_p_decay.mag();daughter.decay_pt=daughter_p_decay.perp(); daughter.decay_eta=daughter_p_decay.pseudoRapidity();
          daughter.decay_phi=daughter_p_decay.phi();
          daughter.decay_px=daughter_p_decay.x();daughter.decay_py=daughter_p_decay.y();daughter.decay_pz=daughter_p_decay.z();

          
          /*
          double x=daughter.decay_px*fK.decay_Vx + daughter.decay_py*fK.decay_Vy;
          double y=daughter.decay_px*fK.decay_Vy - daughter.decay_py*fK.decay_Vx;
         daughter.phi_wrt_mother=TMath::ATan2(y,x); 
          */
          daughter.phi_wrt_Vr=decayVtx.angle(daughter_p_decay);


          cout<<" dp_Decay="<<daughter.dp_Decay<<" best now="<<best_dt.dp_Decay<<endl;
          if (bestDca>daughter.DecayDca_mu) bestDca=daughter.DecayDca_mu;
          if (best_dt.dp_Decay>daughter.dp_Decay){
            cout<<" new best"<<endl;
            daughter.isBest=1;
            //now swap the best
            TDaughter tmpd=best_dt;   best_dt=daughter;  daughter=tmpd;
            StThreeVectorD tmpp=best_dt_p_PVX; best_dt_p_PVX=daughter_p_decay; daughter_p_decay=tmpp;
            daughter.isBest=0;
            cout<<"best so far .."<< best_dt.dp_Decay<<endl;
          }
          if (best_dt.DecayDca_mu>bestDca) best_dt.isBest=2; //candidate based on dca would be different
         // now in daughter there is the second best

         
                
         } //loop over primary    

         //save the best
         
         daughter= best_dt;
        /* skipping so far

        //now save the best
        if (best_dt.dp_Decay<= 1000.){
         //if (best_dt.DecayDca_mu<= 20.){
           Float_t toFill[] = { runId,eventId,Vz, mother_PID, mother_pt, mother_px,mother_py,mother_pz, mother_eta, mother_phi, 
            mother_pt_PVX, mother_px_PVX, mother_py_PVX, mother_pz_PVX, mother_eta_PVX, mother_phi_PVX, mother_m,  mother_PV_l, mother_PV_dl, mother_isMc,
            decay_Vr, decay_Vx, decay_Vy, decay_Vz,
            best_dt.index, best_dt.pdg, best_dt.idTruth, best_dt.qaTruth, best_dt.nhits,best_dt.nhits_pos,best_dt.lastPointR,best_dt.nhits_dEdx, best_dt.dEdx, best_dt.p,best_dt.pt, 
            best_dt.eta,best_dt.phi,best_dt.px,best_dt.py,best_dt.pz,best_dt.dp_PVX,
            0,best_dt_p_PVX.perp(),best_dt_p_PVX.pseudoRapidity(),best_dt_p_PVX.phi(),best_dt_p_PVX.x(),best_dt_p_PVX.y(),best_dt_p_PVX.z(),best_dt.dp_Decay,0,
            best_dt.DecayDca_KF,best_dt.DecayDca_mu,best_dt.PvtxDca_KF,best_dt.PvtxDca_mu,best_dt.PvtxDca_official,best_dt.isBest};
            hKaonNtuple->Fill(toFill);        
        }
        */
       } //if (!fIsPico)


     if (!fIsPicoAnalysis){ //for MuDstAnalysis
      cout<<"matching.. from muDst"<<endl;
      StMuTrack *mutrack;
      TDaughter best_dt;
      StThreeVectorD best_dt_p_PVX;
      float bestDca=10e20;
      best_dt.dp_Decay=10e20;
      best_dt.DecayDca_mu=10e20;
      //loop over global
      TDaughter &daughter =fK.daughter(4);
      for (UInt_t k = 0; k < fMuDst->numberOfGlobalTracks(); k++) {
       mutrack = (StMuTrack *) fMuDst->array(muGlobal)->UncheckedAt(k);
      //lglob
      //for (UInt_t k = 0; k < fMuDst->numberOfPrimaryTracks(); k++) {
      //    track = (StMuTrack *) fMuDst->array(muPrimary)->UncheckedAt(k);
          if (! mutrack) continue;
         
          TDaughter &daughter =fK.daughter(4);
          daughter.Clear();
         
         
          //DCA from MuDST
          daughter.PvtxDca_official=mutrack->dcaGlobal().mag();
          if (daughter.PvtxDca_official>10 ) continue;//junk


          daughter.index=4;
          daughter.id=mutrack->id();

          daughter.charge=mutrack->charge();
          daughter.pt = mutrack->pt();
          daughter.p  = mutrack->p().mag();
          daughter.eta = mutrack->eta();
          daughter.phi = mutrack->phi();
          daughter.px = mutrack->p().x();
          daughter.py = mutrack->p().y();
          daughter.pz = mutrack->p().z();
          daughter.nhits=mutrack->nHitsFit(); 
          daughter.nhits_dEdx=mutrack->nHitsDedx();
          daughter.nhits_pos=mutrack->nHitsPoss(); 
          daughter.lastPointR=mutrack->lastPoint().perp(); //from PV to last hist
          daughter.dEdx=mutrack->dEdx(); 
          daughter.idTruth =mutrack->idTruth();
          daughter.qaTruth =mutrack->qaTruth();


           StPhysicalHelixD helix = mutrack->helix();
           TVector3 decayVtx_(fK.decay_Vx, fK.decay_Vy, fK.decay_Vz); //for picoDst analyses
           StThreeVectorD decayVtx(fK.decay_Vx, fK.decay_Vy, fK.decay_Vz); 
           double pathlength = helix.pathLength(decayVtx, true ); // false- do not scan periods
           //daughter.pdg=pathlength;// ok, that's dirty...
           daughter.DecayDca_mu=helix.distance(decayVtx);//(helix.at(pathlength)-decayVtx).mag();
            
         
           //momentum at decay vtx - from global track!
           StThreeVectorD daughter_p_decay= helix.momentumAt(pathlength,fMuDst->event()->runInfo().magneticField()*kilogauss);
           daughter.dp_Decay=(daughter_p_decay-mother_p_decay).mag();

        
          //pathlength = helix.pathLength(pVtx, true );
          //daughter.PvtxDca_mu=(helix.at(pathlength)-pVtx).mag();
          StThreeVectorD pVtx(fK.EvInfo().Vx, fK.EvInfo().Vy, fK.EvInfo().Vz);
          daughter.PvtxDca_mu=helix.distance(pVtx); 
             //daughter.DecayDca_mu=fabs(helix.geometricSignedDistance(v));
          StThreeVectorD daughter_p_PVX= helix.momentumAt(pathlength,fMuDst->event()->runInfo().magneticField()*kilogauss);
          daughter.dp_PVX=(daughter_p_PVX-mother_p_PVX).mag();

          if (daughter.dp_Decay>2.) continue; //junk                
          if (daughter.DecayDca_mu>2.) continue; //junk

          cout<<" match ok"<<endl;
        

           //if (daughter.dp_Decay>1000.) continue; //junk from initialization best_dt.dp_Decay=10e20;
          
          //daughter.decay_dl=daugh.GetDeviationFromVertex(particle);

          daughter.decay_p=daughter_p_decay.mag();daughter.decay_pt=daughter_p_decay.perp(); daughter.decay_eta=daughter_p_decay.pseudoRapidity();
          daughter.decay_phi=daughter_p_decay.phi();
          daughter.decay_px=daughter_p_decay.x();daughter.decay_py=daughter_p_decay.y();daughter.decay_pz=daughter_p_decay.z();

          /*
          double x=daughter.decay_px*fK.decay_Vx + daughter.decay_py*fK.decay_Vy;
          double y=daughter.decay_px*fK.decay_Vy - daughter.decay_py*fK.decay_Vx;
         daughter.phi_wrt_mother=TMath::ATan2(y,x); 
        ` */
          daughter.phi_wrt_Vr=decayVtx.angle(daughter_p_decay);


         cout<<" dp_Decay="<<daughter.dp_Decay<<" best now="<<best_dt.dp_Decay<<endl;
          if (bestDca>daughter.DecayDca_mu) bestDca=daughter.DecayDca_mu;
          if (best_dt.dp_Decay>daughter.dp_Decay){
            cout<<" new best"<<endl;
                        daughter.isBest=1;
            //now swap the best
            TDaughter tmpd=best_dt;   best_dt=daughter;  daughter=tmpd;
            StThreeVectorD tmpp=best_dt_p_PVX; best_dt_p_PVX=daughter_p_decay; daughter_p_decay=tmpp;
            daughter.isBest=0;
            cout<<"best so far .."<< best_dt.dp_Decay<<endl;
          }
          if (best_dt.DecayDca_mu>bestDca) best_dt.isBest=2; //candidate based on dca would be different
          
          //daughter is now second best                   
         


        /*
          Float_t toFill[] = { runId,eventId,Vz, mother_PID, mother_pt, mother_px,mother_py,mother_pz, mother_eta, mother_phi, 
            mother_pt_PVX, mother_px_PVX, mother_py_PVX, mother_pz_PVX, mother_eta_PVX, mother_phi_PVX, mother_m,  mother_PV_l, mother_PV_dl, mother_isMc,
            decay_Vr, decay_Vx, decay_Vy, decay_Vz,
            daughter.index, daughter.pdg, daughter.idTruth, daughter.qaTruth,daughter.nhits,daughter.nhits_pos,daughter.lastPointR,daughter.nhits_dEdx, daughter.dEdx, daughter.p,daughter.pt, 
            daughter.eta,daughter.phi,daughter.px,daughter.py,daughter.pz,daughter.dp_PVX,
            0,daughter_p_decay.perp(),daughter_p_decay.pseudoRapidity(),daughter_p_decay.phi(),daughter_p_decay.x(),daughter_p_decay.y(),daughter_p_decay.z(),daughter.dp_Decay,0,
            daughter.DecayDca_KF,daughter.DecayDca_mu,daughter.PvtxDca_KF,daughter.PvtxDca_mu,daughter.PvtxDca_official,0};
            hKaonNtuple->Fill(toFill);        
            */
            } //loop over primary     
  
       daughter= best_dt;
       
          /* skipping so far
        //now save the best
        if (best_dt.dp_Decay<= 1000.){
         //if (best_dt.DecayDca_mu<= 20.){
           Float_t toFill[] = { runId,eventId,Vz, mother_PID, mother_pt, mother_px,mother_py,mother_pz, mother_eta, mother_phi, 
            mother_pt_PVX, mother_px_PVX, mother_py_PVX, mother_pz_PVX, mother_eta_PVX, mother_phi_PVX, mother_m,  mother_PV_l, mother_PV_dl,  mother_isMc,
            decay_Vr, decay_Vx, decay_Vy, decay_Vz,
            best_dt.index, best_dt.pdg,best_dt.idTruth, best_dt.qaTruth, best_dt.nhits,best_dt.nhits_pos,best_dt.lastPointR,best_dt.nhits_dEdx, best_dt.dEdx, best_dt.p,best_dt.pt, 
            best_dt.eta,best_dt.phi,best_dt.px,best_dt.py,best_dt.pz,best_dt.dp_PVX,
            0,best_dt_p_PVX.perp(),best_dt_p_PVX.pseudoRapidity(),best_dt_p_PVX.phi(),best_dt_p_PVX.x(),best_dt_p_PVX.y(),best_dt_p_PVX.z(),best_dt.dp_Decay,0,
            best_dt.DecayDca_KF,best_dt.DecayDca_mu,best_dt.PvtxDca_KF,best_dt.PvtxDca_mu,best_dt.PvtxDca_official,best_dt.isBest};
            hKaonNtuple->Fill(toFill);   
            *
        }*/
       } //if (!fIsPico)

}

void StKFParticleAnalysisMaker::GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle)
{
  if(particle.NDaughters() == 1)
  {
     int trackId = particle.DaughterIds()[0];   
    iDaughterTrack++;
  }
  else if(particle.NDaughters() > 1)
  {
    int order[4] = {0, 1, 2, 3};
    if( particle.GetPDG() == -421 || particle.GetPDG() == -411 || particle.GetPDG() == -431 ||   
        particle.GetPDG() == -429 || particle.GetPDG() == -4122) 
    { 
      order[0] = 1; 
      order[1] = 0; 
    }
    
    for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
    {
      const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
      KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
      //set pdg for correct order of cuts
      if(particle.GetPDG() == 521 && daughter.GetPDG() == -1) daughter.SetPDG(-421);
      if(particle.GetPDG() ==-521 && daughter.GetPDG() == -1) daughter.SetPDG( 421);
      if(particle.GetPDG() == 511 && daughter.GetPDG() == -1) daughter.SetPDG(-411);
      if(particle.GetPDG() ==-511 && daughter.GetPDG() == -1) daughter.SetPDG( 411);
        
      GetDaughterParameters(iReader, iDaughterTrack, iDaughterParticle, daughter);
    }
    
    
    KFParticleSIMD tempSIMDParticle(particle);
    float_v l,dl;
    KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
    
    tempSIMDParticle.SetProductionVertex(pv);
    iDaughterParticle++;
  }
}

/*
void StKFParticleAnalysisMaker::GetParticleParameters(const int iReader, KFParticle& particle)
{
  bool isBMeson = abs(particle.GetPDG()) == 511 || abs(particle.GetPDG()) == 521;
//   if( !isBMeson ) return;
  
  int iDaughterTrack = 0;
  int iDaughterParticle = 0;
  GetDaughterParameters(iReader, iDaughterTrack, iDaughterParticle, particle);

  int nDaughterParticleCut = 0;
  if(isBMeson) nDaughterParticleCut += 3;
  nDaughterParticleCut += fDaughterNames[iReader].size()*fNTrackTMVACuts;
  
  fTMVAParticleParameters[iReader][nDaughterParticleCut]   = particle.Chi2()/particle.NDF();  
  
  KFParticleSIMD tempSIMDParticle(particle);
  float_v l,dl;
  KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 1] = l[0]/dl[0];
  
  tempSIMDParticle.SetProductionVertex(pv);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 2] = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);

  if(fIsPicoAnalysis)
    fTMVAParticleParameters[iReader][nDaughterParticleCut + 3] = fPicoDst->event()->refMult();
  else
    fTMVAParticleParameters[iReader][nDaughterParticleCut + 3] = fMuDst->event()->refMult();
}
*/
Int_t StKFParticleAnalysisMaker::Finish() 
{
  
  if (fKaonFile) {
    fKaonFile->cd();
    fKaonTree->Write();
    fKaonFile->Close(); delete fKaonFile;}
  
  if(fStoreCandidates)
  {
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    fCandidateFile->cd();
    fCandidatesTree->Write();
    gFile = curFile;
    gDirectory = curDirectory;
  }
  
  return kStOK;
}

long StKFParticleAnalysisMaker::GetUniqueEventId(const int iRun, const int iEvent) const
{
  long id = 1000000000;
  return id*(iRun%1000) + iEvent;
}


void StKFParticleAnalysisMaker::AddDecayToReconstructionList( int iDecay ) { fDecays.push_back(iDecay); }

void StKFParticleAnalysisMaker::AddCandidateToStore(int pdg) {
  KFPartEfficiencies parteff;
  const int particleIndex = parteff.GetParticleIndex(pdg);
  fIsStoreCandidate[particleIndex] = true;
  fStoreCandidates = true;
}
