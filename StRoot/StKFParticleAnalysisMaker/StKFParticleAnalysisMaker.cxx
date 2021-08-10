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
//--- TMVA classes ---
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
//--- StRefMult class ---
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"


ClassImp(StKFParticleAnalysisMaker);

//________________________________________________________________________________
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name) : StMaker(name), fNTrackTMVACuts(0), fIsPicoAnalysis(true), fdEdXMode(1), 
  fStoreTmvaNTuples(false), fProcessSignal(false), fCollectTrackHistograms(false), fCollectPIDHistograms(false),fTMVAselection(false), 
  fFlowAnalysis(false), fFlowChain(NULL), fFlowRunId(-1), fFlowEventId(-1), fCentrality(-1), fFlowFiles(), fFlowMap(), 
  fRunCentralityAnalysis(0), fRefmultCorrUtil(0), fCentralityFile(""), fAnalyseDsPhiPi(false), fDecays(0), fIsProduce3DEfficiencyFile(false), f3DEfficiencyFile(""), 
  fStoreCandidates(false), fPartcileCandidate(), fIsStoreCandidate(KFPartEfficiencies::nParticles, false), fCandidateFile(nullptr), fCandidatesTree(nullptr),
  fKaonAnalysis(false), fKaonFileName("kaons.root"), hKaonNtuple(NULL)
{
  memset(mBeg,0,mEnd-mBeg+1);
  
  fNTuplePDG[0] = 421;
  fNTuplePDG[1] = 411;
  fNTuplePDG[2] = 431;
  fNTuplePDG[3] = 4122;
  fNTuplePDG[4] = 426;
  fNTuplePDG[5] = 429;
  fNTuplePDG[6] = 521;
  fNTuplePDG[7] = 511;
  
  fNtupleNames[0] = "D0"; 
  fNtupleNames[1] = "DPlus"; 
  fNtupleNames[2] = "Ds"; 
  fNtupleNames[3] = "Lc";
  fNtupleNames[4] = "D0KK";
  fNtupleNames[5] = "D04";
  fNtupleNames[6] = "BPlus";
  fNtupleNames[7] = "B0";
  
  vector<TString> trackCutNames;
  trackCutNames.push_back("pt_");
  trackCutNames.push_back("chi2Primary_");
  trackCutNames.push_back("dEdXPi_");
  trackCutNames.push_back("dEdXK_");
  trackCutNames.push_back("dEdXP_");
  trackCutNames.push_back("ToFPi_");
  trackCutNames.push_back("ToFK_");
  trackCutNames.push_back("ToFP_");
  fNTrackTMVACuts = trackCutNames.size();
  
  fDaughterNames[0].push_back("K");     fDaughterNames[0].push_back("Pi");                                                                              //D0 -> Kpi
  fDaughterNames[1].push_back("K");     fDaughterNames[1].push_back("Pi1");    fDaughterNames[1].push_back("Pi2");                                      //D+ -> Kpipi
  fDaughterNames[2].push_back("KPlus"); fDaughterNames[2].push_back("KMinus"); fDaughterNames[2].push_back("Pi");                                       //Ds -> KKpi
  fDaughterNames[3].push_back("K");     fDaughterNames[3].push_back("Pi");     fDaughterNames[3].push_back("P");                                        //Lc -> pKpi
  fDaughterNames[4].push_back("KPlus"); fDaughterNames[4].push_back("KMinus");                                                                          //D0 -> KK
  fDaughterNames[5].push_back("K");     fDaughterNames[5].push_back("Pi1");    fDaughterNames[5].push_back("Pi2");  fDaughterNames[5].push_back("Pi3"); //D0 -> Kpipipi
  fDaughterNames[6].push_back("PiD");   fDaughterNames[6].push_back("KD");     fDaughterNames[6].push_back("Pi");                                       //B+ -> D0_bpi
  fDaughterNames[7].push_back("Pi1D");  fDaughterNames[7].push_back("KD");     fDaughterNames[7].push_back("Pi2D"); fDaughterNames[7].push_back("Pi");  //B0 -> D-pi+

  for(int iDecay=0; iDecay<fNNTuples; iDecay++)
  {
    for(unsigned int iDaughter=0; iDaughter<fDaughterNames[iDecay].size(); iDaughter++)
    {
      for(int iTrackTMVACut=0; iTrackTMVACut<fNTrackTMVACuts; iTrackTMVACut++)
      {
        if(iDaughter==0 && iTrackTMVACut==0)
          fNtupleCutNames[iDecay] = trackCutNames[iTrackTMVACut];  
        else
          fNtupleCutNames[iDecay] += trackCutNames[iTrackTMVACut];
        fNtupleCutNames[iDecay] += fDaughterNames[iDecay][iDaughter];
        fNtupleCutNames[iDecay] += ":";
      }
    }
    if(iDecay<6)
      fNtupleCutNames[iDecay] += "Chi2NDF:LdL:Chi2Topo:refMult";
    else if(iDecay>=6 && iDecay<8)
    {
      fNtupleCutNames[iDecay] += "Chi2NDF_D:LdL_D:Chi2Topo_D:Chi2NDF:LdL:Chi2Topo:refMult";
    } 
    
    SetTMVABins(iDecay);
  }
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
  
  if(fTMVAselection || fStoreTmvaNTuples)
  {
    for(int iReader=0; iReader<fNNTuples; iReader++)
    {
      TString cutName;
      int firstSymbolOfCutName = 0;
      
      int nCuts = 0;
      while(fNtupleCutNames[iReader].Tokenize(cutName,firstSymbolOfCutName,":"))
        nCuts++;
      fTMVAParticleParameters[iReader].resize(nCuts);
    }
  }
  
  if(fTMVAselection)
  {
    for(int iReader=0; iReader<fNNTuples; iReader++)
    {
      const int nCentralityBins = fTMVACentralityBins[iReader].size() - 1;
      const int nPtBins = fTMVAPtBins[iReader].size() - 1;
      
      for(int iCentralityBin=0; iCentralityBin<nCentralityBins; iCentralityBin++)
      {
        for(int iPtBin=0; iPtBin<nPtBins; iPtBin++)
        {
          fTMVAReader[iReader][iCentralityBin][iPtBin] = new TMVA::Reader("Silent");

          TString cutName;
          int firstSymbolOfCutName = 0;      
          unsigned int iCut = 0;
          while(fNtupleCutNames[iReader].Tokenize(cutName,firstSymbolOfCutName,":"))
          {
            fTMVAReader[iReader][iCentralityBin][iPtBin] -> AddVariable( cutName.Data(), &fTMVAParticleParameters[iReader][iCut] );
            iCut++;
            if(iCut == (fTMVAParticleParameters[iReader].size()-1)) break;
          }
          
          fTMVAReader[iReader][iCentralityBin][iPtBin] -> BookMVA("BDT", fTMVACutFile[iReader][iCentralityBin][iPtBin].Data());
        }
      }
    }
  }
      
  //Create file with NTuples for cut optimization
  if(fStoreTmvaNTuples)
  {  
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    for(int iNtuple=0; iNtuple<fNNTuples; iNtuple++)
    {
      TString SignalPrefix = "_Signal";
      if(!fProcessSignal) SignalPrefix = "_BG";
      TString currentNTupleFileName = fNtupleNames[iNtuple]+SignalPrefix+TString(".root");
      fNTupleFile[iNtuple] = new TFile(currentNTupleFileName.Data(),"RECREATE");
      fCutsNTuple[iNtuple] = new TNtuple(fNtupleNames[iNtuple].Data(), fNtupleNames[iNtuple].Data(), fNtupleCutNames[iNtuple].Data());
    }
    gFile = curFile;
    gDirectory = curDirectory;
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

  fRefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P16id();
  fRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
  fRefmultCorrUtil->readScaleForWeight("/gpfs01/star/pwg/pfederic/qVectors/StRoot/StRefMultCorr/macros/weight_grefmult_VpdnoVtx_Vpd5_Run16.txt"); //for new StRefMultCorr, Run16, SL16j
  
  //Initialise the chain with files containing centrality and reaction plane
  if(fFlowAnalysis)
  {
    std::cout << "StKFParticleAnalysisMaker: run flow analysis. Flow file list:"<<std::endl;
    
    fFlowChain = new TChain("mTree");
    for(unsigned int iFlowFile=0; iFlowFile<fFlowFiles.size(); iFlowFile++)
    {
      std::cout << "      " << fFlowFiles[iFlowFile] << std::endl;
      fFlowChain->Add(fFlowFiles[iFlowFile].Data());
    }
    
    fFlowChain->SetBranchStatus("*",0);
    fFlowChain->SetBranchAddress("runid",   &fFlowRunId);   fFlowChain->SetBranchStatus("runid", 1);
    fFlowChain->SetBranchAddress("eventid", &fFlowEventId); fFlowChain->SetBranchStatus("eventid", 1);
    fFlowChain->SetBranchAddress("cent", &fCentrality);  fFlowChain->SetBranchStatus("cent", 1);
    
    std::cout << "StKFParticleAnalysisMaker: number of entries in the flow chain" << fFlowChain->GetEntries() << std::endl;
    for(int iEntry=0; iEntry<fFlowChain->GetEntries(); iEntry++)
    {
      fFlowChain->GetEvent(iEntry);
      fFlowMap[GetUniqueEventId(fFlowRunId, fFlowEventId)] = iEntry;
    }
  }

  if(fKaonAnalysis)
  {
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    
    fKaonFile = new TFile(fKaonFileName, "RECREATE");
    hKaonNtuple = new TNtuple("kaons","kaons","RunId:eventId:Vz:mother_PID:mother_pt:mother_px:mother_py:mother_pz:mother_eta:mother_phi:mother_pt_PVX:mother_px_PVX:mother_py_PVX:mother_pz_PVX:mother_eta_PVX:mother_phi_PVX:mother_m:mother_PV_l:mother_PV_dl:\
decay_Vr:decay_Vx:decay_Vy:decay_Vz:\
daughter_index:daughter_pdg:daughter_nhits:daughter_nhitsPoss:daughter_lastHitR:daughter_nhits_dEdx:daughter_dEdx:daughter_p:daughter_pt:daughter_eta:daughter_phi:daughter_px:daughter_py:daughter_pz:dp_PVX:\
daughter_decay_dl:daughter_DecayPt:daughter_DecayEta:daughter_DecayPhi:daughter_DecayPx:daughter_DecayPy:daughter_DecayPz:dp_Decay:dp_Decay_KF:\
daughter_DecayDca_KF:daughter_DecayDca_mu:daughter_PvtxDca_KF:daughter_PvtxDca_mu:daughter_PvtxDca_official:isBest");
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
  if(!fIsPicoAnalysis && fProcessSignal) storeMCHistograms = true;
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
    for(unsigned int iTrack = 0; iTrack < fPicoDst->numberOfTracks(); iTrack++) 
    {
      StPicoTrack *gTrack = fPicoDst->track(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  else
  {
    for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfGlobalTracks(); iTrack++) 
    {
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
  
//   fStKFParticleInterface->SetTriggerMode();
//   fStKFParticleInterface->SetSoftKaonPIDMode();
//   fStKFParticleInterface->SetSoftTofPidMode();
//   fStKFParticleInterface->SetChiPrimaryCut(10);
//   
//   fStKFParticleInterface->SetPtCutCharm(0.5);
//   fStKFParticleInterface->SetChiPrimaryCutCharm(8);
//   fStKFParticleInterface->SetLdLCutCharmManybodyDecays(3);
//   fStKFParticleInterface->SetChi2TopoCutCharmManybodyDecays(10);
//   fStKFParticleInterface->SetChi2CutCharmManybodyDecays(3);
//   fStKFParticleInterface->SetLdLCutCharm2D(3);
//   fStKFParticleInterface->SetChi2TopoCutCharm2D(10);
//   fStKFParticleInterface->SetChi2CutCharm2D(3);
  
  vector<int> triggeredTracks;
  bool isGoodEvent = false;
  
  //Process the event
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
  
  if(isGoodEvent)
  {
    int centralityBin = -1;
    float centralityWeight = 0.;
    
    if(fRunCentralityAnalysis)
    {
      fRefmultCorrUtil->init(fPicoDst->event()->runId());
      if(! (fRefmultCorrUtil->isBadRun(fPicoDst->event()->runId())) )
      {
        fRefmultCorrUtil->initEvent(fPicoDst->event()->grefMult(), fPicoDst->event()->primaryVertex().z(), fPicoDst->event()->ZDCx()) ;
        centralityBin = fRefmultCorrUtil->getCentralityBin9();
        centralityWeight = fRefmultCorrUtil->getWeight();
      }
//       refmultCor = fRefmultCorrUtil->getRefMultCorr();
    }
    
    if(fTMVAselection)
    {
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
              
        for(int iReader=0; iReader<fNNTuples; iReader++)
        {
          if( abs(particle.GetPDG()) == fNTuplePDG[iReader] )
          {
            GetParticleParameters(iReader, particle);
            
            const int iTMVACentralityBin = GetTMVACentralityBin(iReader, centralityBin);
            const int iTMVAPtBin = GetTMVAPtBin(iReader, particle.GetPt());
            
            if(iTMVACentralityBin<0 || iTMVAPtBin<0) 
            {
              fStKFParticleInterface->RemoveParticle(iParticle);
              continue;
            }
            
            if(fTMVAReader[iReader][iTMVACentralityBin][iTMVAPtBin]->EvaluateMVA("BDT") < fTMVACut[iReader][iTMVACentralityBin][iTMVAPtBin])
              fStKFParticleInterface->RemoveParticle(iParticle);
            
            if(fAnalyseDsPhiPi && abs(fStKFParticleInterface->GetParticles()[iParticle].GetPDG()) == 431)
            {              
              KFParticle phi;
              if(particle.GetPDG() == 431)
                phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
              else
                phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
              phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[2]];
              float mass = 0.f, dmass = 0.f;
              phi.GetMass(mass, dmass);
              if( fabs(mass - 1.01946) > 0.015)
                fStKFParticleInterface->RemoveParticle(iParticle);
            }
          }
        }
      }      
    }

#if 1
    //clean clusters
    int nTracks = 0;

    std::vector<bool> isValidTrack(maxGBTrackIndex+1);
    for(int iTrack=0; iTrack<maxGBTrackIndex+1; iTrack++)
      isValidTrack[iTrack] = true;

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
        }
        
        if(cluster.NDaughters() > 5) {
          
//           float l, dl;
//           particle.GetDistanceToVertexLine(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex(), l, dl);
//           std::cout << "event " << iEvent << "   pdg: " << particle.GetPDG() << " l " << l <<" l/dl " << (l/dl) << "  killed by cluster" << std::endl;
        
          fStKFParticleInterface->RemoveParticle(iParticle);
          for(unsigned int iTrackId=0; iTrackId<trackIds.size(); iTrackId++)
            isValidTrack[trackIds[iTrackId]] = false;
        }
      }
    }
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


//rotational background
#if 0
    const int nParticles0 = fStKFParticleInterface->GetParticles().size();
    for(int iParticle=0; iParticle<nParticles0; iParticle++) {
      KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
      if(particle.GetPDG()==3006 || 
         particle.GetPDG()==3007 || 
         particle.GetPDG()==3012 || 
         particle.GetPDG()==3013)
      {
        KFParticle pion     = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
        KFParticle fragment = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
        KFParticle proton   = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[2]];
        
        KFParticle ppi;
        ppi += pion;
        ppi += proton;
        
        fStKFParticleInterface->RemoveParticle(iParticle);
        
        float vtx[3] = { ppi.X(), ppi.Y(), ppi.Z() };
//         const int nRotate = 3;
//         for(int i=1; i<nRotate+1; i++) 
        {
          KFParticle fragment_copy = fragment;
          KFParticle ppi_copy = ppi;
          
          fragment_copy.TransportToPoint(vtx);
//           fragment_copy.Rotate(2.*TMath::Pi()/(nRotate+1)*i, particle);
          fragment_copy.Rotate(TMath::Pi(), particle);
          fragment_copy.SetId(fStKFParticleInterface->GetParticles().size());
          fStKFParticleInterface->AddParticle(fragment_copy);         
          
          ppi_copy += fragment_copy;
          ppi_copy.SetPDG(particle.GetPDG());
          ppi_copy.SetId(fStKFParticleInterface->GetParticles().size());
          
          ppi_copy.CleanDaughtersId();
          ppi_copy.AddDaughterId(pion.Id());
          ppi_copy.AddDaughterId(fragment_copy.Id());
          ppi_copy.AddDaughterId(proton.Id());
          
          fStKFParticleInterface->AddParticle(ppi_copy);          
        }
      }
    }
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
    }
#endif

    //clean H3L, H4L, Ln, Lnn
    for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
    {
      KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
//       if( abs(particle.GetPDG())==3003 || abs(particle.GetPDG())==3103 || abs(particle.GetPDG())==3004 || abs(particle.GetPDG())==3005)
      if((abs(particle.GetPDG()) > 3002) && (abs(particle.GetPDG()) < 3200))
      {        
        for(int iD=0; iD<particle.NDaughters(); iD++)
        {
          const int daughterId = particle.DaughterIds()[iD];
          const KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterId];
          if(abs(daughter.GetPDG())==211 && daughter.GetP() > 0.7)
            fStKFParticleInterface->RemoveParticle(iParticle);
          if(abs(daughter.GetPDG())!=211 && daughter.GetP() < 0.5) //TODO remove me
            fStKFParticleInterface->RemoveParticle(iParticle);
        }
        
//         float l = sqrt(particle.X()*particle.X() + particle.Y()*particle.Y() + particle.Z()*particle.Z());
//         float r = sqrt(particle.X()*particle.X() + particle.Y()*particle.Y());
//         if(r > 50)// || (r>2.5 && r<3.6) || (r>7.5&&r<8.8))
//           fStKFParticleInterface->RemoveParticle(iParticle);
      }

#if 0
      if( (abs(particle.GetPDG()) == 3006) || 
          (abs(particle.GetPDG()) == 3007) ||
          (abs(particle.GetPDG()) == 3012) ||
          (abs(particle.GetPDG()) == 3013) ||
          (abs(particle.GetPDG()) == 3014) ||
          (abs(particle.GetPDG()) == 3015) ||
          (abs(particle.GetPDG()) == 3017) ||
          (abs(particle.GetPDG()) == 3018) ||
          (abs(particle.GetPDG()) == 3020) ||
          (abs(particle.GetPDG()) == 3021) ||
          (abs(particle.GetPDG()) == 3023) ||
          (abs(particle.GetPDG()) == 3024) ||
          (abs(particle.GetPDG()) == 3026) ||
          (abs(particle.GetPDG()) == 3027) ||
          (abs(particle.GetPDG()) == 3028)
        )
      {
//         //clean gamma and clones
//         for(int iD0=0; iD0<particle.NDaughters(); iD0++)
//         {
//           for(int iD1=0; iD1<iD0; iD1++)
//           {
//             int index0 = particle.DaughterIds()[iD0];
//             int index1 = particle.DaughterIds()[iD1];
//             KFParticle d0 = fStKFParticleInterface->GetParticles()[index0];
//             KFParticle d1 = fStKFParticleInterface->GetParticles()[index1];
//             float vertex[3] = {particle.GetX(), particle.GetY(), particle.GetZ()};
//             d0.TransportToPoint(vertex);
//             d1.TransportToPoint(vertex);
//             float qtAlpha[2];
//             KFParticle::GetArmenterosPodolanski(d0, d1, qtAlpha );
//             if(qtAlpha[0] < 0.005)
//               fStKFParticleInterface->RemoveParticle(iParticle);
//           }
//         }

        if(particle.NDaughters() == 3) {
          const KFParticle d[3] = {
            fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]],
            fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]],
            fStKFParticleInterface->GetParticles()[particle.DaughterIds()[2]]
          };
          KFParticle v[3];
          int index[3][2] = { {1,2}, {0,2}, {0,1} }; 
          
          bool ok = true;
          for(int iD=0; iD<3; iD++){
            
            const KFParticle* vd[2] = {&d[index[iD][0]], &d[index[iD][1]]};
            v[iD].Construct(vd, 2);

            float q1q2 = vd[0]->Px()*vd[1]->Px() + vd[0]->Py()*vd[1]->Py() + vd[0]->Pz()*vd[1]->Pz();
            float q12  = vd[0]->Px()*vd[0]->Px() + vd[0]->Py()*vd[0]->Py() + vd[0]->Pz()*vd[0]->Pz();
            float q22  = vd[1]->Px()*vd[1]->Px() + vd[1]->Py()*vd[1]->Py() + vd[1]->Pz()*vd[1]->Pz();
            ok &= q1q2 > -q12;
            ok &= q1q2 > -q22;
            

            float p1p2 = d[iD].Px()*v[iD].Px() + d[iD].Py()*v[iD].Py() + d[iD].Pz()*v[iD].Pz();
            float p12  = d[iD].Px()*d[iD].Px() + d[iD].Py()*d[iD].Py() + d[iD].Pz()*d[iD].Pz();
            float p22  = v[iD].Px()*v[iD].Px() + v[iD].Py()*v[iD].Py() + v[iD].Pz()*v[iD].Pz();
            ok &= p1p2 > -p12;
            ok &= p1p2 > -p22;
            
            v[iD] += d[iD];
            ok &= v[iD].Chi2()/float(v[iD].NDF()) < 3;
            
            float m=0.f, dm=1e6f;
            ok &= (v[iD].GetMass(m, dm) == 0);
            
            float l=0.f, dl=1e6f;
            v[iD].GetDistanceToVertexLine(particle, l, dl);
            ok &= l/dl < 3.f;
          }
          
          if(!ok)
            fStKFParticleInterface->RemoveParticle(iParticle);
        }
      }
#endif
      
    }
    
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
    
    if(fFlowAnalysis)
    {
      if(fIsPicoAnalysis) 
      {
        runId   = fPicoDst->event()->runId();
        eventId = fPicoDst->event()->eventId();
      }
      else
      {
        runId   = fMuDst->event()->runId();
        eventId = fMuDst->event()->eventId();
      }
    
      long entryId = GetUniqueEventId(runId, eventId);
      std::map<long,int>::iterator flowMapIterator = fFlowMap.find(entryId);
      if (flowMapIterator != fFlowMap.end())
      {
        fFlowChain->GetEvent(fFlowMap[GetUniqueEventId(runId, eventId)]);
        centralityBin = fCentrality;
      }
    }
    
    centralityWeight = 1;
    
    fStKFParticlePerformanceInterface->SetMCTracks(mcTracks);
    fStKFParticlePerformanceInterface->SetMCIndexes(mcIndices);    
    fStKFParticlePerformanceInterface->SetCentralityBin(centralityBin);
    fStKFParticlePerformanceInterface->SetCentralityWeight(centralityWeight);
    Int_t nevent = 100000;
    fStKFParticlePerformanceInterface->SetPrintEffFrequency(nevent);
    fStKFParticlePerformanceInterface->PerformanceAnalysis();
    
    if(fStoreTmvaNTuples)
    {
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        KFParticle particle;
        bool isMCParticle = fStKFParticlePerformanceInterface->GetParticle(particle, iParticle);
              
        if( !( (fProcessSignal && isMCParticle) || (!fProcessSignal && !isMCParticle) ) ) continue;
                  
        for(int iNTuple=0; iNTuple<fNNTuples; iNTuple++)
        {
          if( particle.GetPDG() == fNTuplePDG[iNTuple] )
          {
            GetParticleParameters(iNTuple, particle);
            fCutsNTuple[iNTuple]->Fill(fTMVAParticleParameters[iNTuple].data());
          }
        }
      }
    }
   if(fKaonAnalysis) Fill_KaonNtuples();
  
  } //is good event
  
  return kStOK;
}

void StKFParticleAnalysisMaker::Fill_KaonNtuples() {
  cout<<"StKFParticleAnalysisMaker::Fill_KaonNtuples() "<<endl;
 typedef struct{Float_t id=0,index=0,p=0,pt=0,eta=0,phi=0, px=0,py=0,pz=0,dp_Decay=0,dp_PVX=0,nhits=0, nhits_dEdx=0,nhits_pos=0, 
  lastPointR=0,pdg=0,DecayDca_KF=0,DecayDca_mu=0,PvtxDca_KF=0,PvtxDca_official=0,PvtxDca_mu=0,dEdx=0,isBest=0;} _daughter;
  const int nPIDS=4;
  int particlesPDG[nPIDS] = {100321, 200321,-100321, -200321}; //K->3pi with found K (200321), K->3pi only (100321)
  int nDaughters[nPIDS] = {3, 4, 3, 4};

//Note There can be mutiple primary vertexes. I should take only the  one that was intarfaced to KFP
// the 3pi should actually point "somewhere" around this vertex - I do not want potential kaons from other vertexes
// and then macth the primary kaon

 const KFParticleTopoReconstructor *topoRec=fStKFParticleInterface->GetTopoReconstructor();
 // cout<<" Primary vertices:"<<endl<<"in KF: "<<topoRec->NPrimaryVertices()<< " , ID="<<topoRec->GetPrimVertex().Id()<<endl;
  //cout<<"in MuDst: "<<fMuDst->numberOfPrimaryVertices()<<endl;
  int npt=fStKFParticleInterface->GetParticles().size();
  //cout<<"fStKFParticleInterface->GetParticles().size()="<<npt<<endl;
  /*test 
  for(unsigned int iParticle=0; iParticle<npt; iParticle++) {
    cout<<"iParticle="<<iParticle<<" ID="<<fStKFParticleInterface->GetParticles()[iParticle].Id()<<endl;
  }  */

//primary vertex
  Float_t Vz=topoRec->GetPrimVertex().GetZ();
  Float_t Vx=topoRec->GetPrimVertex().GetX();
  Float_t Vy=topoRec->GetPrimVertex().GetY();
  StThreeVectorD pVtx(Vx, Vy, Vz);

  for(unsigned int iParticle=0; iParticle<npt; iParticle++) {
    //cout<<"i="<<iParticle<<" KFP_iD="<<fStKFParticleInterface->GetParticles()[iParticle].Id()
    //<<" PDG="<<fStKFParticleInterface->GetParticles()[iParticle].GetPDG()<<" pt="<<
    //fStKFParticleInterface->GetParticles()[iParticle].GetPt()<<endl;
    //cout<< fStKFParticleInterface->GetParticles()[iParticle]<<endl;
    for(int iPDG = 0; iPDG < nPIDS; iPDG++){
      if(fStKFParticleInterface->GetParticles()[iParticle].GetPDG() != particlesPDG[iPDG]) continue;
        KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];

        Float_t runId;
        Float_t eventId;
        if(fIsPicoAnalysis) 
        {
          runId   = fPicoDst->event()->runId();
          eventId = fPicoDst->event()->eventId();
        }
        else
        {
          runId   = fMuDst->event()->runId();
          eventId = fMuDst->event()->eventId();
        }
      

        Float_t mother_PID=particle.GetPDG();
        Float_t mother_m = particle.GetMass();
       
        //decay point
        Float_t decay_Vx=particle.GetX();
        Float_t decay_Vy=particle.GetY();
        Float_t decay_Vz=particle.GetZ();
        Float_t decay_Vr=particle.GetR();
        if (decay_Vr<50) continue; // skip decays out of TPC
        StThreeVectorD decayVtx(decay_Vx, decay_Vy, decay_Vz); 

         //momentum at decay point
        Float_t mother_pt = particle.GetPt();
        Float_t mother_px = particle.GetPx();
        Float_t mother_py = particle.GetPy();
        Float_t mother_pz = particle.GetPz();
        StThreeVectorD mother_p_decay(mother_px,mother_py,mother_pz);
        Float_t mother_eta = particle.GetEta();
        Float_t mother_phi = particle.GetPhi();          

        // move to primary vertex
        particle.TransportToPoint(topoRec->GetPrimVertex().Parameters());
        //momentum at primary vertex
        Float_t mother_px_PVX = particle.GetPx();
        Float_t mother_py_PVX = particle.GetPy();
        Float_t mother_pz_PVX = particle.GetPz();
        Float_t mother_pt_PVX = particle.GetPt();
        Float_t mother_eta_PVX = particle.GetEta();
        Float_t mother_phi_PVX = particle.GetPhi();
        StThreeVectorD mother_p_PVX(mother_px_PVX,mother_py_PVX,mother_pz_PVX);
      
        //TODO!!
        //distance to PV - usefull to see if it comes from the primary vertex...I coudl have seom decays of secondary kaons
        //float_v lCandidate, dlCandidate;
         //KFParticleSIMD part(particle);
        //KFParticleSIMD pvt(topoRec->GetPrimVertex());
        //daugh.GetDistanceToVertexLine(pvt, lCandidate, dlCandidate);
       
        Float_t mother_PV_l=particle.GetDistanceFromVertex(topoRec->GetPrimVertex());
        Float_t mother_PV_dl=particle.GetDeviationFromVertex(topoRec->GetPrimVertex());
     
        
        if(particle.NDaughters() != nDaughters[iPDG]) {
            cout << "Wrong number of daughters! for"<< particlesPDG[iPDG] << " expected "<< nDaughters[iPDG]<<" but found "<< particle.NDaughters()<<endl;
             continue;
        } 
/*
      cout<<"---test------"<<endl;
      cout<<" fMuDst->numberOfGlobalTracks()="<<fMuDst->numberOfGlobalTracks()<<" fMuDst->globalTracks(fMuDst->numberOfGlobalTracks()-1).Id()"
      <<fMuDst->globalTracks(fMuDst->numberOfGlobalTracks()-1)->id()<<endl;
      cout<<"fStKFParticleInterface->GetParticles().size()="<<npt<<endl<<endl;
     
      cout<<" mother particle PDG="<<particle.GetPDG()<<" with #Daughters="<<particle.NDaughters()<<endl;
      cout<<" mother particle  pt="<<particle.GetPt()<<"  R="<<particle.GetR()<<" phi="<<particle.GetPhi()<<" eta="<<particle.GetEta()<<endl;
  */    
      //cout<<" daugterIds[0]="<<particle.DaughterIds()[0]<<endl;;
      //cout<<" daugterIds[1]="<<particle.DaughterIds()[1]<<endl;;
      //cout<<" daugterIds[2]="<<particle.DaughterIds()[2]<<endl;;
      //cout<<" StKFParticleInterface->TrackIdToI().size()="<<fStKFParticleInterface->TrackIdToI().size()<<endl;
      //cout<<" TrackIdToI[particle.DaughterIds()[0]]="<<fStKFParticleInterface->TrackIdToI()[particle.DaughterIds()[0]]<<endl;

      for(int iD=0; iD<particle.NDaughters(); iD++) {
              _daughter daughter;
             // cout<<"filling...daughter index#="<<iD<<endl;               
              daughter.index=iD;
              int id_kf=particle.DaughterIds()[iD];
              KFParticle daugh=(fStKFParticleInterface->GetParticles())[id_kf];
              daughter.id=daugh.DaughterIds()[0]; //id of parent MuDST tracks
              const int iMuTrack = fStKFParticleInterface->TrackIdToI()[daughter.id];
              StMuTrack *mutrack = (StMuTrack *) fMuDst->globalTracks(iMuTrack); 
         

              daughter.pt = daugh.GetPt();
              daughter.p  = daugh.GetP();
              daughter.pdg = daugh.GetPDG();
              daughter.eta = daugh.GetEta();
              daughter.phi = daugh.GetPhi();
              daughter.px = daugh.GetPx();
              daughter.py = daugh.GetPy();
              daughter.pz = daugh.GetPz();
              daughter.nhits=mutrack->nHitsFit(); 
              daughter.nhits_dEdx=mutrack->nHitsDedx(); 
              daughter.nhits_pos=mutrack->nHitsPoss(); 
              daughter.lastPointR=mutrack->lastPoint().perp(); //from PV to last hist
              daughter.dEdx=mutrack->dEdx();
             

              //DCA's - mainly useful for mother track - iD==3
              //a) to decay Vtx: from KFP, from MuTrack
                 //from KFP
              daughter.DecayDca_KF = daugh.GetDistanceFromVertex(particle);
              KFParticle tmp=daugh;
              tmp.TransportToParticle(particle);
              StThreeVectorD daugher_p_decay_KF(tmp.Px(),tmp.Py(),tmp.Pz());
              double dp_decay_KF=(daugher_p_decay_KF-mother_p_decay).mag();
              //muDST
              StPhysicalHelixD helix = mutrack->helix();
              double pathlength = helix.pathLength(decayVtx, true ); // false- do not scan periods
              daughter.pdg=pathlength; //temporary
              daughter.DecayDca_mu=helix.distance(decayVtx);//(helix.at(pathlength)-decayVtx).mag();
              //momentum at decay vtx - from global track!
              StThreeVectorD daughter_p_decay= helix.momentumAt(pathlength,fMuDst->event()->runInfo().magneticField()*kilogauss);
              daughter.dp_Decay=(daughter_p_decay-mother_p_decay).mag();

            
              
              //B)to prim Vtx:
              //from KFP
              daughter.PvtxDca_KF = daugh.GetDistanceFromVertex(topoRec->GetPrimVertex());
              //from MuTrack 
              pathlength = helix.pathLength(pVtx, true );
              daughter.PvtxDca_mu=(helix.at(pathlength)-pVtx).mag();
              //daughter.PvtxDca_mu=helix.distance(pVtx); o
              StThreeVectorD daughter_p_PVX= helix.momentumAt(pathlength,fMuDst->event()->runInfo().magneticField()*kilogauss);
              double dp_PVX=(daughter_p_PVX-mother_p_PVX).mag();
       
              //DCA from MuDST
              daughter.PvtxDca_official=mutrack->dcaGlobal().mag();                

              //KFP probabilities
              Float_t mother_PV_l=particle.GetDistanceFromVertex(topoRec->GetPrimVertex());
              Float_t mother_PV_dl=particle.GetDeviationFromVertex(topoRec->GetPrimVertex());
              Float_t daughter_decay_dl=daugh.GetDeviationFromVertex(particle);

              Float_t toFill[] = { runId,eventId,Vz, mother_PID, mother_pt, mother_px,mother_py,mother_pz, mother_eta, mother_phi, 
                mother_pt_PVX, mother_px_PVX, mother_py_PVX, mother_pz_PVX, mother_eta_PVX, mother_phi_PVX, mother_m,  mother_PV_l, mother_PV_dl,
                decay_Vr, decay_Vx, decay_Vy, decay_Vz,
                daughter.index, daughter.pdg, daughter.nhits,daughter.nhits_pos,daughter.lastPointR,daughter.nhits_dEdx, daughter.dEdx, daughter.p,daughter.pt, daughter.eta,daughter.phi,daughter.px,daughter.py,daughter.pz,dp_PVX,
                daughter_decay_dl,daughter_p_decay.perp(),daughter_p_decay.pseudoRapidity(),daughter_p_decay.phi(),daughter_p_decay.x(),daughter_p_decay.y(),daughter_p_decay.z(),daughter.dp_Decay,dp_decay_KF,
                daughter.DecayDca_KF,daughter.DecayDca_mu,daughter.PvtxDca_KF,daughter.PvtxDca_mu,daughter.PvtxDca_official,0};
                hKaonNtuple->Fill(toFill);
             } //NDaughters loop

      // matching of tracks from MuDST to the found decay vertex
     if(fabs(fStKFParticleInterface->GetParticles()[iParticle].GetPDG())!=100321) continue;// run only decay vertices
     if(decay_Vr <50) continue; //there can be no kaon associated

      StMuTrack *mutrack;
      _daughter best_dt;
      StThreeVectorD best_dt_p_PVX;
      best_dt.dp_Decay=10e20;
      best_dt.DecayDca_mu=10e20;
      //loop over global
      for (UInt_t k = 0; k < fMuDst->numberOfGlobalTracks(); k++) {
           mutrack = (StMuTrack *) fMuDst->array(muGlobal)->UncheckedAt(k);
      //loop over primary
      //for (UInt_t k = 0; k < fMuDst->numberOfPrimaryTracks(); k++) {
      //    track = (StMuTrack *) fMuDst->array(muPrimary)->UncheckedAt(k);
           if (! mutrack) continue;
           _daughter daughter;
         
          //DCA from MuDST
          daughter.PvtxDca_official=mutrack->dcaGlobal().mag();
          if (daughter.PvtxDca_official>10 ) continue;//junk


          daughter.index=4;
          daughter.id=mutrack->id();

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

           StPhysicalHelixD helix = mutrack->helix();
           double pathlength = helix.pathLength(decayVtx, true ); // false- do not scan periods
           daughter.pdg=pathlength;// ok, that's dirty...
           daughter.DecayDca_mu=helix.distance(decayVtx);//(helix.at(pathlength)-decayVtx).mag();
            
         
           //momentum at decay vtx - from global track!
           StThreeVectorD daughter_p_decay= helix.momentumAt(pathlength,fMuDst->event()->runInfo().magneticField()*kilogauss);
           daughter.dp_Decay=(daughter_p_decay-mother_p_decay).mag();

        
          pathlength = helix.pathLength(pVtx, true );
          daughter.PvtxDca_mu=(helix.at(pathlength)-pVtx).mag();
              //daughter.PvtxDca_mu=helix.distance(pVtx); 
             //daughter.DecayDca_mu=fabs(helix.geometricSignedDistance(v));
          StThreeVectorD daughter_p_PVX= helix.momentumAt(pathlength,fMuDst->event()->runInfo().magneticField()*kilogauss);
          daughter.dp_PVX=(daughter_p_PVX-mother_p_PVX).mag();

          if (daughter.dp_Decay>2.) continue; //junk                
          if (daughter.DecayDca_mu>2.) continue; //junk

          
          if (best_dt.dp_Decay>daughter.dp_Decay){
            daughter.isBest=1;
            if (best_dt.DecayDca_mu>daughter.DecayDca_mu) daughter.isBest=2;
            _daughter tmpd=best_dt;   best_dt=daughter;  daughter=tmpd;
            StThreeVectorD tmpp=best_dt_p_PVX; best_dt_p_PVX=daughter_p_decay; daughter_p_decay=tmpp;
          }

         
          if (daughter.dp_Decay>1000.) continue; //junk from initialization best_dt.dp_Decay=10e20;
          
          Float_t toFill[] = { runId,eventId,Vz, mother_PID, mother_pt, mother_px,mother_py,mother_pz, mother_eta, mother_phi, 
            mother_pt_PVX, mother_px_PVX, mother_py_PVX, mother_pz_PVX, mother_eta_PVX, mother_phi_PVX, mother_m,  mother_PV_l, mother_PV_dl,
            decay_Vr, decay_Vx, decay_Vy, decay_Vz,
            daughter.index, daughter.pdg, daughter.nhits,daughter.nhits_pos,daughter.lastPointR,daughter.nhits_dEdx, daughter.dEdx, daughter.p,daughter.pt, 
            daughter.eta,daughter.phi,daughter.px,daughter.py,daughter.pz,daughter.dp_PVX,
            0,daughter_p_decay.perp(),daughter_p_decay.pseudoRapidity(),daughter_p_decay.phi(),daughter_p_decay.x(),daughter_p_decay.y(),daughter_p_decay.z(),daughter.dp_Decay,0,
            daughter.DecayDca_KF,daughter.DecayDca_mu,daughter.PvtxDca_KF,daughter.PvtxDca_mu,daughter.PvtxDca_official,0};
            hKaonNtuple->Fill(toFill);        
            } //loop over primary     

        //now save the best
        if (best_dt.dp_Decay<= 1000.){
         //if (best_dt.DecayDca_mu<= 20.){
           Float_t toFill[] = { runId,eventId,Vz, mother_PID, mother_pt, mother_px,mother_py,mother_pz, mother_eta, mother_phi, 
            mother_pt_PVX, mother_px_PVX, mother_py_PVX, mother_pz_PVX, mother_eta_PVX, mother_phi_PVX, mother_m,  mother_PV_l, mother_PV_dl,
            decay_Vr, decay_Vx, decay_Vy, decay_Vz,
            best_dt.index, best_dt.pdg, best_dt.nhits,best_dt.nhits_pos,best_dt.lastPointR,best_dt.nhits_dEdx, best_dt.dEdx, best_dt.p,best_dt.pt, 
            best_dt.eta,best_dt.phi,best_dt.px,best_dt.py,best_dt.pz,best_dt.dp_PVX,
            0,best_dt_p_PVX.perp(),best_dt_p_PVX.pseudoRapidity(),best_dt_p_PVX.phi(),best_dt_p_PVX.x(),best_dt_p_PVX.y(),best_dt_p_PVX.z(),best_dt.dp_Decay,0,
            best_dt.DecayDca_KF,best_dt.DecayDca_mu,best_dt.PvtxDca_KF,best_dt.PvtxDca_mu,best_dt.PvtxDca_official,best_dt.isBest};
            hKaonNtuple->Fill(toFill);        
        }

    }//nPIDs loop
       
  }//FKparticles loop
}

void StKFParticleAnalysisMaker::GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle)
{
  if(particle.NDaughters() == 1)
  {
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts]   = particle.GetPt();
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+1] = particle.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    int trackId = particle.DaughterIds()[0];
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+2]   = fStKFParticleInterface->GetdEdXNSigmaPion(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+3]   = fStKFParticleInterface->GetdEdXNSigmaKaon(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+4]   = fStKFParticleInterface->GetdEdXNSigmaProton(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+5]   = fStKFParticleInterface->GetTofNSigmaPion(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+6]   = fStKFParticleInterface->GetTofNSigmaKaon(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+7]   = fStKFParticleInterface->GetTofNSigmaProton(trackId);
    
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
    
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3] = particle.Chi2()/particle.NDF();  
    
    KFParticleSIMD tempSIMDParticle(particle);
    float_v l,dl;
    KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3 + 1] = l[0]/dl[0];
    
    tempSIMDParticle.SetProductionVertex(pv);
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3 + 2] = 
      double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
    
    iDaughterParticle++;
  }
}

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

Int_t StKFParticleAnalysisMaker::Finish() 
{
  if(fStoreTmvaNTuples)
  {
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    for(int iNtuple=0; iNtuple<fNNTuples; iNtuple++)
    {
      fNTupleFile[iNtuple]->cd();
      fCutsNTuple[iNtuple]->Write();
    }
    gFile = curFile;
    gDirectory = curDirectory;
  }

  if (fKaonFile) {
    fKaonFile->cd();
    hKaonNtuple->Write();
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

int StKFParticleAnalysisMaker::GetTMVACentralityBin(int iReader, int centrality)
{
  for(unsigned int iBin=0; iBin<fTMVACentralityBins[iReader].size()-1; iBin++)
    if(centrality >= fTMVACentralityBins[iReader][iBin] && centrality < fTMVACentralityBins[iReader][iBin+1])
      return iBin;
  return -1;
}

int StKFParticleAnalysisMaker::GetTMVAPtBin(int iReader, double pt)
{
  for(unsigned int iBin=0; iBin<fTMVAPtBins[iReader].size()-1; iBin++)
    if(pt >= fTMVAPtBins[iReader][iBin] && pt < fTMVAPtBins[iReader][iBin+1])
      return iBin;
  return -1;
}

void StKFParticleAnalysisMaker::SetTMVACentralityBins(int iReader, TString bins)
{
  fTMVACentralityBins[iReader].clear();
  TString value; int firstSymbol = 0;      
  while(bins.Tokenize(value,firstSymbol,":"))
    fTMVACentralityBins[iReader].push_back(value.Atoi());
}

void StKFParticleAnalysisMaker::SetTMVAPtBins(int iReader, TString bins)
{
  fTMVAPtBins[iReader].clear();
  TString value; int firstSymbol = 0;      
  while(bins.Tokenize(value,firstSymbol,":"))
    fTMVAPtBins[iReader].push_back(value.Atof());
}

void StKFParticleAnalysisMaker::SetTMVABins(int iReader, TString centralityBins, TString ptBins)
{
  SetTMVACentralityBins(iReader, centralityBins);
  SetTMVAPtBins(iReader, ptBins);
  
  const int nCentralityBins = fTMVACentralityBins[iReader].size() - 1;
  const int nPtBins = fTMVAPtBins[iReader].size() - 1;
  
  fTMVACutFile[iReader].resize(nCentralityBins);
  fTMVACut[iReader].resize(nCentralityBins);
  fTMVAReader[iReader].resize(nCentralityBins);
  
  for(int iCentralityBin=0; iCentralityBin<nCentralityBins; iCentralityBin++)
  {
    fTMVACutFile[iReader][iCentralityBin].resize(nPtBins);
    fTMVACut[iReader][iCentralityBin].resize(nPtBins);
    fTMVAReader[iReader][iCentralityBin].resize(nPtBins);
  }
}

void StKFParticleAnalysisMaker::AddDecayToReconstructionList( int iDecay ) { fDecays.push_back(iDecay); }

void StKFParticleAnalysisMaker::AddCandidateToStore(int pdg) {
  KFPartEfficiencies parteff;
  const int particleIndex = parteff.GetParticleIndex(pdg);
  fIsStoreCandidate[particleIndex] = true;
  fStoreCandidates = true;
}
