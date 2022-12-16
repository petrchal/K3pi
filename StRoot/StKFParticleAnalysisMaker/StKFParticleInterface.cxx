#include "StKFParticleInterface.h"

#include "KFParticleTopoReconstructor.h"
#include "KFMCTrack.h"

#include "TMath.h"
#include "TArrayD.h"
#include "TH1F.h"
#include "TH2F.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoHelix.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StPicoEvent/StPicoArrays.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "StBichsel/Bichsel.h"
#include "StBichsel/StdEdxModel.h"
#include "StProbPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"

// #define USETOF

double R(double x,double y){
  return sqrt(x*x+y*y);
}

ClassImp(StKFParticleInterface);
StKFParticleInterface *StKFParticleInterface::fgStKFParticleInterface = 0;
StKFParticleInterface::StKFParticleInterface(): 
  fKFParticleTopoReconstructor(0), fParticles(0), fParticlesPdg(0), fNHftHits(0), fTrackIdToI(0),
  fCollectTrackHistograms(false), fCollectPIDHistograms(false),
  fStrictTofPID(true), fCleanKaonsWitTof(true), fAllIsKaonPID(false),fdEdXMode(1), fTriggerMode(false),
  fChiPrimaryCut(18.6), fChiPrimaryMaxCut(2e4), fIsFixedTarget(false),fCleanLowPVTrackEvents(false), fUseHFTTracksOnly(false)
{
  fKFParticleTopoReconstructor = new KFParticleTopoReconstructor();
  fgStKFParticleInterface = this;
  // set default cuts
  SetPrimaryProbCut(0.0001); // 0.01% to consider primary track as a secondary;
}

StKFParticleInterface::~StKFParticleInterface()
{  
  if(fKFParticleTopoReconstructor) delete fKFParticleTopoReconstructor;
  fgStKFParticleInterface = 0;
}

void StKFParticleInterface::SetField(float field) 
{ 
  if(fKFParticleTopoReconstructor)
    fKFParticleTopoReconstructor->SetField(field); 
}
void StKFParticleInterface::SetBeamLine(KFParticle& p)
{ 
  if(fKFParticleTopoReconstructor)
    fKFParticleTopoReconstructor->SetBeamLine(p);
}

void StKFParticleInterface::InitParticles()
{ 
  fKFParticleTopoReconstructor->Init( fParticles, &fParticlesPdg, &fNHftHits );
  Int_t NPV =  fKFParticleTopoReconstructor->NPrimaryVertices();
  fKFParticleTopoReconstructor->GetKFParticleFinder()->Init(NPV);
  fKFParticleTopoReconstructor->FillPVIndices();
}

#ifdef __kfpAtFirstHit__
void StKFParticleInterface::InitParticlesAtFirstAndLastHit()
{ 

  //cout<<"StKFParticleInterface::InitParticlesAtFirstAndLastHit()"<<endl;
  KFPTrackVector ftracksAtFirstHit;
  KFPTrackVector ftracksAtLastHit;
  int nTracks = fParticlesAtLastHit.size();
  //cout<<"  fParticles.size()="<<fParticles.size()<<endl;
  //cout<<"  fParticlesAtLastHit.size()="<<nTracks<<endl;
  ftracksAtFirstHit.Resize(nTracks);
  ftracksAtLastHit.Resize(nTracks);
  for(int iTr=0; iTr<nTracks; iTr++)
  {  
    int trackPDG = fParticlesPdg[iTr];
    int npixelhits = fNHftHits[iTr];
    
    for(int iP=0; iP<6; iP++) 
    {
      ftracksAtFirstHit.SetParameter(fParticles[iTr].Parameters()[iP], iP, iTr);
      ftracksAtLastHit.SetParameter(fParticlesAtLastHit[iTr].Parameters()[iP], iP, iTr); 
    }
    for(int iC=0; iC<21; iC++) 
    {
      ftracksAtFirstHit.SetCovariance(fParticles[iTr].CovarianceMatrix()[iC], iC, iTr);
      ftracksAtLastHit.SetCovariance(fParticlesAtLastHit[iTr].CovarianceMatrix()[iC], iC, iTr); 
    }
    ftracksAtFirstHit.SetId(fParticles[iTr].Id(), iTr);
    ftracksAtFirstHit.SetPDG(trackPDG, iTr);
    ftracksAtFirstHit.SetQ(fParticles[iTr].Q(), iTr);
    ftracksAtFirstHit.SetPVIndex(-1, iTr);
    ftracksAtFirstHit.SetNPixelHits(npixelhits,iTr);
    
    ftracksAtLastHit.SetId(fParticlesAtLastHit[iTr].Id(), iTr);
    ftracksAtLastHit.SetPDG(trackPDG, iTr);
    ftracksAtLastHit.SetQ(fParticlesAtLastHit[iTr].Q(), iTr);
    ftracksAtLastHit.SetPVIndex(-1, iTr);
    ftracksAtLastHit.SetNPixelHits(npixelhits,iTr);
    //cout<<" i="<<iTr<<" fParticles[i].Id()="<<fParticles[iTr].Id()<<" pt="<<fParticles[iTr].GetPt()<<endl;
    //cout<<"  at first (fParticles[i])         : pt="<<fParticles[iTr].GetPt()<<" R="<<R(fParticles[iTr].GetX(),fParticles[iTr].GetY())<<" phi="<<fParticles[iTr].GetPhi()<<";   "<<fParticles[iTr]<<endl;
    //cout<<"  at last: (fParticlesAtLastHit[i]): pt="<<fParticlesAtLastHit[iTr].GetPt()<<" R="<<R(fParticlesAtLastHit[iTr].GetX(),fParticlesAtLastHit[iTr].GetY())<<" phi="<<fParticlesAtLastHit[iTr].GetPhi()<<";   "<<fParticlesAtLastHit[iTr]<<endl;
  }
  fKFParticleTopoReconstructor->Init( ftracksAtFirstHit, ftracksAtLastHit );
  //cout<<"------------------po topo init-------"<<endl;
 
  Int_t NPV =  fKFParticleTopoReconstructor->NPrimaryVertices();
  /*
  for(unsigned int iParticle=0; GetParticles().size(); iParticle++) {
    KFParticle part=GetParticles()[iParticle];
    cout<<"i="<<iParticle<<" KFP_iD="<<part.Id()<<" Mu_iD="<<part.DaughterIds()[0]
    <<" PDG="<<part.GetPDG()<<" pt="<< part.GetPt()<<" R="<<R(part.GetX(),part.GetY())<<" phi="<< part.GetPhi()<<endl;
   }
  */
  fKFParticleTopoReconstructor->GetKFParticleFinder()->Init(NPV);
  fKFParticleTopoReconstructor->FillPVIndices();
}
#endif

void StKFParticleInterface::ReconstructParticles()
{ 
  fKFParticleTopoReconstructor->SortTracks();
  fKFParticleTopoReconstructor->ReconstructParticles();
  
//   static int iEvent=0;
//   iEvent++;
//   std::cout << "Event " << iEvent << ": init " << fKFParticleTopoReconstructor->StatTime( 0 ) 
//             << " pv " << fKFParticleTopoReconstructor->StatTime( 1 )
//             << " sort " << fKFParticleTopoReconstructor->StatTime( 2 )
//             << " particles " << fKFParticleTopoReconstructor->StatTime( 3 ) <<  std::endl;
}

void StKFParticleInterface::ReconstructTopology()
{ 
  fKFParticleTopoReconstructor->Init( fParticles, &fParticlesPdg );
  fKFParticleTopoReconstructor->ReconstructPrimVertex(0);
  fKFParticleTopoReconstructor->SortTracks();
  fKFParticleTopoReconstructor->ReconstructParticles();
  
//   static int iEvent=0;
//   iEvent++;
//   std::cout << "Event " << iEvent << ": init " << fKFParticleTopoReconstructor->StatTime( 0 ) 
//             << " pv " << fKFParticleTopoReconstructor->StatTime( 1 )
//             << " sort " << fKFParticleTopoReconstructor->StatTime( 2 )
//             << " particles " << fKFParticleTopoReconstructor->StatTime( 3 ) <<  std::endl;
}

void StKFParticleInterface::AddPV(const KFVertex &pv, const vector<int> &tracks) { 
  fKFParticleTopoReconstructor->AddPV(pv, tracks);
  fKFParticleTopoReconstructor->FillPVIndices();
}
void StKFParticleInterface::CleanPV() {
  fKFParticleTopoReconstructor->CleanPV();
}

void StKFParticleInterface::AddPV(const KFVertex &pv) { 
  fKFParticleTopoReconstructor->AddPV(pv);
}

void StKFParticleInterface::AddParticle(const KFParticle &p) { 
  fKFParticleTopoReconstructor->AddParticle(p);
}

void StKFParticleInterface::AddCandidate(const KFParticle& candidate, int iPV) {
  fKFParticleTopoReconstructor->AddCandidate(candidate, iPV);
}
void StKFParticleInterface::AddDecayToReconstructionList(Int_t pdg) {
   fKFParticleTopoReconstructor->GetKFParticleFinder()->AddDecayToReconstructionList(pdg);
}
std::vector<KFParticle> const &StKFParticleInterface::GetParticles() const { return fKFParticleTopoReconstructor->GetParticles(); }
void StKFParticleInterface::RemoveParticle(const int iParticle) { fKFParticleTopoReconstructor->RemoveParticle(iParticle); }
const std::vector<KFParticle>* StKFParticleInterface::GetSecondaryCandidates() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryCandidates();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryK0() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryK0();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryLambda() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryLambda();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryAntiLambda() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryAntiLambda();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryGamma() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryGamma();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryPi0() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryPi0();                           } // Get secondary particles with the mass constraint
const std::vector< std::vector<KFParticle> >* StKFParticleInterface::GetPrimaryCandidates() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetPrimaryCandidates();                } // Get primary particles with the mass constraint
const std::vector< std::vector<KFParticle> >* StKFParticleInterface::GetPrimaryTopoCandidates() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetPrimaryTopoCandidates();        } // Get primary particles with the topologigal constraint
const std::vector< std::vector<KFParticle> >* StKFParticleInterface::GetPrimaryTopoMassCandidates() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetPrimaryTopoMassCandidates();} // Get primary particles with the topologigal and mass constraint

KFParticleFinder*  StKFParticleInterface::GetKFParticleFinder() { return fKFParticleTopoReconstructor->GetKFParticleFinder(); }
void StKFParticleInterface::SetMaxDistanceBetweenParticlesCut(float cut) { GetKFParticleFinder()->SetMaxDistanceBetweenParticlesCut(cut); }
void StKFParticleInterface::SetLCut(float cut)                           { GetKFParticleFinder()->SetLCut(cut); }
void StKFParticleInterface::SetChiPrimaryCut2D(float cut)                { GetKFParticleFinder()->SetChiPrimaryCut2D(cut); }
void StKFParticleInterface::SetChi2Cut2D(float cut)                      { GetKFParticleFinder()->SetChi2Cut2D(cut); }
void StKFParticleInterface::SetLdLCut2D(float cut)                       { GetKFParticleFinder()->SetLdLCut2D(cut); }
void StKFParticleInterface::SetLdLCutXiOmega(float cut)                  { GetKFParticleFinder()->SetLdLCutXiOmega(cut); }
void StKFParticleInterface::SetChi2TopoCutXiOmega(float cut)             { GetKFParticleFinder()->SetChi2TopoCutXiOmega(cut); }
void StKFParticleInterface::SetChi2CutXiOmega(float cut)                 { GetKFParticleFinder()->SetChi2CutXiOmega(cut); }
void StKFParticleInterface::SetChi2TopoCutResonances(float cut)          { GetKFParticleFinder()->SetChi2TopoCutResonances(cut); }
void StKFParticleInterface::SetChi2CutResonances(float cut)              { GetKFParticleFinder()->SetChi2CutResonances(cut); }
void StKFParticleInterface::SetPtCutLMVM(float cut)                      { GetKFParticleFinder()->SetPtCutLMVM(cut); }
void StKFParticleInterface::SetPCutLMVM(float cut)                       { GetKFParticleFinder()->SetPCutLMVM(cut); }
void StKFParticleInterface::SetPtCutJPsi(float cut)                      { GetKFParticleFinder()->SetPtCutJPsi(cut); }
void StKFParticleInterface::SetPtCutCharm(float cut)                     { GetKFParticleFinder()->SetPtCutCharm(cut); }
void StKFParticleInterface::SetChiPrimaryCutCharm(float cut)             { GetKFParticleFinder()->SetChiPrimaryCutCharm(cut); }
void StKFParticleInterface::SetLdLCutCharmManybodyDecays(float cut)      { GetKFParticleFinder()->SetLdLCutCharmManybodyDecays(cut); }
void StKFParticleInterface::SetChi2TopoCutCharmManybodyDecays(float cut) { GetKFParticleFinder()->SetChi2TopoCutCharmManybodyDecays(cut); }
void StKFParticleInterface::SetChi2CutCharmManybodyDecays(float cut)     { GetKFParticleFinder()->SetChi2CutCharmManybodyDecays(cut); }
void StKFParticleInterface::SetLdLCutCharm2D(float cut)                  { GetKFParticleFinder()->SetLdLCutCharm2D(cut); }
void StKFParticleInterface::SetChi2TopoCutCharm2D(float cut)             { GetKFParticleFinder()->SetChi2TopoCutCharm2D(cut); }
void StKFParticleInterface::SetChi2CutCharm2D(float cut)                 { GetKFParticleFinder()->SetChi2CutCharm2D(cut); }
  
double StKFParticleInterface::InversedChi2Prob(double p, int ndf) const
{
  double epsilon = 1.e-14;
  double chi2Left = 0.f;
  double chi2Right = 10000.f;
  
  double probLeft = p - TMath::Prob(chi2Left, ndf);
  
  double chi2Centr = (chi2Left+chi2Right)/2.f;
  double probCentr = p - TMath::Prob( chi2Centr, ndf);
  
  while( TMath::Abs(chi2Right-chi2Centr)/chi2Centr > epsilon )
  {
    if(probCentr * probLeft > 0.f)
    {
      chi2Left = chi2Centr;
      probLeft = probCentr;
    }
    else
    {
      chi2Right = chi2Centr;
    }
    
    chi2Centr = (chi2Left+chi2Right)/2.f;
    probCentr = p - TMath::Prob( chi2Centr, ndf);
  }
  
  return chi2Centr;
}

void StKFParticleInterface::SetPrimaryProbCut(float prob)
{ 
  fKFParticleTopoReconstructor->SetChi2PrimaryCut( InversedChi2Prob(prob, 2) );
}

void StKFParticleInterface::CollectTrackHistograms()
{
  TDirectory *dirs[7] = {0};
  dirs[0] = TDirectory::CurrentDirectory(); assert(dirs[0]);
  dirs[0]->cd();
  if (! dirs[0]->GetDirectory("Tracks")) {
    dirs[0]->mkdir("Tracks");
  }
  dirs[1] = dirs[0]->GetDirectory("Tracks"); assert(dirs[1]);
  dirs[1]->cd();
  
  fTrackHistograms2D[0] = (TH2F *)   dirs[1]->Get("hdEdX");
  if (! fTrackHistograms2D[0]) fTrackHistograms2D[0] = new TH2F("hdEdX", "hdEdX", 1000, 0, 10, 200, 0, 200);

  fTrackHistograms2D[1] = (TH2F *)   dirs[1]->Get("hdEdXPos");
  if (! fTrackHistograms2D[1]) fTrackHistograms2D[1] = new TH2F("hdEdXPos", "hdEdXPos", 1000, 0, 10, 200, 0, 200);
  
  fTrackHistograms2D[2] = (TH2F *)   dirs[1]->Get("hdEdXNeg");
  if (! fTrackHistograms2D[2]) fTrackHistograms2D[2] = new TH2F("hdEdXNeg", "hdEdXNeg", 1000, 0, 10, 200, 0, 200);
  
  fTrackHistograms2D[3] = (TH2F *)   dirs[1]->Get("hdEdXwithToF");
  if (! fTrackHistograms2D[3]) fTrackHistograms2D[3] = new TH2F("hdEdXwithToF", "hdEdXwithToF", 1000, 0, 10, 200, 0, 200);
  
  fTrackHistograms2D[4] = (TH2F *)   dirs[1]->Get("hTofPID");
  if (! fTrackHistograms2D[4]) fTrackHistograms2D[4] = new TH2F("hTofPID", "hTofPID", 300, 0, 15, 1100, -1, 10);

  fTrackHistograms[0] = (TH1F *)   dirs[1]->Get("hNHFTHits");
  if (! fTrackHistograms[0]) fTrackHistograms[0] = new TH1F("hNHFTHits", "hNHFTHits",11, -0.5, 10.5);
  
  fTrackHistograms[1] = (TH1F *)   dirs[1]->Get("hPVError");
  if (! fTrackHistograms[1]) fTrackHistograms[1] = new TH1F("hPVError", "hPVError", 10000, 0, 1);

  fTrackHistograms2D[5] = (TH2F *)   dirs[1]->Get("hPVErrorVsNTracks");
  if (! fTrackHistograms2D[5]) fTrackHistograms2D[5] = new TH2F("hPVErrorVsNTracks", "hPVErrorVsNTracks", 5000, 0, 5000, 5000, 0, 0.5);

  fTrackHistograms2D[6] = (TH2F *)   dirs[1]->Get("hPVErrorVsNPVTracks");
  if (! fTrackHistograms2D[6]) fTrackHistograms2D[6] = new TH2F("hPVErrorVsNPVTracks", "hPVErrorVsNPVTracks", 5000, 0, 5000, 5000, 0, 0.5);

  fTrackHistograms[2] = (TH1F *)   dirs[1]->Get("hPrimaryRatio");
  if (! fTrackHistograms[2]) fTrackHistograms[2] = new TH1F("hPrimaryRatio", "hPrimaryRatio", 100, 0, 1);
  
  fTrackHistograms2D[7] = (TH2F *)   dirs[1]->Get("hSecondaryVsPrimaryTracks");
  if (! fTrackHistograms2D[7]) fTrackHistograms2D[7] = new TH2F("hSecondaryVsPrimaryTracks", "hSecondaryVsPrimaryTracks", 1000, 0, 1000, 1000, 0, 1000);

  dirs[0]->cd();
  
  fCollectTrackHistograms = true;
}

void StKFParticleInterface::CollectPIDHistograms()
{
  TDirectory *dirs[7] = {0};
  dirs[0] = TDirectory::CurrentDirectory(); assert(dirs[0]);
  dirs[0]->cd();
  if (! dirs[0]->GetDirectory("Tracks")) {
    dirs[0]->mkdir("Tracks");
  }
  dirs[1] = dirs[0]->GetDirectory("Tracks"); assert(dirs[1]);
  dirs[1]->cd();
  
  int pdgTrackHisto[NTrackHistoFolders] = { 11, -11, 13, -13, 211, -211, 321, -321, 2212, -2212, 
                                            1000010020, -1000010020, 1000010030, -1000010030, 
                                            1000020030, -1000020030, 1000020040, -1000020040,
                                            1000020060, -1000020060, 1000030060, -1000030060,
                                            1000030070, -1000030070, 1000040070, -1000040070 };
  TString trackFolderName[NTrackHistoFolders] = {"e-", "e+", "mu-", "mu+", "pi+", "pi-", "K+", "K-", "p", "p-", 
                                                 "d", "d-", "t", "t-", "He3", "He3-", "He4", "He4-", 
                                                 "He6", "He6-", "Li6", "Li6-", "Li7", "Li7-", "Be7", "Be7-"};
                    
  for(int iTrackHisto=0; iTrackHisto<NTrackHistoFolders; iTrackHisto++)
  {
    if (!dirs[1]->GetDirectory(trackFolderName[iTrackHisto].Data()))
      dirs[1]->mkdir(trackFolderName[iTrackHisto].Data());
    
    dirs[2] = dirs[1]->GetDirectory(trackFolderName[iTrackHisto].Data()); assert(dirs[2]);
    dirs[2]->cd();
    
    fTrackPdgToHistoIndex[ pdgTrackHisto[iTrackHisto] ] = iTrackHisto;
    
    fHistodEdXTracks[iTrackHisto] = (TH2F *)   dirs[2]->Get("hdEdX");
    if (! fHistodEdXTracks[iTrackHisto]) fHistodEdXTracks[iTrackHisto] = new TH2F("hdEdX", "hdEdX", 1000, 0, 10, 200, 0, 200);

    fHistodEdXwithToFTracks[iTrackHisto] = (TH2F *)   dirs[2]->Get("hdEdXwithToF");
    if (! fHistodEdXwithToFTracks[iTrackHisto]) fHistodEdXwithToFTracks[iTrackHisto] = new TH2F("hdEdXwithToF", "hdEdXwithToF", 1000, 0, 10, 200, 0, 200);
  
    fHistoTofPIDTracks[iTrackHisto] = (TH2F *)   dirs[2]->Get("hTofPID");
    if (! fHistoTofPIDTracks[iTrackHisto]) fHistoTofPIDTracks[iTrackHisto] = new TH2F("hTofPID", "hTofPID", 300, 0, 15, 1100, -1, 10);
  
    fHistoMomentumTracks[iTrackHisto] = (TH1F *)   dirs[2]->Get("hMomentum");
    if (! fHistoMomentumTracks[iTrackHisto]) fHistoMomentumTracks[iTrackHisto] = new TH1F("hMomentum", "hMomentum", 1000, 0, 10);
    
    fHistodEdXPull[iTrackHisto] = (TH2F *)   dirs[2]->Get("hdEdXPull");
    if (! fHistodEdXPull[iTrackHisto]) fHistodEdXPull[iTrackHisto] = new TH2F("hdEdXPull", "hdEdXPull", 2000, 0, 10, 600, -30, 30);
    
    fHistodEdXZ[iTrackHisto] = (TH2F *)   dirs[2]->Get("hdEdXZ");
    if (! fHistodEdXZ[iTrackHisto]) fHistodEdXZ[iTrackHisto] = new TH2F("hdEdXZ", "hdEdXZ", 2000, -5, 5, 280, -1, 6);
    
    dirs[1]->cd();
  }
  
  dirs[0]->cd();
  
  fCollectPIDHistograms = true;
}

bool StKFParticleInterface::IsGoodPV(const KFVertex& pv)
{
  bool isGoodPV = (pv.X() > -0.3) && (pv.X() < -0.1) &&
                  (pv.Y() > -0.27) && (pv.Y() < -0.1);
  return isGoodPV;
}

bool StKFParticleInterface::GetTrack(const StDcaGeometry& dcaG, KFPTrack& track, int q, int index)
{
  Double_t xyzp[6], CovXyzp[21];
  dcaG.GetXYZ(xyzp,CovXyzp);
  
  bool goodTrack=1;
  for(int iPar=0; iPar<6; iPar++)
    goodTrack = goodTrack && finite(xyzp[iPar]);
  for(int iC=0; iC<21; iC++)
    goodTrack = goodTrack && finite(CovXyzp[iC]);
  goodTrack &= goodTrack && CovXyzp[0]  >=0.f && CovXyzp[0]  < 100.f;
  goodTrack &= goodTrack && CovXyzp[2]  >=0.f && CovXyzp[2]  < 100.f;
  goodTrack &= goodTrack && CovXyzp[5]  >=0.f && CovXyzp[5]  < 100.f;
  goodTrack &= goodTrack && CovXyzp[9]  >=0.f && CovXyzp[9]  < 1.f;
  goodTrack &= goodTrack && CovXyzp[14] >=0.f && CovXyzp[14] < 1.f;
  goodTrack &= goodTrack && CovXyzp[20] >=0.f && CovXyzp[20] < 1.f;
  if(!goodTrack) return false;
  
  track.SetParameters(xyzp);
  track.SetCovarianceMatrix(CovXyzp);
  track.SetNDF(1);
  //    track.SetChi2(GlobalTracks_mChiSqXY[k]);
  track.SetID(index);

  track.SetCharge(q);
  return true;
}

std::vector<int> StKFParticleInterface::GetTofPID(double m2, double p, int q, const int trackId)
{
  static const int order = 4;
  static const double parMean[6][order+1] = { { 0.02283190,-0.01482910, 0.01883130,-0.01824250, 0.00409811  }, //pi+
                                              { 0.24842500,-0.00699781,-0.00991387, 0.01327170,-0.00694824  }, //K+
                                              { 0.863211  , 0.0264171 ,-0.0230833 , 0.00239637, 0.000262309 }, //p
                                              { 0.0224095 ,-0.0123235 , 0.0145216 ,-0.0149944 , 0.00325952  }, //pi-
                                              { 0.250696  ,-0.0151308 , 0.00437457, 0.00516669,-0.00529184  }, //K-
                                              { 0.886912  ,-0.0298543 , 0.0449904 ,-0.0286879 , 0.00541963  }};//p-
  static const double parSigma[6][order+1] = { { 0.0112498,-0.0400571, 0.0733615,-0.0316505, 0.00629469 }, //pi+
                                               { 0.0154830,-0.0396312, 0.0719647,-0.0290683, 0.00637164 }, //K+
                                               { 0.114465 ,-0.287213 , 0.356536 ,-0.169257 , 0.0299844  }, //p
                                               { 0.0111682,-0.0394877, 0.0718342,-0.0302914, 0.00587317 }, //pi-
                                               { 0.0157322,-0.0402606, 0.0716639,-0.0272101, 0.00564467 }, //K-
                                               { 0.0899438,-0.211922 , 0.273122 ,-0.129597 , 0.0231844  }};//p-
  double pMax = 2.;
  double nSigmas[3];
  for(int iHypothesys = 0; iHypothesys<3; iHypothesys++)
  {
    double x = p;
    if(x>=pMax) x = pMax;
    
    int iSet = iHypothesys;
    if(q<0)
      iSet += 3;
    double mean = 0;
    for(int iTerm=0; iTerm<=order; iTerm++)
      mean += parMean[iSet][iTerm]*TMath::Power(x,iTerm);  
    
    double sigma = 0;
    for(int iTerm=0; iTerm<=order; iTerm++)
      sigma += parSigma[iSet][iTerm]*TMath::Power(x,iTerm);  
    
    nSigmas[iHypothesys] = fabs((m2 - mean)/sigma);
    fTrackPidTof[iHypothesys][trackId] = nSigmas[iHypothesys];
  }
  
  double minNSigma = nSigmas[0];
  int minHypothesis = 0;
  for(int iHypothesys=1; iHypothesys<3; iHypothesys++)
  {
    if(minNSigma > nSigmas[iHypothesys]) 
    {
      minNSigma = nSigmas[iHypothesys];
      minHypothesis = iHypothesys;
    }
  }

  int pdgHypothesis[3] = {211, 321, 2212};
  vector<int> tofPID;
  
  if(fStrictTofPID)
  {
    if(minNSigma < 3)
      tofPID.push_back(pdgHypothesis[minHypothesis]*q);
  }
  else
  {    
    for(int iHypothesys=0; iHypothesys<3; iHypothesys++)
      if(nSigmas[iHypothesys] < 3)
        tofPID.push_back(pdgHypothesis[iHypothesys]*q);
  }
  
  return tofPID;
}

std::vector<int> StKFParticleInterface::GetPID(double m2, double p, int q, double dEdX, double dEdXPull[8], bool isTofm2, const int trackId)
{
  vector<int> ToFPDG;
  if(isTofm2)
    ToFPDG = GetTofPID(m2, p, q, trackId);
  
  for(int iPdg=0; iPdg<3; iPdg++)
    fTrackPidTpc[iPdg][trackId] = dEdXPull[iPdg+1];
  
  vector<int> dEdXPDG;
  float nSigmaCut = 3.f; //TODO
  
  bool checkKTof = false;
  if(fCleanKaonsWitTof)
    checkKTof = (p > 0.5) && (p < 2.);
  bool checkKHasTof = 0;
  for(unsigned int iTofPDG=0; iTofPDG<ToFPDG.size(); iTofPDG++)
    if(abs(ToFPDG[iTofPDG]) == 321)
      checkKHasTof = 1;

  if(dEdXPull[2] < nSigmaCut) dEdXPDG.push_back(211*q);  
  if (fAllIsKaonPID) dEdXPDG.push_back(321*q);
  else                                       
    if(dEdXPull[3] < 2.f && ((checkKTof && checkKHasTof) || !checkKTof) ) dEdXPDG.push_back(321*q);
  if(dEdXPull[4] < nSigmaCut)                                           dEdXPDG.push_back(2212*q); 
      
  vector<int> totalPDG;
  if(!isTofm2)
    totalPDG = dEdXPDG;
  else
  {
    for(unsigned int iPDG=0; iPDG<dEdXPDG.size(); iPDG++)
      for(unsigned int iTofPDG=0; iTofPDG<ToFPDG.size(); iTofPDG++)
        if(dEdXPDG[iPDG] == ToFPDG[iTofPDG])
          totalPDG.push_back(ToFPDG[iTofPDG]);        
  }

  if(dEdXPull[0] < nSigmaCut)  totalPDG.push_back(-11*q);
  if(dEdXPull[1] < nSigmaCut)  totalPDG.push_back(-13*q);
  if(dEdXPull[9] < nSigmaCut  || 
     dEdXPull[10] < nSigmaCut ||
     dEdXPull[11] < nSigmaCut)  totalPDG.push_back(2000003112*q);
  
  {
    
    //d
#ifdef USETOF
    if(isTofm2 && (m2>3 && m2<4.2))
#endif
    {
      if(p<1.5){
        double lowerParameters[4] = {7.11737e+00,-1.31428e+00, 1.96720e-01, 6.47905e-02};
        double lowerDBound = lowerParameters[0]*TMath::Power(p, lowerParameters[1] + lowerParameters[2]*log(p) + lowerParameters[3]*log(p)*log(p));
   
        double upperParameters[4] = {1.39824e+01,-1.53919e+00, 1.10899e-01, 9.82910e-02};
        double upperDBound = upperParameters[0]*TMath::Power(p, upperParameters[1] + upperParameters[2]*log(p) + upperParameters[3]*log(p)*log(p));
        
        if(dEdX > lowerDBound && dEdX < upperDBound)
          totalPDG.push_back(1000010020*q); 
      }
      else if(dEdXPull[5] < nSigmaCut && dEdX < 8.) 
        totalPDG.push_back(1000010020*q); 
    }
    
    //t
#ifdef USETOF
    if(isTofm2 && (m2>6.8 && m2<9.1))
#endif
    {
      if(p<2.5) {
        double lowerParameters[4] = {1.38117e+01,-1.67910e+00,-4.52185e-03, 9.21224e-02};
        double lowerTBound = lowerParameters[0]*TMath::Power(p, lowerParameters[1] + lowerParameters[2]*log(p) + lowerParameters[3]*log(p)*log(p));

        double upperParameters[4] = {2.29456e+01,-1.41456e+00, 1.04286e-01, 1.26818e-01};
        double upperTBound = upperParameters[0]*TMath::Power(p, upperParameters[1] + upperParameters[2]*log(p) + upperParameters[3]*log(p)*log(p));
        
        if(dEdX > lowerTBound && dEdX < upperTBound)
          totalPDG.push_back(1000010030*q); 
      }
      else if(dEdXPull[6] < nSigmaCut && dEdX < 8.) 
        totalPDG.push_back(1000010030*q);
    }
    
    //He3   
    if(p<3.0)
    {
      double lowerParameters[4] = {2.27715e+01,-1.36600e+00, 3.01143e-01, 2.38046e-01};
      double lowerHe3Bound = lowerParameters[0]*TMath::Power(p, lowerParameters[1] + lowerParameters[2]*log(p) + lowerParameters[3]*log(p)*log(p));
      
      double upperParameters[4] = {3.33751e+01,-1.22800e+00, 2.98371e-01, 1.82920e-01};
      double upperHe3Bound = upperParameters[0]*TMath::Power(p, upperParameters[1] + upperParameters[2]*log(p) + upperParameters[3]*log(p)*log(p));
      
      if(dEdX > lowerHe3Bound && dEdX < upperHe3Bound) 
      {
//         if(p<1.0)
//         {
//           if( isTofm2 && (m2>1.) && (m2<3.) )
//             totalPDG.push_back(1000020030*q);
//         }
//         else
//         {
          if( !isTofm2 || (isTofm2 && (m2>1.) && (m2<3.) ) )
            totalPDG.push_back(1000020030*q);
//         }
      }
    }
    else if(p>=3.0 && dEdX > 11. && dEdX < 18.)
    {
      if(dEdXPull[7] < nSigmaCut) 
        if( !isTofm2 || (isTofm2 && (m2>1.) && (m2<3.) ) )
          totalPDG.push_back(1000020030*q);
    }
    //He4
    if(p<4.0)
    {
      double lowerParameters[4] = {3.45107e+01, -1.22371e+00,  1.89140e-01,  1.07000e-01};
      double lowerHe4Bound = lowerParameters[0]*TMath::Power(p, lowerParameters[1] + lowerParameters[2]*log(p) + lowerParameters[3]*log(p)*log(p));
      
      double upperParameters[4] = {5.09065e+01, -1.36283e+00,  1.90657e-01,  1.98235e-01};
      double upperHe4Bound = upperParameters[0]*TMath::Power(p, upperParameters[1] + upperParameters[2]*log(p) + upperParameters[3]*log(p)*log(p));
      
      if(dEdX > lowerHe4Bound && dEdX < upperHe4Bound) 
      {
//         if(p<1.0)
//         {
//           if( isTofm2 && (m2>3) && (m2<4.2) )
//             totalPDG.push_back(1000020040*q);
//         }
//         else
//         {
          if( !isTofm2 || (isTofm2 && (m2>3) && (m2<4.2) ) )
            totalPDG.push_back(1000020040*q);
//         }
      }
    }
    else if(p>=4.0 && dEdX > 11. && dEdX < 18.)
    {
      if(dEdXPull[8] < nSigmaCut) 
        if( !isTofm2 || (isTofm2 && (m2>3) && (m2<4.2) ) )
          totalPDG.push_back(1000020040*q);
    } 
// }
  }
  
  //He6
  if(p<4.5) {
    double loweParameters[4] = {5.63562e+01,-1.29479e+00, 2.27883e-01,-4.10513e-02};
    double lowerHe6Bound = loweParameters[0]*TMath::Power(p, loweParameters[1] + loweParameters[2]*log(p) + loweParameters[3]*log(p)*log(p));

    double upperParameters[4] = {7.52605e+01,-1.42948e+00, 5.87043e-01,-2.13013e-01};
    double upperHe6Bound = upperParameters[0]*TMath::Power(p, upperParameters[1] + upperParameters[2]*log(p) + upperParameters[3]*log(p)*log(p));
    
    if(dEdX > lowerHe6Bound && dEdX < upperHe6Bound) {
      if(p>3.0) {
        if(isTofm2 && m2>6.5)
          totalPDG.push_back(1000020060*q);
      }
      else{
        if(!isTofm2 || (isTofm2 && m2>6.5))
          totalPDG.push_back(1000020060*q);
      }
    }
  }
  else if(p>=4.5 && dEdX > 11. && dEdX < 18.)
    if(isTofm2 && m2>6.5)
      totalPDG.push_back(1000020060*q);

  //Li6
  if(p<4){
    double loweParameters[4] = {7.30295e+01,-1.08787e+00, 6.87593e-02, 1.14228e-01};
    double lowerLi6Bound = loweParameters[0]*TMath::Power(p, loweParameters[1] + loweParameters[2]*log(p) + loweParameters[3]*log(p)*log(p));

    double upperParameters[4] = {9.35347e+01,-1.25594e+00, 2.91456e-01, 9.52847e-02};
    double upperLi6Bound = upperParameters[0]*TMath::Power(p, upperParameters[1] + upperParameters[2]*log(p) + upperParameters[3]*log(p)*log(p));
    
    if(dEdX > lowerLi6Bound && dEdX < upperLi6Bound){
      if(p<1.){
        if( isTofm2 && (m2>2.8) && (m2<4.2) )
          totalPDG.push_back(1000030060*q);
      }
      else{
        if( !isTofm2 || (isTofm2 && (m2>2.8) && (m2<4.2) ) )
          totalPDG.push_back(1000030060*q);
      }
    }
  }
  else if(p>=4.0 && dEdX > 25. && dEdX < 37.)
    if( !isTofm2 || (isTofm2 && (m2>2.8) && (m2<4.2)) )
      totalPDG.push_back(1000030060*q);

  //Li7
  if(p<4){
    double loweParameters[4] = {9.30989e+01,-1.22084e+00, 3.73173e-01,-1.12695e-01};
    double lowerLi6Bound = loweParameters[0]*TMath::Power(p, loweParameters[1] + loweParameters[2]*log(p) + loweParameters[3]*log(p)*log(p));

    double upperParameters[4] = {1.14003e+02,-1.33179e+00, 4.19395e-01,-3.20841e-02};
    double upperLi6Bound = upperParameters[0]*TMath::Power(p, upperParameters[1] + upperParameters[2]*log(p) + upperParameters[3]*log(p)*log(p));
    
    if(dEdX > lowerLi6Bound && dEdX < upperLi6Bound){
      if(p<1.){
        if( isTofm2 && (m2>4) && (m2<6) )
          totalPDG.push_back(1000030070*q);
      }
      else{
        if( !isTofm2 || (isTofm2 && (m2>4) && (m2<6) ) )
          totalPDG.push_back(1000030070*q);
      }
    }
  }
  else if(p>=4.0 && dEdX > 25. && dEdX < 37.)
    if( !isTofm2 || (isTofm2 && (m2>4) && (m2<6)) )
      totalPDG.push_back(1000030070*q);
    
  //Be7
  if(p<4){
    double loweParameters[4] = { 1.08163e+02,-1.08057e+00, 2.34159e-01, 1.98949e-02};
    double lowerLi6Bound = loweParameters[0]*TMath::Power(p, loweParameters[1] + loweParameters[2]*log(p) + loweParameters[3]*log(p)*log(p));

    double upperParameters[4] = {1.37012e+02,-1.14016e+00, 3.73116e-01,-1.85678e-02};
    double upperLi6Bound = upperParameters[0]*TMath::Power(p, upperParameters[1] + upperParameters[2]*log(p) + upperParameters[3]*log(p)*log(p));
    
    if(dEdX > lowerLi6Bound && dEdX < upperLi6Bound){
      if(p<1.){
        if( isTofm2 && (m2>2) && (m2<4) )
          totalPDG.push_back(1000040070*q);
      }
      else{
        if( !isTofm2 || (isTofm2 && (m2>2) && (m2<4) ) )
          totalPDG.push_back(1000040070*q);
      }
    }
  }
  else if(p>=4.0 && dEdX > 40. && dEdX < 55.)
    if( !isTofm2 || (isTofm2 && (m2>2) && (m2<4)) )
      totalPDG.push_back(1000040070*q);
    
    
    
  if(totalPDG.size() == 0)
    totalPDG.push_back(-1);
  
  return totalPDG;
}

void StKFParticleInterface::AddTrackToParticleList(const KFPTrack& track, int nHftHitsInTrack, int index, const std::vector<int>& totalPDG, KFVertex& pv, 
  std::vector<int>& primaryTrackList, std::vector<int>& nHftHits, std::vector<int>& particlesPdg, std::vector<KFParticle>& particles, int& nPartSaved,
  const KFPTrack* trackAtLastHit, std::vector<KFParticle>* particlesAtLastHit)
{

   for(unsigned int iPDG=0; iPDG<totalPDG.size(); iPDG++)
  {
    if( fTriggerMode && (nHftHitsInTrack < 3) ) continue;
    int pdg = totalPDG[iPDG];
    
    KFPTrack trackPDG = track;

    float r = sqrt(trackPDG.GetX()*trackPDG.GetX() + trackPDG.GetY()*trackPDG.GetY());
    float dx = (trackPDG.GetX() - pv.X());
    float dy = (trackPDG.GetY() - pv.Y());
    float dz = (trackPDG.GetZ() - pv.Z());
    float l = sqrt(dx*dx + dy*dy + dz*dz);
    //TODO uncomment!!!
//     if((abs(pdg) > 1000000000) && (!(r<5. && l<10.))) continue;
//     if((abs(pdg) > 1000) && (!(r<5. && l<10.))) continue;

    //correct for the charge of ions
    const int index2[9] = { 6,7,8, 10,11,12, 15,16,17 };
    const int index4[6] = { 9, 13,14, 18,19,20 };
    if(abs(pdg) == 1000020030 || abs(pdg) == 1000020040 || abs(pdg) == 1000020060) {
      trackPDG.SetCharge( trackPDG.Charge()*2.f );
      trackPDG.SetPx( trackPDG.GetPx()*2.f );
      trackPDG.SetPy( trackPDG.GetPy()*2.f );
      trackPDG.SetPz( trackPDG.GetPz()*2.f );
      for(int iIndex=0; iIndex<9; iIndex++){
        const int iC = index2[iIndex];
        trackPDG.SetCovariance( iC, trackPDG.GetCovariance(iC)*2.f );
      }
      for(int iIndex=0; iIndex<6; iIndex++){
        const int iC = index4[iIndex];
        trackPDG.SetCovariance( iC, trackPDG.GetCovariance(iC)*4.f );
      }
    }
    if(abs(pdg) == 1000030060 || abs(pdg) == 1000030070) {
      trackPDG.SetCharge( trackPDG.Charge()*3.f );
      trackPDG.SetPx( trackPDG.GetPx()*3.f );
      trackPDG.SetPy( trackPDG.GetPy()*3.f );
      trackPDG.SetPz( trackPDG.GetPz()*3.f );
      for(int iIndex=0; iIndex<9; iIndex++) {
        const int iC = index2[iIndex];
        trackPDG.SetCovariance( iC, trackPDG.GetCovariance(iC)*3.f );
      }
      for(int iIndex=0; iIndex<6; iIndex++) {
        const int iC = index4[iIndex];
        trackPDG.SetCovariance( iC, trackPDG.GetCovariance(iC)*9.f );
      }
    }
    if(abs(pdg) == 1000040070) {
      trackPDG.SetCharge( trackPDG.Charge()*4.f );
      trackPDG.SetPx( trackPDG.GetPx()*4.f );
      trackPDG.SetPy( trackPDG.GetPy()*4.f );
      trackPDG.SetPz( trackPDG.GetPz()*4.f );
      for(int iIndex=0; iIndex<9; iIndex++) {
        const int iC = index2[iIndex];
        trackPDG.SetCovariance( iC, trackPDG.GetCovariance(iC)*4.f );
      }
      for(int iIndex=0; iIndex<6; iIndex++) {
        const int iC = index4[iIndex];
        trackPDG.SetCovariance( iC, trackPDG.GetCovariance(iC)*16.f );
      }
    }

//     if(abs(pdg) == 211)
//     {
//       trackPDG.SetCharge( -trackPDG.Charge() );
//       pdg = -pdg;
//     }
    
    //FIXME uncomment
//     nHftHits[nPartSaved] = nHftHitsInTrack;
    nHftHits[nPartSaved] = 0;
    
    KFParticle particle(trackPDG, pdg);
    //cout<<particle<<endl;
    float chiPrim = particle.GetDeviationFromVertex(pv);
    //cout<<"  chiPrim="<<chiPrim<<"  ? "<<fChiPrimaryCut<<"  triggermode="<<fTriggerMode<<endl;
    if(chiPrim < fChiPrimaryCut)
    {
      if(fTriggerMode) continue;
      primaryTrackList.push_back(nPartSaved);
    }

    if(fTriggerMode && chiPrim > fChiPrimaryMaxCut) continue;
//     if(chiPrim > fChiPrimaryMaxCut) continue;
//     if( chiPrim > 1.e6 ) continue;
    particle.SetId(index);
    particles[nPartSaved] = particle;

#ifdef __kfpAtFirstHit__
    if(trackAtLastHit && particlesAtLastHit){
      KFPTrack trackPDGAtLastHit = *trackAtLastHit;

      if(abs(pdg) == 1000020030 || abs(pdg) == 1000020040 || abs(pdg) == 1000020060) {
        trackPDGAtLastHit.SetCharge( trackPDGAtLastHit.Charge()*2.f );
        trackPDGAtLastHit.SetPx( trackPDGAtLastHit.GetPx()*2.f );
        trackPDGAtLastHit.SetPy( trackPDGAtLastHit.GetPy()*2.f );
        trackPDGAtLastHit.SetPz( trackPDGAtLastHit.GetPz()*2.f );
        for(int iIndex=0; iIndex<9; iIndex++){
          const int iC = index2[iIndex];
          trackPDGAtLastHit.SetCovariance( iC, trackPDGAtLastHit.GetCovariance(iC)*2.f );
        }
        for(int iIndex=0; iIndex<6; iIndex++){
          const int iC = index4[iIndex];
          trackPDGAtLastHit.SetCovariance( iC, trackPDGAtLastHit.GetCovariance(iC)*4.f );
        }
      }
      if(abs(pdg) == 1000030060 || abs(pdg) == 1000030070) {
        trackPDGAtLastHit.SetCharge( trackPDGAtLastHit.Charge()*3.f );
        trackPDGAtLastHit.SetPx( trackPDGAtLastHit.GetPx()*3.f );
        trackPDGAtLastHit.SetPy( trackPDGAtLastHit.GetPy()*3.f );
        trackPDGAtLastHit.SetPz( trackPDGAtLastHit.GetPz()*3.f );
        for(int iIndex=0; iIndex<9; iIndex++) {
          const int iC = index2[iIndex];
          trackPDGAtLastHit.SetCovariance( iC, trackPDGAtLastHit.GetCovariance(iC)*3.f );
        }
        for(int iIndex=0; iIndex<6; iIndex++) {
          const int iC = index4[iIndex];
          trackPDGAtLastHit.SetCovariance( iC, trackPDGAtLastHit.GetCovariance(iC)*9.f );
        }
      }
      if(abs(pdg) == 1000040070) {
        trackPDGAtLastHit.SetCharge( trackPDGAtLastHit.Charge()*4.f );
        trackPDGAtLastHit.SetPx( trackPDGAtLastHit.GetPx()*4.f );
        trackPDGAtLastHit.SetPy( trackPDGAtLastHit.GetPy()*4.f );
        trackPDGAtLastHit.SetPz( trackPDGAtLastHit.GetPz()*4.f );
        for(int iIndex=0; iIndex<9; iIndex++) {
          const int iC = index2[iIndex];
          trackPDGAtLastHit.SetCovariance( iC, trackPDGAtLastHit.GetCovariance(iC)*4.f );
        }
        for(int iIndex=0; iIndex<6; iIndex++) {
          const int iC = index4[iIndex];
          trackPDGAtLastHit.SetCovariance( iC, trackPDGAtLastHit.GetCovariance(iC)*16.f );
        }
      }
          
      particlesAtLastHit->at(nPartSaved) = KFParticle(trackPDGAtLastHit, pdg);
      particlesAtLastHit->at(nPartSaved).SetId(index);
    }
#endif

    particlesPdg[nPartSaved] = pdg;

    nPartSaved++;
  }
}

void StKFParticleInterface::FillPIDHistograms(StPicoTrack *gTrack, const std::vector<int>& pdgVector, const bool isTofm2, float m2tof)
{
  float momentum = gTrack->gPtot();
  for(unsigned int iPdg = 0; iPdg<pdgVector.size(); iPdg++)
  {
    int pdg = pdgVector[iPdg];
    const int iTrackHisto = fTrackPdgToHistoIndex[pdg];
    if( ! (iTrackHisto < 0 || iTrackHisto >= NTrackHistoFolders) )
    {
      fHistoMomentumTracks[iTrackHisto] -> Fill(momentum);
      fHistodEdXTracks[iTrackHisto] -> Fill(momentum, gTrack->dEdx());
      if(isTofm2)
      {
        fHistodEdXwithToFTracks[iTrackHisto] -> Fill(momentum, gTrack->dEdx());
        fHistoTofPIDTracks[iTrackHisto] -> Fill(momentum, m2tof);
        
        if(abs(pdg)==211)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(0.139570, fdEdXMode, 1));
          float betaGamma = TMath::Log10(momentum/0.139570);
          float z = gTrack->dEdxPull(0.139570, fdEdXMode, 1)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
          
          betaGamma = TMath::Log10(momentum/5.485799e-4);
          z = gTrack->nSigmaElectron()*gTrack->dEdxError();
          fHistodEdXZ[0]->Fill(betaGamma, z);
        }
        if(abs(pdg)==321)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(0.493677, fdEdXMode, 1));
          float betaGamma = TMath::Log10(momentum/0.493677);
          float z = gTrack->dEdxPull(0.493677, fdEdXMode, 1)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
        }
        if(abs(pdg)==2212)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(0.938272, fdEdXMode, 1));
          float betaGamma = TMath::Log10(momentum/0.938272);
          float z = gTrack->dEdxPull(0.938272, fdEdXMode, 1)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
        }
        if(abs(pdg)==1000010020)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(1.876124, fdEdXMode, 1));
          float betaGamma = TMath::Log10(momentum/1.876124);
          float z = gTrack->dEdxPull(1.876124, fdEdXMode, 1)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
        }
        if(abs(pdg)==1000010030)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(2.809432, fdEdXMode, 1));
          float betaGamma = TMath::Log10(momentum/2.809432);
          float z = gTrack->dEdxPull(2.809432, fdEdXMode, 1)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
        }
        if(abs(pdg)==1000020030)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(2.809413, fdEdXMode, 2));
          float betaGamma = TMath::Log10(momentum/2.809413);
          float z = gTrack->dEdxPull(2.809413, fdEdXMode, 2)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
        }
        if(abs(pdg)==1000020040)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(3.728400, fdEdXMode, 2));
          float betaGamma = TMath::Log10(momentum/3.728400);
          float z = gTrack->dEdxPull(3.728400, fdEdXMode, 2)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
        }
      }
    }
  }
}

void StKFParticleInterface::FillPIDHistograms(StMuTrack *gTrack, const std::vector<int>& pdgVector, const bool isTofm2, float m2tof)
{
  float momentum = gTrack->p().mag();
  for(unsigned int iPdg = 0; iPdg<pdgVector.size(); iPdg++)
  {
    int pdg = pdgVector[iPdg];
    const int iTrackHisto = fTrackPdgToHistoIndex[pdg];
    if( ! (iTrackHisto < 0 || iTrackHisto >= NTrackHistoFolders) )
    {
      fHistoMomentumTracks[iTrackHisto] -> Fill(momentum);
      fHistodEdXTracks[iTrackHisto] -> Fill(momentum, gTrack->dEdx()*1.e6);
      if(isTofm2)
      {
        fHistodEdXwithToFTracks[iTrackHisto] -> Fill(momentum, gTrack->dEdx()*1.e6);
        fHistoTofPIDTracks[iTrackHisto] -> Fill(momentum, m2tof);
        
        if(abs(pdg)==211)
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(0.139570, fdEdXMode, 1));
        if(abs(pdg)==321)
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(0.493677, fdEdXMode, 1));
        if(abs(pdg)==2212)
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(0.938272, fdEdXMode, 1));
        if(abs(pdg)==1000010020)
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(1.876124, fdEdXMode, 1));
        if(abs(pdg)==1000010030)
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(2.809432, fdEdXMode, 1));
        if(abs(pdg)==1000020030)
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(2.809413, fdEdXMode, 2));
        if(abs(pdg)==1000020040)
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(3.728400, fdEdXMode, 2));
      }
    }
  }
}

bool StKFParticleInterface::OpenCharmTrigger() 
{
  bool triggerDMesons = false;
  if(fKFParticleTopoReconstructor->NPrimaryVertices() == 0) return false;
    
  for(unsigned int iParticle=0; iParticle<GetParticles().size(); iParticle++)
  {
    KFParticle particle = GetParticles()[iParticle];
    
    if( abs(particle.GetPDG()) == 421 ||
        abs(particle.GetPDG()) == 429 || 
        abs(particle.GetPDG()) == 420 || 
        abs(particle.GetPDG()) == 411 || 
        abs(particle.GetPDG()) == 431 || 
        abs(particle.GetPDG()) == 4122 ||
        abs(particle.GetPDG()) == 426 )
    {
      KFParticleSIMD tempSIMDPart(particle);
      float_v l,dl;
      KFParticleSIMD pv(fKFParticleTopoReconstructor->GetPrimVertex());
      tempSIMDPart.GetDistanceToVertexLine(pv, l, dl);
      
      if(abs(particle.GetPDG()) == 411)
        triggerDMesons = (l[0] < 0.4);
      else    
        triggerDMesons = (l[0] < 0.2);
    }
  }
  
  return triggerDMesons;
}

void StKFParticleInterface::OpenCharmTriggerCompression(int nTracksTriggered, int nTracksInEvent, bool triggerDMesons) 
{
  static int nTriggeredEvents = 0;
  static int nTracksInEventTriggered = 0;
  static int nTracksInEventTotal = 0;
  static int nEvents = 0;
  nEvents++;
  nTracksInEventTotal += nTracksInEvent;
  if(triggerDMesons)
  {
    nTriggeredEvents++;
    nTracksInEventTriggered += nTracksTriggered;
    std::cout << "N Events " << nEvents << "    N triggered events " << nTriggeredEvents << "    ratio " << (double(nEvents)/double(nTriggeredEvents)) << std::endl;
    std::cout << "N Tracks " << nTracksInEventTotal << "    N triggered events " << nTracksInEventTriggered << "    ratio " << (double(nTracksInEventTotal)/double(nTracksInEventTriggered)) << std::endl;
  }
}

void StKFParticleInterface::ResizeTrackPidVectors(const int nTracks)
{
  for(int iHypothesis=0; iHypothesis<3; iHypothesis++)
  {
    fTrackPidTof[iHypothesis].clear();
    fTrackPidTof[iHypothesis].resize(nTracks, -1);
    
    fTrackPidTpc[iHypothesis].clear();
    fTrackPidTpc[iHypothesis].resize(nTracks, -1);
  }
}

bool StKFParticleInterface::ProcessEvent(StPicoDst* picoDst, std::vector<int>& triggeredTracks)
{
  cout<<"StKFParticleInterface::ProcessEvent - picoDst"<<endl;
  triggeredTracks.resize(0);
  
  //read PV from pico Event
  KFVertex primaryVertex;
  vector<int> primaryTrackList;
    
  StPicoEvent* picoEvent = picoDst->event();
  if(!picoEvent) return 0;
  
  const TVector3 picoPV = picoEvent->primaryVertex();
  const TVector3 picoPVError = picoEvent->primaryVertexError();
  
  KFPVertex primVtx_tmp;
  primVtx_tmp.SetXYZ(picoPV.x(), picoPV.y(), picoPV.z());
  double dx = picoPVError.x();
  double dy = picoPVError.y();
  double dz = picoPVError.z();
  primVtx_tmp.SetCovarianceMatrix( dx*dx, 0, dy*dy, 0, 0, dz*dz );
  primaryVertex = KFVertex(primVtx_tmp);

//   if(!IsGoodPV(primaryVertex)) return 0;
  
  Int_t nGlobalTracks = picoDst->numberOfTracks( );
  
  fParticles.resize(nGlobalTracks*10);
  fNHftHits.resize(nGlobalTracks*10);
  fParticlesPdg.resize(nGlobalTracks*10);
  int nPartSaved = 0;
  int nUsedTracks = 0;
  
  std::vector<float> m2TofArray(nGlobalTracks*2, -1.0e6f);

  fTrackIdToI.clear();
  fTrackIdToI.resize(nGlobalTracks*2);
 

  for (Int_t iTrack = 0; iTrack < nGlobalTracks; iTrack++) 
  {
    StPicoTrack *gTrack = picoDst->track(iTrack);
    if (! gTrack)            continue;
    if (! gTrack->charge())  continue;
    //if (  gTrack->nHitsFit() < 15) continue;
    if (  gTrack->nHitsFit() < 10) continue; //changed by Petr

    // removed - Petr ..since dEdx is not simulated well ... this should be done offline
    //if (  gTrack->dEdxError() < 0.04 || gTrack->dEdxError() > 0.12 ) continue; 

    const int id = gTrack->id();
    if(id >= fTrackIdToI.size()) fTrackIdToI.resize(id+1);
    fTrackIdToI[id] = iTrack;
    
    int nHftHitsInTrack = 0;
    if(gTrack->hasPxl1Hit()) nHftHitsInTrack++;
    if(gTrack->hasPxl2Hit()) nHftHitsInTrack++;
    if(gTrack->hasIstHit()) nHftHitsInTrack++;
//       if(gTrack->hasSstHit()) nHftHitsInTrack++;
    
    //FIXME temporary solution!!!
    nHftHitsInTrack = gTrack->nHitsFit();
    
    if(fCollectTrackHistograms) fTrackHistograms[0]->Fill(nHftHitsInTrack);
    
    // if(fUseHFTTracksOnly && nHftHitsInTrack < 3) continue;
    if(fUseHFTTracksOnly && !gTrack->hasIstHit()) continue;
    
    StPicoTrackCovMatrix *cov = picoDst->trackCovMatrix(iTrack);
    const StDcaGeometry dcaG = cov->dcaGeometry();
    Int_t q = 1; if (gTrack->charge() < 0) q = -1;
    KFPTrack track;
    if( !GetTrack(dcaG, track, q, id) ) continue;
    
    if(fCollectTrackHistograms)
    {
      fTrackHistograms2D[0]->Fill(track.GetP(), gTrack->dEdx());
      if(q>0) fTrackHistograms2D[1]->Fill(track.GetP(), gTrack->dEdx());
      else    fTrackHistograms2D[2]->Fill(track.GetP(), gTrack->dEdx());  
    }
    
    double m2tof = -1.e6;
    bool isTofm2 = false;
#ifdef USETOF
    if(gTrack->bTofPidTraitsIndex() >= 0)
    {
      const StPicoBTofPidTraits* btofPid = picoDst->btofPidTraits(gTrack->bTofPidTraitsIndex());
      double betaTof2 = btofPid->btofBeta() * btofPid->btofBeta();
      if(fabs(betaTof2) > 1.e-6)
      {
        m2tof = track.GetP()*track.GetP()*(1./betaTof2 - 1.);
        
        if(index >= m2TofArray.size()) m2TofArray.resize(index+1);
        m2TofArray[index] = m2tof;
        
        isTofm2 = true;
      }
//       else
//       {
//         const TVector3 & tofPoint  = btofPid->btofHitPos();
//         StPicoPhysicalHelix innerHelix = gTrack->helix(picoEvent->bField());
//         double lengthTof = fabs( innerHelix.pathLength( tofPoint ));
//         
//         double timeTof = btofPid->btof();
//         if(timeTof > 0. && lengthTof > 0.)
//         {
//           m2tof = track.GetP()*track.GetP()*(1./((lengthTof/timeTof/29.9792458)*(lengthTof/timeTof/29.9792458))-1.);
//           isTofm2 = true;
//         }
//       }
      
      if(fCollectTrackHistograms)
      {
        fTrackHistograms2D[3]->Fill(track.GetP(), gTrack->dEdx());
        fTrackHistograms2D[4]->Fill(track.GetP(), m2tof);
      }
    }
#endif
    double dEdXPull[12] = { fabs(gTrack->dEdxPull(5.109989461E-04, fdEdXMode, 1)), //0 - e
                            fabs(gTrack->dEdxPull(0.1056583745, fdEdXMode, 1)),    //1 - mu
                            fabs(gTrack->dEdxPull(0.139570, fdEdXMode, 1)),        //2 - pi
                            fabs(gTrack->dEdxPull(0.493677, fdEdXMode, 1)),        //3 - K
                            fabs(gTrack->dEdxPull(0.938272, fdEdXMode, 1)),        //4 - p
                            fabs(gTrack->dEdxPull(1.876124, fdEdXMode, 1)),        //5 - d
                            fabs(gTrack->dEdxPull(2.809432, fdEdXMode, 1)),        //6 - t
                            fabs(gTrack->dEdxPull(2.809413, fdEdXMode, 2)),        //7 - He3
                            fabs(gTrack->dEdxPull(3.728400, fdEdXMode, 2)),        //8 - He4
                            fabs(gTrack->dEdxPull(1.197449, fdEdXMode, 1)),        //9 - Sigma
                            fabs(gTrack->dEdxPull(1.32171, fdEdXMode, 1)),         //10- Xi
                            fabs(gTrack->dEdxPull(1.67245, fdEdXMode, 1))};        //11- Omega
    
    vector<int> totalPDG = GetPID(m2tof, track.GetP(), q, gTrack->dEdx(), dEdXPull, isTofm2, id);
    
    int nPartSaved0 = nPartSaved;
    unsigned int nPrimaryTracks = primaryTrackList.size();
    AddTrackToParticleList(track, nHftHitsInTrack, id, totalPDG, primaryVertex, primaryTrackList, fNHftHits, fParticlesPdg, fParticles, nPartSaved); 
    
    if(nPartSaved > nPartSaved0) 
      triggeredTracks.push_back(iTrack);
    
    //fill PID histograms if they are created
    //Histograms are filled for secondary tracks only!!!
    if(fCollectPIDHistograms && (nPrimaryTracks == primaryTrackList.size())) 
    {
      vector<int> selectedPdg;
      for(int iPdg = nPartSaved0; iPdg<nPartSaved; iPdg++)
        selectedPdg.push_back(fParticlesPdg[iPdg]);
      FillPIDHistograms(gTrack, selectedPdg, isTofm2, m2tof);
    }

    nUsedTracks++;
  }
  
  fParticles.resize(nPartSaved);
  fParticlesPdg.resize(nPartSaved);
  fNHftHits.resize(nPartSaved);

  if(nUsedTracks==0) return 0;
  const int nPrimary = primaryTrackList.size();
  if(fCollectTrackHistograms)
  {
    fTrackHistograms[2]->Fill( double(nPrimary)/double(nUsedTracks) );
    fTrackHistograms2D[7]->Fill( nPrimary, (nUsedTracks - nPrimary) );
  }
  
  if( fCleanLowPVTrackEvents && ( 10*nPrimary < (nUsedTracks - nPrimary) ) ) return 0;  
  if( fCleanLowPVTrackEvents && sqrt(dx*dx + dy*dy) > 0.45 ) return 0;
  
  if (fIsFixedTarget){
    if( fCleanLowPVTrackEvents && !(primaryVertex.Z() > 196. && primaryVertex.Z() < 204.) ) return 0;
    if( fCleanLowPVTrackEvents && !(primaryVertex.X() > -2.5 && primaryVertex.X() < 2.  ) ) return 0;
    if( fCleanLowPVTrackEvents && !(primaryVertex.Y() > -4.  && primaryVertex.Y() < 0.  ) ) return 0;

  //   if( fCleanLowPVTrackEvents && ( (nUsedTracks - nPrimary) > (50. + 0.75*nPrimary) ) ) return 0;
    if( fCleanLowPVTrackEvents && ( (nUsedTracks - nPrimary) > 150 ) ) return 0;
  }
  else
    { //collider
      if( fCleanLowPVTrackEvents && primaryVertex.GetR() > 2.5 ) return 0;
   // if( fCleanLowPVTrackEvents && fabs(primaryVertex.Z()) > 75. ) return 0;
    }
  
  const Double_t field = picoEvent->bField();  
  SetField(field);

  CleanPV();
  InitParticles();

  //read PV
//   vector<int> emptyPrimaryList(0);
//   AddPV(primaryVertex, emptyPrimaryList);
  AddPV(primaryVertex, primaryTrackList);
  if(fCollectTrackHistograms)
  {
    fTrackHistograms[1]->Fill(sqrt(dx*dx + dy*dy));
    fTrackHistograms2D[5]->Fill( nPartSaved, sqrt(dx*dx + dy*dy) );
    fTrackHistograms2D[6]->Fill( primaryTrackList.size(), sqrt(dx*dx + dy*dy) );
  }  
  //reconstruct short-lived particles
  ReconstructParticles();

  return 1;
}

bool StKFParticleInterface::ProcessEvent(StMuDst* muDst, vector<KFMCTrack>& mcTracks, vector<int>& mcIndices,
 StKFParticleAnalysisMaker::cProcessSignal processSignal)
{  
  cout<<"StKFParticleInterface::ProcessEvent - MuDst, processSignal="<<processSignal<<endl;
  cout<<" muDst->numberOfGlobalTracks()="<< muDst->numberOfGlobalTracks()<<endl;
  cout<<" muDst->numberOfMcTracks()="<< muDst->numberOfMcTracks()<<endl;
 
  mcTracks.resize(muDst->numberOfMcTracks());
  for (unsigned int iMCTrack=0; iMCTrack<muDst->numberOfMcTracks(); iMCTrack++) 
  {
    StMuMcTrack *mcTrack = muDst->MCtrack(iMCTrack);
    if (! mcTrack) continue;    
    KFMCTrack &mcTrackKF = mcTracks[iMCTrack];
    mcTrack->FillKFMCTrack(mcTrackKF);
    mcTrackKF.SetNMCPixelPoints(mcTrack->No_ist_hit() + mcTrack->No_ssd_hit() + mcTrack->No_pix_hit());
  }
  
  //read PV
  KFVertex primaryVertex;
  vector<int> primaryTrackList;

  float bestRank=-1000000;
  int bestPV=-100000;
  double dx = 0., dy = 0., dz = 0.;
  cout<<" mu - # primary vertices="<<muDst->numberOfPrimaryVertices()<<endl;
  for(unsigned int iPV=0; iPV<muDst->numberOfPrimaryVertices(); iPV++) 
  {
    StMuPrimaryVertex *Vtx = muDst->primaryVertex(iPV);
    if(!Vtx) continue;
    cout<<" rank for i="<<iPV<<" is "<<Vtx->ranking()<<endl;
    //cout<<*Vtx<<endl;
    if (bestRank < Vtx->ranking()) {
      bestRank = Vtx->ranking();
      bestPV = iPV;
    }
    else continue;
    //convert StMuPrimaryVertex to KFVertex
    KFPVertex kfVertex;
    kfVertex.SetXYZ(Vtx->position().x(), Vtx->position().y(), Vtx->position().z());
    dx = Vtx->posError().x();
    dy = Vtx->posError().y();
    dz = Vtx->posError().z();
    kfVertex.SetCovarianceMatrix( dx*dx, 0, dy*dy, 0, 0, dz*dz );
    UShort_t noTracks = Vtx->noTracks();
    kfVertex.SetNContributors(noTracks);
    kfVertex.SetChi2(Vtx->chiSquared());
    primaryVertex = KFVertex(kfVertex);
    primaryVertex.SetId(bestPV); 
  } 
   cout<<"  bestPV="<<bestPV<<endl;
   cout<<"   "<<primaryVertex<<endl;
//   if(!IsGoodPV(primaryVertex)) return 0;

  Int_t nGlobalTracks = muDst->numberOfGlobalTracks();
  fParticles.resize(nGlobalTracks*10);
#ifdef __kfpAtFirstHit__
  fParticlesAtLastHit.resize(nGlobalTracks*10);
#endif
  fNHftHits.resize(nGlobalTracks*10);
  fParticlesPdg.resize(nGlobalTracks*10);
  int nPartSaved = 0;
  int nUsedTracks = 0;
  
  fTrackIdToI.clear();
  fTrackIdToI.resize(nGlobalTracks*2);
 
  for (Int_t iTrack = 0; iTrack < nGlobalTracks; iTrack++) 
  {
    StMuTrack *gTrack = muDst->globalTracks(iTrack);
    if (! gTrack)            continue;
    if (! gTrack->charge())  continue;
    if (  gTrack->flag() < 100 ||  gTrack->flag()%100 == 11) continue; // bad fit or short track pointing to EEMC
    if (  gTrack->flag() > 1000) continue;  // pile up track in TPC
    //if (  gTrack->nHitsFit() < 15) continue;
    if (  gTrack->nHitsFit() < 10) continue; //changed by Petr

    // removed - Petr ..since dEdx is not simulated well ... this should be done offline
    //if (  gTrack->probPidTraits().dEdxErrorFit() < 0.04 || gTrack->probPidTraits().dEdxErrorFit() > 0.12 ) continue;
    
    int nHftHitsInTrack = gTrack->nHitsFit(kIstId) + gTrack->nHitsFit(kSsdId) + gTrack->nHitsFit(kPxlId);
    if(fCollectTrackHistograms) fTrackHistograms[0]->Fill(nHftHitsInTrack);
    if(fUseHFTTracksOnly && nHftHitsInTrack < 3) continue;
    
    const int id = gTrack->id();
    if(id >= fTrackIdToI.size()) fTrackIdToI.resize(id+1);
    fTrackIdToI[id] = iTrack;

    int mcIndex = gTrack->idTruth()-1;
    //This should be equivalne to check IdTruth<10000, above 10000 are real tracks
    if(mcIndex >= int(mcTracks.size()))
      mcIndex = -1;

    if(mcIndex > -1) { //mc(matched) tracks
      if(processSignal==StKFParticleAnalysisMaker::kRealTracksOnly) continue;
    }
    else //real (not matched) track
     if(processSignal==StKFParticleAnalysisMaker::kMcTracksOnly) continue; 
    
    Int_t q = 1; if (gTrack->charge() < 0) q = -1;
    //cout<<" iMuDst="<<iTrack<<" id="<<id<<" mu_pt="<<gTrack->pt()<<endl;
#ifdef __kfpAtFirstHit__
    KFPTrack track = gTrack->kfpTrackAtFirstHit();
    KFPTrack trackAtLastHit  = gTrack->kfpTrackAtLastHit();
    //cout<<"  KF at first point: pt="<<track.GetPt()<<"  R="<<R(track.GetX(),track.GetY())<<endl;
    //cout<<track<<endl;
    //cout<<"  KF at last  point: pt="<<trackAtLastHit.GetPt()<<"  R="<<R(trackAtLastHit.GetX(),trackAtLastHit.GetY())<<endl;
    //cout<<trackAtLastHit<<endl;
   #else
    Int_t dcaGeometryIndex = gTrack->index2Cov();
    if (dcaGeometryIndex < 0) continue;
    StDcaGeometry *dcaG = StMuDst::instance()->covGlobTracks(dcaGeometryIndex);
    if (! dcaG) continue;
      
    KFPTrack track;

    if( !GetTrack(*dcaG, track, q, id) ) continue;
    
    //cout<<endl; 
#endif

    if(fCollectTrackHistograms)
    {
      fTrackHistograms2D[0]->Fill(track.GetP(), gTrack->dEdx()*1.e6);
      if(q>0) fTrackHistograms2D[1]->Fill(track.GetP(), gTrack->dEdx()*1.e6);
      else    fTrackHistograms2D[2]->Fill(track.GetP(), gTrack->dEdx()*1.e6);  
    }
    
    const StMuBTofPidTraits &btofPid = gTrack->btofPidTraits();
    double timeTof = btofPid.timeOfFlight();
    double lengthTof = btofPid.pathLength();
//     if(lengthTof < 0.)
//     {
//       const StThreeVectorF & tofPoint  = btofPid.position();
//       const StThreeVectorF & dcaPoint  = gTrack->dca(bestPV);
//       StPhysicalHelixD innerHelix = gTrack->helix();
//       double dlDCA = fabs( innerHelix.pathLength( StThreeVector<double>(dcaPoint.x(), dcaPoint.y(), dcaPoint.z()) ) );
//       StPhysicalHelixD outerHelix = gTrack->outerHelix();
//       double dlTOF = fabs( outerHelix.pathLength( StThreeVector<double>(tofPoint.x(), tofPoint.y(), tofPoint.z()) ) );
//       
//       double l = gTrack->length();
//       lengthTof = l + dlDCA + dlTOF;
//     }
//     if(lengthTof < 0.)
//     {
//       const StThreeVectorF & tofPoint  = btofPid.position();
//       StPhysicalHelixD innerHelix = dcaG->helix();
//       StMuPrimaryVertex *pv = muDst->primaryVertex(bestPV);
//       lengthTof = fabs( innerHelix.pathLength( StThreeVector<double>(tofPoint.x(), tofPoint.y(), tofPoint.z()) ) ) + 
//                   fabs( innerHelix.pathLength( StThreeVector<double>(pv->position().x(), 
//                                                                      pv->position().y(), 
//                                                                      pv->position().z()) ) );
//     }
    double m2tof = -1.e6;
    bool isTofm2 = false;
#ifdef USETOF
    if(timeTof > 0. && lengthTof > 0.)
    {
      m2tof = track.GetP()*track.GetP()*(1./((lengthTof/timeTof/29.9792458)*(lengthTof/timeTof/29.9792458))-1.);
      isTofm2 = true;
      
      if(fCollectTrackHistograms)
      {
        fTrackHistograms2D[3]->Fill(track.GetP(), gTrack->dEdx()*1.e6);
        fTrackHistograms2D[4]->Fill(track.GetP(), m2tof);
      }
    }
#endif
    double dEdXPull[12] = { fabs(gTrack->nSigmaElectron()),                     //0 - e
                            fabs(gTrack->dEdxPull(0.1056583745, fdEdXMode, 1)), //1 - mu
                            fabs(gTrack->nSigmaPion()),                         //2 - pi
                            fabs(gTrack->nSigmaKaon()),                         //3 - K
                            fabs(gTrack->nSigmaProton()),                       //4 - p
                            fabs(gTrack->dEdxPull(1.876124, fdEdXMode, 1)),     //5 - d
                            fabs(gTrack->dEdxPull(2.809432, fdEdXMode, 1)),     //6 - t
                            fabs(gTrack->dEdxPull(2.809413, fdEdXMode, 2)),     //7 - He3
                            fabs(gTrack->dEdxPull(3.728400, fdEdXMode, 2)),     //8 - He4
                            fabs(gTrack->dEdxPull(1.197449, fdEdXMode, 1)),     //9 - Sigma
                            fabs(gTrack->dEdxPull(1.32171, fdEdXMode, 1)),      //10- Xi
                            fabs(gTrack->dEdxPull(1.67245, fdEdXMode, 1))       //11- Omega
    };    
    
    vector<int> totalPDG = GetPID(m2tof, track.GetP(), q, gTrack->dEdx()*1.e6, dEdXPull, isTofm2, id);
    
    int nPartSaved0 = nPartSaved;
    unsigned int nPrimaryTracks = primaryTrackList.size();
    //cout<<" adding track .."<<iTrack<<endl;
#ifdef __kfpAtFirstHit__
    AddTrackToParticleList(track, nHftHitsInTrack, id, totalPDG, primaryVertex, primaryTrackList, fNHftHits, fParticlesPdg, fParticles, nPartSaved, &trackAtLastHit, &fParticlesAtLastHit);
  #else
     AddTrackToParticleList(track, nHftHitsInTrack, id, totalPDG, primaryVertex, primaryTrackList, fNHftHits, fParticlesPdg, fParticles, nPartSaved);
#endif

    
    if(nPartSaved > nPartSaved0)
      mcIndices[id] = mcIndex;

    //fill PID histograms if they are created
    //Histograms are filled for secondary tracks only!!!
    if(fCollectPIDHistograms && (nPrimaryTracks == primaryTrackList.size())) 
    {
      vector<int> selectedPdg;
      for(int iPdg = nPartSaved0; iPdg<nPartSaved; iPdg++)
        selectedPdg.push_back(fParticlesPdg[iPdg]);
      FillPIDHistograms(gTrack, selectedPdg, isTofm2, m2tof);
    }
    
    nUsedTracks++;
  }
  cout<<"  nUsedTracks="<<nUsedTracks<<endl;
  cout<<"  nPartSaved="<<nPartSaved<<endl;

  fParticles.resize(nPartSaved);
#ifdef __kfpAtFirstHit__
  fParticlesAtLastHit.resize(nPartSaved);
#endif
  fParticlesPdg.resize(nPartSaved);
  fNHftHits.resize(nPartSaved);
  //cout<<"chp0"<<endl;
 
  if(nUsedTracks==0) return 0;
  const int nPrimary = primaryTrackList.size();
  if(fCollectTrackHistograms)
  {
    fTrackHistograms[2]->Fill( double(nPrimary)/double(nUsedTracks) );
    fTrackHistograms2D[7]->Fill( nPrimary, (nUsedTracks - nPrimary) );
  }
  
  //cout<<"nPrimary="<<nPrimary<<" nUsedTracks=="<<nUsedTracks<<endl;
  if( fCleanLowPVTrackEvents && ( 10*nPrimary < (nUsedTracks - nPrimary) ) ) return 0;  
  //cout<<"chp0.1"<<endl;
  if( fCleanLowPVTrackEvents && sqrt(dx*dx + dy*dy) > 0.45 ) return 0;
  //cout<<"chp1"<<endl;
  if (fIsFixedTarget){
    if( fCleanLowPVTrackEvents && !(primaryVertex.Z() > 196. && primaryVertex.Z() < 204.) ) return 0;
    if( fCleanLowPVTrackEvents && !(primaryVertex.X() > -2.5 && primaryVertex.X() < 2.  ) ) return 0;
    if( fCleanLowPVTrackEvents && !(primaryVertex.Y() > -4.  && primaryVertex.Y() < 0.  ) ) return 0;
   //   if( fCleanLowPVTrackEvents && ( (nUsedTracks - nPrimary) > (50. + 0.75*nPrimary) ) ) return 0;
    if( fCleanLowPVTrackEvents && ( (nUsedTracks - nPrimary) > 150 ) ) return 0;
  }
  else
  {// collider
    if( fCleanLowPVTrackEvents && primaryVertex.GetR() > 2.5 ) return 0;
  //   if( fCleanLowPVTrackEvents && fabs(primaryVertex.Z()) > 75. ) return 0;
  }

 // cout<<"chp2"<<endl;
  const Double_t field = muDst->event()->magneticField();
  SetField(field);

  CleanPV();
#ifdef __kfpAtFirstHit__
   InitParticlesAtFirstAndLastHit();
#else
  InitParticles();
#endif
 //cout<<"chp3"<<endl;
 
  //read PV
  AddPV(primaryVertex, primaryTrackList);
  if(fCollectTrackHistograms)
  {
    fTrackHistograms[1]->Fill(sqrt(dx*dx + dy*dy));
    fTrackHistograms2D[5]->Fill( nPartSaved, sqrt(dx*dx + dy*dy) );
    fTrackHistograms2D[6]->Fill( primaryTrackList.size(), sqrt(dx*dx + dy*dy) );
  }  
  //reconstruct short-lived particles
  //cout<<"------------------pred reco particle----------"<<endl;
  //tohle je prazdne
  /*
  Int_t NPV =  fKFParticleTopoReconstructor->NPrimaryVertices();
  for(unsigned int iParticle=0; iParticle<GetParticles().size(); iParticle++) {
    cout<<"i="<<iParticle<<" KFP_iD="<<GetParticles()[iParticle].Id()<<" Mu_iD="<<GetParticles()[iParticle].DaughterIds()[0]
    <<" PDG="<<GetParticles()[iParticle].GetPDG()<<" pt="<< GetParticles()[iParticle].GetPt()<<" R="<<R(GetParticles()[iParticle].GetX(),GetParticles()[iParticle].GetY())<<" phi="<< GetParticles()[iParticle].GetPhi()<<endl;
    cout<<"  "<<GetParticles()[iParticle]<<endl;
   }
 cout<<"-----------------pred reco particle..DONE----------"<<endl;
 */ 
   //cout<<"chp4"<<endl;
 
   ReconstructParticles();
   
  //cout<<"-----------------Po reco particle----------"<<endl;
  int NPV =  fKFParticleTopoReconstructor->NPrimaryVertices();
  int npt=GetParticles().size();
  //cout<<"GetParticles().size()="<<npt<<endl;
  for(unsigned int iParticle=0; iParticle<npt; iParticle++) {
    //cout<<"i="<<iParticle<<" KFP_iD="<<GetParticles()[iParticle].Id()<<" Mu_iD="<<GetParticles()[iParticle].DaughterIds()[0]   
    //<<" PDG="<<GetParticles()[iParticle].GetPDG()<<" pt="<< GetParticles()[iParticle].GetPt()<<" R="<<R(GetParticles()[iParticle].GetX(),GetParticles()[iParticle].GetY())<<" phi="<< GetParticles()[iParticle].GetPhi()<<endl;
    //cout<<"  "<<GetParticles()[iParticle]<<endl;
   }
    //cout<<"-----------------Po reco particle..DONE----------"<<endl;
cout<<"StKFParticleInterface::ProcessEvent..DONE"<<endl;
 
  return 1;
}
