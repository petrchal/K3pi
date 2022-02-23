// $Id: StKFParticleAnalysisMaker.h,v 1.16 2014/08/06 11:43:53 jeromel Exp $
/*!
 * \class  StKFParticleAnalysisMaker
 * \author Maksym Zyzak
 * \date   2017/10/17
 * \brief  class for analysis of PicoDst
 */                                                                      
#ifndef STAR_StKFParticleAnalysisMaker
#define STAR_StKFParticleAnalysisMaker
//#define __DEVT__
#ifndef StMaker_H
#include "StMaker.h"
#endif
#include "TMVA/Reader.h"
#include "TClonesArray.h"

#include "KFParticle.h"
#include "K3pi.h"

class StKFParticleInterface;
class StKFParticlePerformanceInterface;
// class KFParticle;
class StPicoDst;
class StMuDst;
class TNtuple;
class TFile;
class TChain;
class TTree;



class StKFParticleAnalysisMaker : public StMaker {
  public: 
  StKFParticleAnalysisMaker(const char *name="KFParticleAnalysis");
  virtual       ~StKFParticleAnalysisMaker();
  virtual Int_t  Init();
  virtual Int_t  InitRun(Int_t runumber);
  void           BookVertexPlots();
  virtual Int_t  Make();
  virtual Int_t  Finish();
  Bool_t         Check();
  void AnalysePicoDst() { fIsPicoAnalysis = true;  }
  void AnalyseMuDst()   { fIsPicoAnalysis = false; }
  static void    PrintMem(const Char_t *opt = "");
  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $ $Id: StKFParticleAnalysisMaker.h,v 1.0 2017/10/07 11:43:53 mzyzak Exp $ built " __DATE__ " " __TIME__ ; 
    return cvs;
  }

  enum cProcessSignal {kAllTracks=0,kMcTracksOnly=1,kRealTracksOnly=2};
  //note: MC means mathced to MC tracks ... with propper TruthId
  //note: real ... without matching to MC tracks
  //so this does not adress ghosting, etc...
  void SetProcessSignal(cProcessSignal v) { fProcessSignal = v; }
  void CollectTrackHistograms() { fCollectTrackHistograms = true; }
  void CollectPIDHistograms() { fCollectPIDHistograms = true; }
    
  void RunKaonAnalysis()         { fKaonAnalysis = true; }
  void SetKaonFile(TString file) { fKaonFileName = file;}
  
  
  void AddDecayToReconstructionList( int iDecay );
  
  void Produce3DEfficiencyFile() { fIsProduce3DEfficiencyFile = true; }
  void Set3DEfficiency(TString fileName) { f3DEfficiencyFile = fileName; }

  void AddCandidateToStore(int pdg);


 private:
  static const int fNNTuples = 8;
  Char_t                mBeg[1];        //!
  StMuDst                          *fMuDst;
  StPicoDst                        *fPicoDst;                          //!
  StKFParticleInterface            *fStKFParticleInterface;            //!
  StKFParticlePerformanceInterface *fStKFParticlePerformanceInterface; //!
  int fNTuplePDG[fNNTuples];
  Char_t                mEnd[1];        //!
  
  //general setup options - important
  Bool_t fIsPicoAnalysis; //run on picoDst or MuDst
  cProcessSignal fProcessSignal;  //use real or MC matched tacks
  
  //kaon analysis
  Bool_t fKaonAnalysis;
  TK3pi fK; //for global sharing of currently selected particle
  TString fKaonFileName;
  TFile* fKaonFile;
  TTree* fKaonTree;

  
   std::vector<int> fDecays;

  //interafce settings
  Bool_t fCollectTrackHistograms;
  Bool_t fCollectPIDHistograms;
  bool fIsProduce3DEfficiencyFile;
  TString f3DEfficiencyFile;

  void GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle);
  void GetParticleParameters(const int iReader, KFParticle& particle);
  long GetUniqueEventId(const int iRun, const int iEvent) const;
  
  
  //K->3p
  bool FillKFDaughters(KFParticle &particle);
  void Fill_KaonNtuples();
  void StKFParticleAnalysisMaker::MatchMotherKaon(KFParticle& particle);
 
  bool fStoreCandidates;
  KFParticle fPartcileCandidate;
  std::vector<bool> fIsStoreCandidate;
  TFile* fCandidateFile;
  TTree* fCandidatesTree;

 
  public:

  ClassDef(StKFParticleAnalysisMaker,0)   //
};





#endif
// $Log: StKFParticleAnalysisMaker.h,v $
