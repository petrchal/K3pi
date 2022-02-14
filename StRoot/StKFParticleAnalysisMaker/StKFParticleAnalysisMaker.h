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
class StRefMultCorr;



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
  void StoreTMVANtuples() { fStoreTmvaNTuples = true; }
  void CollectTrackHistograms() { fCollectTrackHistograms = true; }
  void CollectPIDHistograms() { fCollectPIDHistograms = true; }
  void UseTMVA() { fTMVAselection = true; }
  void SetTMVABinsD0   (TString centralityBins, TString ptBins) { SetTMVABins(0, centralityBins, ptBins); }
  void SetTMVABinsDPlus(TString centralityBins, TString ptBins) { SetTMVABins(1, centralityBins, ptBins); }
  void SetTMVABinsDs   (TString centralityBins, TString ptBins) { SetTMVABins(2, centralityBins, ptBins); }
  void SetTMVABinsLc   (TString centralityBins, TString ptBins) { SetTMVABins(3, centralityBins, ptBins); }
  void SetTMVABinsD0KK (TString centralityBins, TString ptBins) { SetTMVABins(4, centralityBins, ptBins); }
  void SetTMVABinsD04  (TString centralityBins, TString ptBins) { SetTMVABins(5, centralityBins, ptBins); }
  void SetTMVABinsBPlus(TString centralityBins, TString ptBins) { SetTMVABins(6, centralityBins, ptBins); }
  void SetTMVABinsB0   (TString centralityBins, TString ptBins) { SetTMVABins(7, centralityBins, ptBins); }
  void SetTMVAcutsD0   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[0][iCentralityBin][iPtBin] = file; fTMVACut[0][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsDPlus(TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[1][iCentralityBin][iPtBin] = file; fTMVACut[1][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsDs   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[2][iCentralityBin][iPtBin] = file; fTMVACut[2][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsLc   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[3][iCentralityBin][iPtBin] = file; fTMVACut[3][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsD0KK (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[4][iCentralityBin][iPtBin] = file; fTMVACut[4][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsD04  (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[5][iCentralityBin][iPtBin] = file; fTMVACut[5][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsBPlus(TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[6][iCentralityBin][iPtBin] = file; fTMVACut[6][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsB0   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[7][iCentralityBin][iPtBin] = file; fTMVACut[7][iCentralityBin][iPtBin] = cut; }
   
  void RunKaonAnalysis()         { fKaonAnalysis = true; }
  void SetKaonFile(TString file) { fKaonFileName = file;}
  
  void RunCentralityAnalysis() { fRunCentralityAnalysis = true; }
  void SetCentralityFile(TString file) { fCentralityFile = file; }
  
  void AnalyseDsPhiPi() { fAnalyseDsPhiPi = true; }
  
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
  TNtuple* fCutsNTuple[fNNTuples];
  TFile* fNTupleFile[fNNTuples];
  int fNTuplePDG[fNNTuples];
  TString fNtupleNames[fNNTuples];
  TString fNtupleCutNames[fNNTuples];
  std::vector<TString> fDaughterNames[fNNTuples];
  vector< vector<TString> > fTMVACutFile[fNNTuples];
  vector< vector<double> > fTMVACut[fNNTuples];
  vector< vector<TMVA::Reader*> > fTMVAReader[fNNTuples];
  std::vector<int> fTMVACentralityBins[fNNTuples];
  std::vector<double> fTMVAPtBins[fNNTuples];
  Char_t                mEnd[1];        //!
  std::vector<float> fTMVAParticleParameters[fNNTuples];
  int fNTrackTMVACuts;

  //general setup options - important
  Bool_t fIsPicoAnalysis; //run on picoDst or MuDst
  cProcessSignal fProcessSignal;  //use real or MC matched tacks
  //tmwa related options
  Bool_t fStoreTmvaNTuples;
  Bool_t fCollectTrackHistograms;
  Bool_t fCollectPIDHistograms;
  Bool_t fTMVAselection;
  
  //Centrality
  bool fRunCentralityAnalysis;
  StRefMultCorr *fRefmultCorrUtil;
  TString fCentralityFile; //not used so far
  
  //kaon analysis
  Bool_t fKaonAnalysis;
  TK3pi fK; //for global sharing of currently selected particle
  TString fKaonFileName;
  TFile* fKaonFile;
  TTree* fKaonTree;

  
  bool fAnalyseDsPhiPi;
  std::vector<int> fDecays;

  
  bool fIsProduce3DEfficiencyFile;
  TString f3DEfficiencyFile;

  void GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle);
  void GetParticleParameters(const int iReader, KFParticle& particle);
  long GetUniqueEventId(const int iRun, const int iEvent) const;
  
  //charm TMWA
  int GetTMVACentralityBin(int iReader, int centrality);
  int GetTMVAPtBin(int iReader, double pt);
  void SetTMVACentralityBins(int iReader, TString bins);
  void SetTMVAPtBins(int iReader, TString bins);
  void SetTMVABins(int iReader, TString centralityBins="-1:1000", TString ptBins="-1.:1000.");

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
