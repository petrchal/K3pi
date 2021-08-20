/*
  root.exe -q -b -x 'muMc.C(1e6,"../*MuDst.root")'
*/
void muAnalysis(
     Int_t N = 1e9, char isSim=false, 
     char noPID=false, 
     const Char_t *input = "/gpfs01/star/subsys-tpc/fisyak/Pico/2020/TFG21b/RF/9p8GeV_fixedTarget/032/hlt_21032016_10_01_000.MuDst.root", const Char_t *output = "mu.root"
     ) {
#if !defined(__CINT__)
  std::cout << "This code cannot be compiled" << std::endl;
#else
  //  gSystem->SetFPEMask(kInvalid | kDivByZero | kOverflow );
  gROOT->LoadMacro("lMuDst.C");
  lMuDst(-1,input,"ry2016,picoEvt,RMuDst,mysql,kfpAna,quiet,nodefault",output);
  
  StKFParticleAnalysisMaker* kfpAnalysis = (StKFParticleAnalysisMaker*) StMaker::GetTopChain()->Maker("KFParticleAnalysis");
  kfpAnalysis->AnalyseMuDst();
  kfpAnalysis->RunKaonAnalysis();
  TString fname="kaon_";fname+=output;
  kfpAnalysis->SetKaonFile(fname);
  if (isSim) kfpAnalysis->ProcessSignal();  //enable for simulations

//   kfpAnalysis->CollectPIDHistograms();
  kfpAnalysis->CollectTrackHistograms();

  kfpAnalysis->AddDecayToReconstructionList( 310);
  kfpAnalysis->AddDecayToReconstructionList( 100321);
  kfpAnalysis->AddDecayToReconstructionList(-100321);
  kfpAnalysis->AddDecayToReconstructionList( 200321);
  kfpAnalysis->AddDecayToReconstructionList(-200321);
  kfpAnalysis->AddDecayToReconstructionList( 3122);
  kfpAnalysis->AddDecayToReconstructionList(-3122);
  kfpAnalysis->AddDecayToReconstructionList( 3312);
  kfpAnalysis->AddDecayToReconstructionList(-3312);
  kfpAnalysis->AddDecayToReconstructionList( 3334);
  kfpAnalysis->AddDecayToReconstructionList(-3334);
  
  kfpAnalysis->AddDecayToReconstructionList(   22);
  kfpAnalysis->AddDecayToReconstructionList(  111);
  kfpAnalysis->AddDecayToReconstructionList(  333);
//   kfpAnalysis->AddDecayToReconstructionList(  313);
//   kfpAnalysis->AddDecayToReconstructionList( -313);
//   kfpAnalysis->AddDecayToReconstructionList(  323);
//   kfpAnalysis->AddDecayToReconstructionList( -323);
  kfpAnalysis->AddDecayToReconstructionList( 3324);
  kfpAnalysis->AddDecayToReconstructionList(-3324);
  
  kfpAnalysis->AddDecayToReconstructionList( 3000);
  kfpAnalysis->AddDecayToReconstructionList( 3001);
  kfpAnalysis->AddDecayToReconstructionList( 3003);
  kfpAnalysis->AddDecayToReconstructionList( 3103);
  kfpAnalysis->AddDecayToReconstructionList( 3004);
  kfpAnalysis->AddDecayToReconstructionList( 3005);
  kfpAnalysis->AddDecayToReconstructionList( 3006);
  kfpAnalysis->AddDecayToReconstructionList( 3007);
  kfpAnalysis->AddDecayToReconstructionList( 3012);
  kfpAnalysis->AddDecayToReconstructionList( 3013);
/*
  kfpAnalysis->AddDecayToReconstructionList(-7000211);
  kfpAnalysis->AddDecayToReconstructionList(-7000014);
  kfpAnalysis->AddDecayToReconstructionList( 7000211);
  kfpAnalysis->AddDecayToReconstructionList( 7000014);
  kfpAnalysis->AddDecayToReconstructionList(-7000321);
  kfpAnalysis->AddDecayToReconstructionList(-8000014);
  kfpAnalysis->AddDecayToReconstructionList( 7000321);
  kfpAnalysis->AddDecayToReconstructionList( 8000014);
  kfpAnalysis->AddDecayToReconstructionList( 7003112);
  kfpAnalysis->AddDecayToReconstructionList( 7002112);
  kfpAnalysis->AddDecayToReconstructionList(-7003112);
  kfpAnalysis->AddDecayToReconstructionList(-7002112);
  kfpAnalysis->AddDecayToReconstructionList(-7003222);
  kfpAnalysis->AddDecayToReconstructionList(-8002112);
  kfpAnalysis->AddDecayToReconstructionList( 7003222);
  kfpAnalysis->AddDecayToReconstructionList( 8002112);
  kfpAnalysis->AddDecayToReconstructionList( 7003312);
  kfpAnalysis->AddDecayToReconstructionList( 7003122);
  kfpAnalysis->AddDecayToReconstructionList(-7003312);
  kfpAnalysis->AddDecayToReconstructionList(-7003122);
  kfpAnalysis->AddDecayToReconstructionList( 7003334);
  kfpAnalysis->AddDecayToReconstructionList( 7003322);
  kfpAnalysis->AddDecayToReconstructionList(-7003334);
  kfpAnalysis->AddDecayToReconstructionList(-7003322);
  kfpAnalysis->AddDecayToReconstructionList(-9000321);
  kfpAnalysis->AddDecayToReconstructionList(-9000111);
  kfpAnalysis->AddDecayToReconstructionList( 9000321);
  kfpAnalysis->AddDecayToReconstructionList( 9000111);
  kfpAnalysis->AddDecayToReconstructionList( 8003334);
  kfpAnalysis->AddDecayToReconstructionList( 8003122);
  kfpAnalysis->AddDecayToReconstructionList(-8003334);
  kfpAnalysis->AddDecayToReconstructionList(-8003122);
  kfpAnalysis->AddDecayToReconstructionList(-8003222);
  kfpAnalysis->AddDecayToReconstructionList(-8000111);
  kfpAnalysis->AddDecayToReconstructionList( 8003222);
  kfpAnalysis->AddDecayToReconstructionList( 8000111);
*/
  
  ((StBFChain *) StMaker::GetTopChain())->Init();
  
 
  StKFParticleInterface::instance()->CleanLowPVTrackEvents();  
  
  StKFParticleInterface::instance()->SetSoftKaonPIDMode();
  if (noPID) StKFParticleInterface::instance()->SetAllIsKaonPIDMode();
  StKFParticleInterface::instance()->SetSoftTofPidMode();

  
  StKFParticleInterface::instance()->SetChiPrimaryCut(12);
  
  StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(1);
  StKFParticleInterface::instance()->SetLCut(0.3f);

  
  StKFParticleInterface::instance()->SetChiPrimaryCut2D(0);
  StKFParticleInterface::instance()->SetChi2Cut2D(3);
  StKFParticleInterface::instance()->SetLdLCut2D(5);
  
  
  
  StKFParticleInterface::instance()->SetChi2CutXiOmega(30);
  StKFParticleInterface::instance()->SetChi2TopoCutXiOmega(30);
  StKFParticleInterface::instance()->SetLdLCutXiOmega(50);  
  
  StKFParticleInterface::instance()->SetChi2CutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetLdLCutCharmManybodyDecays(5);

  
  
  ((StBFChain *) StMaker::GetTopChain())->EventLoop(N);
#endif
  
}
