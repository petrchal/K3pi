/*
  root.exe -q -b -x 'muMc.C(1e6,"../*MuDst.root")'
*/
void kfpAnalysis(
     Int_t N = 1e9, 
     char TrkType=0,  //processignal flag: {kAllTracks=0,kMcTracksOnly=1,kRealTracksOnly=2;
     char isPico=false, 
     char isFXT=false, 
     char noPID=true, 
     const Char_t *input = ,
     const Char_t *output,
     const Char_t *triggerSet //force user to set it ...  = "y2019"
     ) {
#if !defined(__CINT__)
  std::cout << "This code cannot be compiled" << std::endl;
#else
  //  gSystem->SetFPEMask(kInvalid | kDivByZero | kOverflow );
  gROOT->LoadMacro("lMuDst.C");
  TString Chain("r");
  Chain += triggerSet;
  if (! isPico) Chain += ",RMuDst";
  //if (! isPico) {Chain += ",RMuDst,PicoWrite";}//  isPico = kTRUE;}
  else          Chain += ",RpicoDst";
  Chain += ",kfpAna,mysql,detDb,nodefault,quiet";
  lMuDst(-1,input,Chain,output);

  StKFParticleAnalysisMaker* kfpAnalysis = (StKFParticleAnalysisMaker*) StMaker::GetTopChain()->Maker("KFParticleAnalysis");
  if (isPico) kfpAnalysis->AnalysePicoDst();
    else kfpAnalysis->AnalyseMuDst();
 
  kfpAnalysis->RunKaonAnalysis();
  TString fname="kaon_";fname+=output;
  kfpAnalysis->SetKaonFile(fname);

  kfpAnalysis->SetProcessSignal(TrkType);  //enable for simulations

//   kfpAnalysis->CollectPIDHistograms();
  kfpAnalysis->CollectTrackHistograms();

  kfpAnalysis->AddDecayToReconstructionList( 310); //Kshort
  kfpAnalysis->AddDecayToReconstructionList( 100321); //K->3pi
  kfpAnalysis->AddDecayToReconstructionList(-100321);
  kfpAnalysis->AddDecayToReconstructionList( 200321);
  kfpAnalysis->AddDecayToReconstructionList(-200321);
  kfpAnalysis->AddDecayToReconstructionList( 3122); //lambda
  kfpAnalysis->AddDecayToReconstructionList(-3122);
  kfpAnalysis->AddDecayToReconstructionList( 3312);//xi
  kfpAnalysis->AddDecayToReconstructionList(-3312);
  kfpAnalysis->AddDecayToReconstructionList( 3334);
  kfpAnalysis->AddDecayToReconstructionList(-3334);//omega
  
  /*
  kfpAnalysis->AddDecayToReconstructionList(   22); //gamma
  kfpAnalysis->AddDecayToReconstructionList(  111); //pi0
  kfpAnalysis->AddDecayToReconstructionList(  333); //phi
//   kfpAnalysis->AddDecayToReconstructionList(  313); //K*
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
*/
  
  ((StBFChain *) StMaker::GetTopChain())->Init();
  
 
  StKFParticleInterface::instance()->CleanLowPVTrackEvents();  
  if (isFXT) StKFParticleInterface::instance()->FixedTarget();
  
  StKFParticleInterface::instance()->SetSoftKaonPIDMode();
  if (noPID) StKFParticleInterface::instance()->SetAllIsKaonPIDMode();
  StKFParticleInterface::instance()->SetSoftTofPidMode();

  
  StKFParticleInterface::instance()->SetChiPrimaryCut(12);
  
  StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(1);
  StKFParticleInterface::instance()->SetLCut(0.3f);

    
  StKFParticleInterface::instance()->SetChiPrimaryCut2D(0);
  StKFParticleInterface::instance()->SetChi2Cut2D(3);
  StKFParticleInterface::instance()->SetLdLCut2D(5);
  /*/
  StKFParticleInterface::instance()->SetChiPrimaryCut2D(0);
  StKFParticleInterface::instance()->SetChi2Cut2D(0);
  StKFParticleInterface::instance()->SetLdLCut2D(0);
  */

  //these should be important for 3pi 
  // 3pi vertex quality - saved into mother_chi2ndf 
  StKFParticleInterface::instance()->SetChi2CutXiOmega(30);//30 - quite wide cut
  //chi2 od the 3pi to primary vertex ..mother_PV_chi2  (5->25)
  StKFParticleInterface::instance()->SetChi2TopoCutXiOmega(30); //30  ..
  //decay length/ its sigma, but that's not what I save....
  StKFParticleInterface::instance()->SetLdLCutXiOmega(50);  
  
  /*
  StKFParticleInterface::instance()->SetChi2CutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetLdLCutCharmManybodyDecays(5);
*/
  /*dont need these for 3pi
  StKFParticleInterface::instance()->SetChi2CutCharmManybodyDecays(0);
  StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(0);
  StKFParticleInterface::instance()->SetLdLCutCharmManybodyDecays(0);
  */
  
  ((StBFChain *) StMaker::GetTopChain())->EventLoop(N);
#endif
  
}
