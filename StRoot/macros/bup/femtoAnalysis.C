/*
  root.exe -q -b -x 'muMc.C(1e6,"../*MuDst.root")'
*/
void femtoAnalysis(Int_t N = 1000000, 
                   const Char_t *input = "/gpfs01/star/pwg/fisyak/Femto/2016/100/17100030/*.femtoDst.root", 
                   const Char_t *output = "femto.root", 
                   int year = 2016,
                   TString flowFileList = "") // list of files with centrality and reaction plane, files should be separated with ";"
{
#if !defined(__CINT__)
  std::cout << "This code cannot be compiled" << std::endl;
#else
  //  gSystem->SetFPEMask(kInvalid | kDivByZero | kOverflow );
  gROOT->LoadMacro("lMuDst.C");
  lMuDst(-1,input,"ry2016,RpicoDst,mysql,kfpAna,quiet,nodefault",output);
  
//   StKFParticleInterface::instance()->SetChiPrimaryCut(10);
//   StKFParticleInterface::instance()->SetPtCutCharm(0.7);
//   StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(2);
//   StKFParticleInterface::instance()->SetSoftKaonPIDMode();
//   StKFParticleInterface::instance()->SetSoftTofPidMode();
  std::cout << "KFParticleAnalysis: running analysis for the year " << year << "." << std::endl; 
  StKFParticleAnalysisMaker* kfpAnalysis = (StKFParticleAnalysisMaker*) StMaker::GetTopChain()->Maker("KFParticleAnalysis");
  kfpAnalysis->AnalyseDsPhiPi();
  if(year == 2016)
  {
    kfpAnalysis->UseTMVA();
    
//     kfpAnalysis->SetTMVABinsD0("0:2:3:4:5:6:7:8:9","-1:1000");
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_0_1_pt0_80_BDT.weights.xml", 0.075, 0);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_2_2_pt0_80_BDT.weights.xml", 0.05,  1);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_3_3_pt0_80_BDT.weights.xml", 0.05,  2);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_4_4_pt0_80_BDT.weights.xml", 0.1,   3);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_5_5_pt0_80_BDT.weights.xml", 0.1,   4);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_6_6_pt0_80_BDT.weights.xml", 0.125, 5);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_7_7_pt0_80_BDT.weights.xml", 0.125, 6);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_8_8_pt0_80_BDT.weights.xml", 0.125, 7);
//     
//     kfpAnalysis->SetTMVABinsD0KK("0:2:3:4:5:6:7:8:9","-1:1000");
//     kfpAnalysis->SetTMVAcutsD0KK("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0KK/weights/TMVAClassification_0_1_pt0_80_BDT.weights.xml", 0.075, 0);
//     kfpAnalysis->SetTMVAcutsD0KK("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0KK/weights/TMVAClassification_2_2_pt0_80_BDT.weights.xml", 0.075, 1);
//     kfpAnalysis->SetTMVAcutsD0KK("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0KK/weights/TMVAClassification_3_3_pt0_80_BDT.weights.xml", 0.1  , 2);
//     kfpAnalysis->SetTMVAcutsD0KK("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0KK/weights/TMVAClassification_4_4_pt0_80_BDT.weights.xml", 0.1  , 3);
//     kfpAnalysis->SetTMVAcutsD0KK("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0KK/weights/TMVAClassification_5_5_pt0_80_BDT.weights.xml", 0.1  , 4);
//     kfpAnalysis->SetTMVAcutsD0KK("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0KK/weights/TMVAClassification_6_6_pt0_80_BDT.weights.xml", 0.125, 5);
//     kfpAnalysis->SetTMVAcutsD0KK("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0KK/weights/TMVAClassification_7_7_pt0_80_BDT.weights.xml", 0.15 , 6);
//     kfpAnalysis->SetTMVAcutsD0KK("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0KK/weights/TMVAClassification_8_8_pt0_80_BDT.weights.xml", 0.15 , 7);
// 
//     kfpAnalysis->SetTMVABinsD04("0:3:6:9","-1:1000");
//     kfpAnalysis->SetTMVAcutsD04  ("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D04/weights/TMVAClassification_0_2_pt20_100_BDT.weights.xml", -0.025, 0);
//     kfpAnalysis->SetTMVAcutsD04  ("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D04/weights/TMVAClassification_3_5_pt20_100_BDT.weights.xml", -0.025, 1);
//     kfpAnalysis->SetTMVAcutsD04  ("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D04/weights/TMVAClassification_6_8_pt20_100_BDT.weights.xml",  0    , 2);
// 
//     kfpAnalysis->SetTMVABinsDPlus("0:2:3:4:5:6:7:8:9","-1:1000");
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/weights/TMVAClassification_0_1_pt0_80_BDT.weights.xml", 0.125, 0);
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/weights/TMVAClassification_2_2_pt0_80_BDT.weights.xml", 0.125, 1);
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/weights/TMVAClassification_3_3_pt0_80_BDT.weights.xml", 0.075, 2);
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/weights/TMVAClassification_4_4_pt0_80_BDT.weights.xml", 0.075, 3);
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/weights/TMVAClassification_5_5_pt0_80_BDT.weights.xml", 0.075, 4);
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/weights/TMVAClassification_6_6_pt0_80_BDT.weights.xml", 0.125, 5);
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/weights/TMVAClassification_7_7_pt0_80_BDT.weights.xml", 0.125, 6);
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/weights/TMVAClassification_8_8_pt0_80_BDT.weights.xml", 0.15 , 7);
// 
//     kfpAnalysis->SetTMVABinsDs("0:3:6:9","-1:1000");
//     kfpAnalysis->SetTMVAcutsDs("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/Ds/weights/TMVAClassification_0_2_pt05_100_BDT.weights.xml", 0.05 , 0);
//     kfpAnalysis->SetTMVAcutsDs("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/Ds/weights/TMVAClassification_3_5_pt05_100_BDT.weights.xml", 0.1  , 1);
//     kfpAnalysis->SetTMVAcutsDs("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/Ds/weights/TMVAClassification_6_8_pt05_100_BDT.weights.xml", 0.175, 2);
// 
//     kfpAnalysis->SetTMVABinsLc("0:3:6:9","-1:1000");
//     kfpAnalysis->SetTMVAcutsLc("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/Lc/weights/TMVAClassification_0_2_pt20_100_BDT.weights.xml", 0.15 , 0);
//     kfpAnalysis->SetTMVAcutsLc("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/Lc/weights/TMVAClassification_3_5_pt20_100_BDT.weights.xml", 0.   , 1);
//     kfpAnalysis->SetTMVAcutsLc("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/Lc/weights/TMVAClassification_6_8_pt20_100_BDT.weights.xml", 0.025, 2);

    
    kfpAnalysis->SetTMVABinsDPlus("0:2:4:6:7:9","1:2:3:5:8");
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_0_1_pt10_20_BDT.weights.xml", -0.025, 0, 0);  
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_0_1_pt20_30_BDT.weights.xml", 0.175 , 0, 1);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_0_1_pt30_50_BDT.weights.xml", 1     , 0, 2);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_0_1_pt30_50_BDT.weights.xml", 1     , 0, 3);
        
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_2_3_pt10_20_BDT.weights.xml", 0.05  , 1, 0);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_2_3_pt20_30_BDT.weights.xml", 0.1   , 1, 1);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_2_3_pt30_50_BDT.weights.xml", -0.025, 1, 2);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_2_3_pt50_80_BDT.weights.xml", 0.175 , 1, 3);
        
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_4_5_pt10_20_BDT.weights.xml", 0.05, 2, 0);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_4_5_pt20_30_BDT.weights.xml", 0.05, 2, 1);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_4_5_pt30_50_BDT.weights.xml", 0   , 2, 2);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_4_5_pt50_80_BDT.weights.xml", 0.15, 2, 3);

    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_6_6_pt10_20_BDT.weights.xml", 0.05 , 3, 0);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_6_6_pt20_30_BDT.weights.xml", 0.075, 3, 1);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_6_6_pt30_50_BDT.weights.xml", 0    , 3, 2);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_6_6_pt50_80_BDT.weights.xml", 0.05 , 3, 3);
    
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_7_8_pt10_20_BDT.weights.xml", 0.1  , 4, 0);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_7_8_pt20_30_BDT.weights.xml", 0.05 , 4, 1);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_7_8_pt30_50_BDT.weights.xml", 0.025, 4, 2);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/DPlus/cent+pt_dependence/weights/TMVAClassification_7_8_pt50_80_BDT.weights.xml", 0125 , 4, 3);


    kfpAnalysis->SetTMVAcutsD0   ("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2016/D0.xml",    0.1);
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2016/DPlus.xml", 0.075);
    kfpAnalysis->SetTMVAcutsDs   ("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2016/Ds.xml",    0.125);
    kfpAnalysis->SetTMVAcutsLc   ("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2016/Lc.xml",    0.);
    kfpAnalysis->SetTMVAcutsD0KK ("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2016/D0KK.xml",  0.125);
    kfpAnalysis->SetTMVAcutsD04  ("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2016/D04.xml",   0.);
    kfpAnalysis->SetTMVAcutsBPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2016/BPlus.xml",-0.1);
    kfpAnalysis->SetTMVAcutsB0   ("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2016/B0.xml",   -0.1);
    
    kfpAnalysis->RunCentralityAnalysis();
    kfpAnalysis->SetCentralityFile("/gpfs01/star/pwg/mzyzak/Femto/Template/Centrality/centrality_2016.txt");
  }
  if(year == 2014)
  {
    kfpAnalysis->UseTMVA();
    kfpAnalysis->SetTMVAcutsD0(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D0.xml",    0.1);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/DPlus.xml", 0.05);
    kfpAnalysis->SetTMVAcutsDs(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/Ds.xml",    0.075);
    kfpAnalysis->SetTMVAcutsLc(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/Lc.xml",    0.1);
    kfpAnalysis->SetTMVAcutsD0KK( "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D0KK.xml",  0.125);
    kfpAnalysis->SetTMVAcutsD04(  "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D04.xml",  -0.05);
    kfpAnalysis->SetTMVAcutsBPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/BPlus.xml",-0.1);
    kfpAnalysis->SetTMVAcutsB0(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/B0.xml",   -0.1);
    
    kfpAnalysis->RunCentralityAnalysis();
    kfpAnalysis->SetCentralityFile("/gpfs01/star/pwg/mzyzak/Femto/Template/Centrality/centrality_2014.txt");
  }
  
  if(!flowFileList.IsNull())
  {
    kfpAnalysis->RunFlowAnalysis();
    
    TString fileName;
    int firstSymbolOfFileName = 0;      
    unsigned int iCut = 0;
    while(flowFileList.Tokenize(fileName,firstSymbolOfFileName,";"))
    {
      kfpAnalysis->AddFlowFile(fileName.Data());
      iCut++;
    }
  }
  else
  {
    std::cout << "Flow file list is empty, no flow or centrality analysis will be run." << std::endl;
  }
  
  chain->Init();

  StKFParticleInterface::instance()->SetSoftKaonPIDMode();
  StKFParticleInterface::instance()->SetSoftTofPidMode();
  StKFParticleInterface::instance()->SetChiPrimaryCut(10);
  
  StKFParticleInterface::instance()->SetPtCutCharm(0.5);
  StKFParticleInterface::instance()->SetChiPrimaryCutCharm(8);
  StKFParticleInterface::instance()->SetLdLCutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(10);
  StKFParticleInterface::instance()->SetChi2CutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetLdLCutCharm2D(3);
  StKFParticleInterface::instance()->SetChi2TopoCutCharm2D(10);
  StKFParticleInterface::instance()->SetChi2CutCharm2D(3);
  
  StKFParticleInterface::instance()->AddDecayToReconstructionList(  421);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( -421);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(  426);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(  429);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( -429);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(  411);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( -411);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(  431);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( -431);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( 4122);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(-4122);
  
//   StKFParticleInterface::instance()->AddDecayToReconstructionList(  500);
//   StKFParticleInterface::instance()->AddDecayToReconstructionList(  501);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(  521);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( -521);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(  529);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( -529);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(  511);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( -511);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(  519);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( -519);

  StPicoDstMaker* maker = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
  if (! maker) return;
  maker->SetStatus("*",1);
  TChain *tree = maker->chain();
  Long64_t nentries = tree->GetEntries();
  if (nentries <= 0) return;
  Long64_t nevent = N;
  nevent = TMath::Min(nevent,nentries);
  cout << nentries << " events in chain " << nevent << " will be read." << endl;
//   new StGoodTrigger("y2014");
//   chain->SetAttr(".Privilege",1,"StPicoDstMaker::*")
  chain->EventLoop(nevent);
#endif
  
}
