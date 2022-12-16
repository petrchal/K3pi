/*
 * This file is part of KFParticle package
 * Copyright (C) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef DO_TPCCATRACKER_EFF_PERFORMANCE

#include "KFParticlePerformanceBase.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"

KFParticlePerformanceBase::KFParticlePerformanceBase():
  fParteff(), fPVeff(), fPVeffMCReconstructable(), outfileName(), histodir(0), fNEvents(0), fStoreMCHistograms(1), 
  fStorePrimSecHistograms(1), fStoreZRHistograms(1), fStore3DEfficiency(0), fHistoDir(0)
{
  /** The default constructor. Initialises all pointers to nullptr.  **/
  for(int iParticle=0; iParticle<KFPartEfficiencies::nParticles; iParticle++)
  {
    for(int iFitQA=0; iFitQA<nFitQA; iFitQA++)
    {
      hFitDaughtersQA         [iParticle][iFitQA] = nullptr;
      hFitQA                  [iParticle][iFitQA] = nullptr;
      hFitQANoConstraint      [iParticle][iFitQA] = nullptr;
      hFitQAMassConstraint    [iParticle][iFitQA] = nullptr;
      hFitQATopoConstraint    [iParticle][iFitQA] = nullptr;
      hFitQATopoMassConstraint[iParticle][iFitQA] = nullptr;
    }
  }

  for(int iParticle=0; iParticle<KFPartEfficiencies::nParticles; iParticle++)
    for(int iQA=0; iQA<nDSToParticleQA; iQA++)
      hDSToParticleQA[iParticle][iQA] = nullptr;
  
  for(int iParameterSet=0; iParameterSet<nParametersSet; iParameterSet++)
  {
    for(int iParticle=0; iParticle<KFPartEfficiencies::nParticles; iParticle++)
    {
      for(int iHisto=0; iHisto<nHistoPartParam; iHisto++)
      {
        hPartParam               [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParamPrimary        [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParamPrimaryMass    [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParamPrimaryTopo    [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParamPrimaryTopoMass[iParameterSet][iParticle][iHisto] = nullptr;
        hPartParamSecondary      [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParamSecondaryMass  [iParameterSet][iParticle][iHisto] = nullptr;
      }
    }
  }

  for(int iParameterSet=0; iParameterSet<nParametersSet; iParameterSet++)
  {
    for(int iParticle=0; iParticle<KFPartEfficiencies::nParticles; iParticle++)
    {
      for(int iHisto=0; iHisto<nHistoPartParam2D; iHisto++)
      {
        hPartParam2D               [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParam2DPrimary        [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParam2DPrimaryMass    [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParam2DPrimaryTopo    [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParam2DPrimaryTopoMass[iParameterSet][iParticle][iHisto] = nullptr;
        hPartParam2DSecondary      [iParameterSet][iParticle][iHisto] = nullptr;
        hPartParam2DSecondaryMass  [iParameterSet][iParticle][iHisto] = nullptr;
      }
    }
  }

  for(int iParticle=0; iParticle<KFPartEfficiencies::nParticles; iParticle++)
    for(int iHisto=0; iHisto<nHistoPartParam3D; iHisto++)
      hPartParam3D[0][iParticle][iHisto] = nullptr;

  for(int iParticle=0; iParticle<KFPartEfficiencies::nParticles; iParticle++)
    for(int iEffSet=0; iEffSet<3; iEffSet++)
      for(int iEff=0; iEff<nPartEfficiency; iEff++)
        hPartEfficiency[iParticle][iEffSet][iEff] = nullptr;
        
  for(int iParticle=0; iParticle<KFPartEfficiencies::nParticles; iParticle++)
    for(int iEffSet=0; iEffSet<3; iEffSet++)
      for(int iEff=0; iEff<nPartEfficiency2D; iEff++)
        hPartEfficiency2D[iParticle][iEffSet][iEff] = nullptr;

  for(int iParticle=0; iParticle<KFPartEfficiencies::nParticles; iParticle++)
    for(int iEffSet=0; iEffSet<4; iEffSet++)
      hPartEfficiencyMulti[iParticle][iEffSet] = nullptr;

  for(int iEffSet=0; iEffSet<2; iEffSet++)
    for(int iHistoPV=0; iHistoPV<nHistosPV; iHistoPV++)
      hPVFitQa[iEffSet][iHistoPV] = nullptr;

  for(int iEffSet1=0; iEffSet1<2; iEffSet1++)
    for(int iEffSet2=0; iEffSet2<2; iEffSet2++)
      for(int iHistoPV=0; iHistoPV<nHistosPV-1; iHistoPV++)
        hPVFitQa2D[iEffSet1][iEffSet2][iHistoPV] = nullptr;

  for(int iParam=0; iParam<nHistosPVParam; iParam++)
  {
    hPVParam      [iParam] = nullptr;
    hPVParamGhost [iParam] = nullptr;
    hPVParamSignal[iParam] = nullptr;
    hPVParamPileup[iParam] = nullptr;
    hPVParamBG    [iParam] = nullptr;
  }
  
  for(int iParam=0; iParam<nHistosPVParam2D; iParam++)
    hPVParam2D[iParam] = nullptr;
  
  for(int iFitQA=0; iFitQA<nFitPVTracksQA; iFitQA++)
    hFitPVTracksQA[iFitQA] = nullptr;
  
  for(int iTP=0; iTP<nHistosTP; iTP++)
    hTrackParameters[iTP] = nullptr;

  for(int iEffSet=0; iEffSet<4; iEffSet++)
    for(int iEff=0; iEff<nPVefficiency; iEff++)
      hPVefficiency[iEffSet][iEff] = nullptr;
}

void KFParticlePerformanceBase::CreateHistos(std::string histoDir, TDirectory* outFile, std::map<int,bool> decays)
{
  /** Creates all histograms. If "outFile" is provided - creates a new ROOT directory and stores all 
   ** histograms there. Otherwise histograms are stored in TDirectory::CurrentDirectory().
   ** \param[in] histoDir - name of the ROOT directory with histograms
   ** \param[in] outFile - pointer to the external ROOT directory or file, where histograms should be stored
   ** \param[in] decays - a list of decays, for which histograms are created, if empty - histograms are
   ** created for all decay channels from the KF Particle Finder reconstruction scheme
   **/
  TDirectory *curdir = gDirectory;
  if (outFile) {
    outFile->cd();
    fHistoDir = outFile;
    if (histoDir != "") {
      fHistoDir = outFile->mkdir( TString(histoDir) );
      fHistoDir->cd();
    }
  } else {
    fHistoDir = TDirectory::CurrentDirectory();
  }
  {
    gDirectory->mkdir("KFParticlesFinder");
    gDirectory->cd("KFParticlesFinder");
    histodir = gDirectory;
    gDirectory->mkdir("Particles");
    gDirectory->cd("Particles");
    for(int iPart=0; iPart<fParteff.nParticles; ++iPart)
    {
      if(!(decays.empty()) && (iPart < fParteff.fFirstStableParticleIndex || iPart > fParteff.fLastStableParticleIndex))
        if(decays.find(fParteff.partPDG[iPart]) == decays.end()) continue;
        
      gDirectory->mkdir(fParteff.partName[iPart].data());
      gDirectory->cd(fParteff.partName[iPart].data());
      {
        if(fStoreMCHistograms)
        {
          TString res = "res";
          TString pull = "pull";

          gDirectory->mkdir("DaughtersQA");
          gDirectory->cd("DaughtersQA");
          {
            TString parName[nFitQA/2] = {"X","Y","Z","Px","Py","Pz","E","M"};
            int nBins = 100;
//             float xMax[nFitQA/2] = {0.15,0.15,0.03,0.01,0.01,0.06,0.06,0.01};
            float xMax[nFitQA/2] = {4.,4.,10.,0.3,0.3,0.3,0.3,0.01};

            for( int iH=0; iH<nFitQA/2; iH++ ){
              hFitDaughtersQA[iPart][iH]   = new TH1F((res+parName[iH]).Data(),
                                                      (GetDirectoryPath()+res+parName[iH]).Data(), 
                                                      nBins, -xMax[iH],xMax[iH]);
              hFitDaughtersQA[iPart][iH+8] = new TH1F((pull+parName[iH]).Data(),
                                                      (GetDirectoryPath()+pull+parName[iH]).Data(), 
                                                      nBins, -6,6);
            }
          }
          gDirectory->cd(".."); //particle directory

          gDirectory->mkdir("DSToParticleQA");
          gDirectory->cd("DSToParticleQA");
          {
            TString parName[3] = {"X","Y","Z"};
            int nBins = 100;
//             float xMax[3] = {0.5, 0.5, 2.};
            float xMax[3] = {4.f, 4.f, 10.f};

            for( int iH=0; iH<3; iH++ ){
              hDSToParticleQA[iPart][iH]   = new TH1F((res+parName[iH]).Data(),
                                                      (GetDirectoryPath()+res+parName[iH]).Data(), 
                                                      nBins, -xMax[iH],xMax[iH]);
              hDSToParticleQA[iPart][iH+3] = new TH1F((pull+parName[iH]).Data(),
                                                      (GetDirectoryPath()+pull+parName[iH]).Data(), 
                                                      nBins, -6,6);
            }
            
            hDSToParticleQA[iPart][6] = new TH1F("r", (GetDirectoryPath()+TString("r")).Data(), 1000, 0.0, 20.0);
          }
          gDirectory->cd(".."); //particle directory
          
          CreateFitHistograms(hFitQA[iPart], iPart);
          CreateEfficiencyHistograms(hPartEfficiency[iPart], hPartEfficiency2D[iPart], hPartEfficiencyMulti[iPart]);
        }
        gDirectory->mkdir("Parameters");
        gDirectory->cd("Parameters");
        {
          const bool drawZR = IsCollectZRHistogram(iPart);
          CreateParameterHistograms(hPartParam[0], hPartParam2D[0], hPartParam3D[0], iPart, drawZR);

          if(IsCollect3DHistogram(iPart))
          {
            gDirectory->mkdir("SignalReco");
            gDirectory->cd("SignalReco");
            {
              CreateParameterHistograms(hPartParam[4], hPartParam2D[4], 0, iPart, drawZR);
            }
            gDirectory->cd(".."); // Parameters
            gDirectory->mkdir("BGReco");
            gDirectory->cd("BGReco");
            {
              CreateParameterHistograms(hPartParam[5], hPartParam2D[5], 0, iPart, drawZR);
            }
            gDirectory->cd(".."); // Parameters
          }
          
          if(fStoreMCHistograms)
          {
            gDirectory->mkdir("Signal");
            gDirectory->cd("Signal");
            {
              CreateParameterHistograms(hPartParam[1], hPartParam2D[1], 0, iPart, drawZR);
            }
            gDirectory->cd(".."); // particle directory / Parameters
            gDirectory->mkdir("Background");
            gDirectory->cd("Background");
            {
              CreateParameterHistograms(hPartParam[2], hPartParam2D[2], 0, iPart, drawZR);
            }
            gDirectory->cd(".."); // particle directory
            gDirectory->mkdir("Ghost");
            gDirectory->cd("Ghost");
            {
              CreateParameterHistograms(hPartParam[3], hPartParam2D[3], 0, iPart, drawZR);
            }
            gDirectory->cd(".."); // Parameters
            gDirectory->mkdir("MCSignal");
            gDirectory->cd("MCSignal");
            {
              CreateParameterHistograms(hPartParam[6], hPartParam2D[6], 0, iPart, drawZR);
            }
            gDirectory->cd(".."); // Parameters
            
            
            bool plotPrimaryHistograms = abs(fParteff.partPDG[iPart]) == 310 ||
                                         abs(fParteff.partPDG[iPart]) == 3122 ||
                                         abs(fParteff.partPDG[iPart]) == 22 ||
                                         abs(fParteff.partPDG[iPart]) == 111 ||
                                         abs(fParteff.partPDG[iPart]) == 3312 ||
                                         abs(fParteff.partPDG[iPart]) == 3334;  
                                          
            bool plotSecondaryHistograms = abs(fParteff.partPDG[iPart]) == 310 ||
                                           abs(fParteff.partPDG[iPart]) == 3122 ||
                                           abs(fParteff.partPDG[iPart]) == 22 ||
                                           abs(fParteff.partPDG[iPart]) == 111;
                                            
            if(fStorePrimSecHistograms && plotPrimaryHistograms)
            {
              gDirectory->mkdir("Primary");
              gDirectory->cd("Primary");
              {
                CreateParameterSubfolder("NoConstraint (1C-Fit)", hPartParamPrimary, hPartParam2DPrimary, hFitQANoConstraint, iPart, true);
                CreateParameterSubfolder("MassConstraint (2C-Fit)", hPartParamPrimaryMass, hPartParam2DPrimaryMass, hFitQAMassConstraint, iPart, true);
                CreateParameterSubfolder("PVConstraint (3C-Fit)", hPartParamPrimaryTopo, hPartParam2DPrimaryTopo, hFitQATopoConstraint, iPart, true);
                CreateParameterSubfolder("PVMassConstraint (4C-Fit)", hPartParamPrimaryTopoMass, hPartParam2DPrimaryTopoMass, hFitQATopoMassConstraint, iPart, true);
              }
              gDirectory->cd(".."); // particle directory / Parameters
            }
            
            if(fStorePrimSecHistograms && plotSecondaryHistograms)
            {
              gDirectory->mkdir("Secondary");
              gDirectory->cd("Secondary");
              {
                CreateParameterSubfolder("NoConstraint (1C-Fit)", hPartParamSecondary, hPartParam2DSecondary, 0, iPart, true);
                CreateParameterSubfolder("MassConstraint (2C-Fit)", hPartParamSecondaryMass, hPartParam2DSecondaryMass, 0, iPart, true);
              }
              gDirectory->cd(".."); // particle directory / Parameters
            }
          }
        }
        gDirectory->cd(".."); //particle directory
      }
      gDirectory->cd(".."); //Particles
    }
    gDirectory->cd(".."); //main
    gDirectory->mkdir("PrimaryVertexQA");
    gDirectory->cd("PrimaryVertexQA");
    {
      struct {
        TString name;
        TString title;
        Int_t n;
        Double_t l,r;
      } Table[nHistosPV]=
      {
        {"PVResX",  "x_{rec}-x_{mc}, cm", 100, -0.1f, 0.1f},
        {"PVResY",  "y_{rec}-y_{mc}, cm", 100, -0.1f, 0.1f},
        {"PVResZ",  "z_{rec}-z_{mc}, cm", 100, -1.f, 1.f},
        {"PVPullX", "Pull X",             100, -6.f, 6.f},
        {"PVPullY", "Pull Y",             100, -6.f, 6.f},
        {"PVPullZ", "Pull Z",             100, -6.f, 6.f},
        {"Lost",    "Lost tracks",        102, -0.01f, 1.01f}        
      };
      
      TString parName[nHistosPVParam] = {"x","y","z","r","Ntracks","Chi2","NDF","Chi2NDF","prob", "PVpurity", 
                                         "ghostTr", "triggerTr", "pileupTr", "bgTr", "dzSamePV"};
      TString parAxisName[nHistosPVParam] = {"x [cm]","y [cm]","z [cm]","r [cm]","N tracks","Chi2","NDF","Chi2NDF","prob","purity",
                                             "ghost tracks [%]", "trigger tracks [%]", "pileup tracks [%]", "bg tracks [%]", "dz [cm]"};
      int nBins[nHistosPVParam] = {1000,1000,1000,1000,1001,10000,1001,10000,100,102,102,102,102,102,1000};
      float xMin[nHistosPVParam] = {-10., -10., -160.,  0,   -0.5,    0.,   -0.5,    0., 0., -0.01, -0.01, -0.01, -0.01, -0.01, 0.};
      float xMax[nHistosPVParam] = { 10.,  10.,  230., 10, 1000.5, 1000., 1000.5, 1000., 1.,  1.01,  1.01,  1.01,  1.01,  1.01, 100.};
      
      TString parName2D[nHistosPVParam2D] = {"xy"};
      TString parXAxisName2D[nHistosPVParam2D] = {"x [cm]"};
      TString parYAxisName2D[nHistosPVParam2D] = {"y [cm]"};
      int nBinsX2D[nHistosPVParam2D] = {1000};
      float xMin2D[nHistosPVParam2D] = {-1.};
      float xMax2D[nHistosPVParam2D] = { 1.};
      int nBinsY2D[nHistosPVParam2D] = {1000};
      float yMin2D[nHistosPVParam2D] = {-1.};
      float yMax2D[nHistosPVParam2D] = { 1.};
      
      for(int iH=0; iH<nHistosPVParam; iH++)
      {
        hPVParam[iH]       = new TH1F(parName[iH].Data(),(GetDirectoryPath()+parName[iH]).Data(),
                                        nBins[iH],xMin[iH],xMax[iH]);
        hPVParam[iH]->GetXaxis()->SetTitle(parAxisName[iH].Data());
      }

      for(int iH=0; iH<nHistosPVParam2D; iH++)
      {
        hPVParam2D[iH]       = new TH2F(parName2D[iH].Data(),(GetDirectoryPath()+parName2D[iH]).Data(),
                                        nBinsX2D[iH],xMin2D[iH],xMax2D[iH],
                                        nBinsY2D[iH],yMin2D[iH],yMax2D[iH]);
        hPVParam2D[iH]->GetXaxis()->SetTitle(parXAxisName2D[iH].Data());
        hPVParam2D[iH]->GetYaxis()->SetTitle(parYAxisName2D[iH].Data());
        hPVParam2D[iH]->GetYaxis()->SetTitleOffset(1.0);
      }
      
      gDirectory->mkdir("Efficiency");
      gDirectory->cd("Efficiency");
      {
        TString effName[nPVefficiency] = {"effVsNMCPVTracks","effVsNMCPV","effVsNMCTracks","effVsNPVTracks","effVsNPV","effVsNTracks"};
        int nBinsEff[nPVefficiency]  = { 100 , 100 ,  100 ,  100 , 100 , 1000 };
        float xMinEff[nPVefficiency] = {   0.,   0.,    0.,    0.,   0.,    0.};
        float xMaxEff[nPVefficiency] = { 100., 100., 1000.,  100., 100., 1000.};

        gDirectory->mkdir("Signal");
        gDirectory->cd("Signal");
        {
          for( int iH=0; iH<nPVefficiency; iH++ ){
            hPVefficiency[0][iH]   = new TProfile( effName[iH].Data(), (GetDirectoryPath()+effName[iH]).Data(), nBinsEff[iH], xMinEff[iH], xMaxEff[iH]);
          }
        }
        gDirectory->cd(".."); //L1

        gDirectory->mkdir("Pileup");
        gDirectory->cd("Pileup");
        {
          for( int iH=0; iH<nPVefficiency; iH++ ){
            hPVefficiency[1][iH]   = new TProfile( effName[iH].Data(), (GetDirectoryPath()+effName[iH].Data()), nBinsEff[iH], xMinEff[iH], xMaxEff[iH]);
          }
        }
        gDirectory->cd(".."); //L1
        
        gDirectory->mkdir("Signal_MCReconstructable");
        gDirectory->cd("Signal_MCReconstructable");
        {
          for( int iH=0; iH<nPVefficiency; iH++ ){
            hPVefficiency[2][iH]   = new TProfile( effName[iH].Data(), (GetDirectoryPath()+effName[iH].Data()), nBinsEff[iH], xMinEff[iH], xMaxEff[iH]);
          }
        }
        gDirectory->cd(".."); //L1
        
        gDirectory->mkdir("Pileup_MCReconstructable");
        gDirectory->cd("Pileup_MCReconstructable");
        {
          for( int iH=0; iH<nPVefficiency; iH++ ){
            hPVefficiency[3][iH]   = new TProfile( effName[iH].Data(), (GetDirectoryPath()+effName[iH].Data()), nBinsEff[iH], xMinEff[iH], xMaxEff[iH]);
          }
        }
        gDirectory->cd(".."); //L1
      }      
      gDirectory->cd(".."); //L1
      
      gDirectory->mkdir("PVTracksQA");
      gDirectory->cd("PVTracksQA");
      {
        TString resTrPV = "resTrPV";
        TString pullTrPV = "pullTrPV";
        TString parNameTrPV[nFitPVTracksQA/2] = {"X","Y","Z","Px","Py","Pz"};
        int nBinsTrPV = 100;
        float xMaxTrPV[nFitPVTracksQA/2] = {0.5,0.5,0.5,0.05,0.05,0.05};

        for( int iH=0; iH<nFitPVTracksQA/2; iH++ ){
          hFitPVTracksQA[iH]   = new TH1F((resTrPV+parNameTrPV[iH]).Data(),
                                                  (GetDirectoryPath()+resTrPV+parNameTrPV[iH]).Data(), 
                                                  nBinsTrPV, -xMaxTrPV[iH],xMaxTrPV[iH]);
          hFitPVTracksQA[iH+nFitPVTracksQA/2] = new TH1F((pullTrPV+parNameTrPV[iH]).Data(),
                                                  (GetDirectoryPath()+pullTrPV+parNameTrPV[iH]).Data(), 
                                                  nBinsTrPV, -6,6);
        }
      }      
      gDirectory->cd(".."); //L1
      
      gDirectory->mkdir("Signal");
      gDirectory->cd("Signal");
      {
        gDirectory->mkdir("FitQA");
        gDirectory->cd("FitQA");
        {
          gDirectory->mkdir("FitQAvcNMCPVTracks");
          gDirectory->cd("FitQAvcNMCPVTracks");
          {
            for(int iHPV=0; iHPV<nHistosPV-1; ++iHPV){
              hPVFitQa2D[0][0][iHPV] = new TH2F(Table[iHPV].name.Data(),(GetDirectoryPath()+Table[iHPV].title).Data(),
                                                500, 0., 5000.,
                                                Table[iHPV].n, Table[iHPV].l, Table[iHPV].r);
            }
          }
          gDirectory->cd(".."); //FitQA

          gDirectory->mkdir("FitQAvcNPVTracks");
          gDirectory->cd("FitQAvcNPVTracks");
          {
            for(int iHPV=0; iHPV<nHistosPV-1; ++iHPV){
              hPVFitQa2D[0][1][iHPV] = new TH2F(Table[iHPV].name.Data(),(GetDirectoryPath()+Table[iHPV].title).Data(),
                                                500, 0., 5000.,
                                                Table[iHPV].n, Table[iHPV].l, Table[iHPV].r);
            }
          }
          gDirectory->cd(".."); //FitQA
          
          for(int iHPV=0; iHPV<nHistosPV; ++iHPV){
            hPVFitQa[0][iHPV] = new TH1F(Table[iHPV].name.Data(),(GetDirectoryPath()+Table[iHPV].title).Data(),
                                         Table[iHPV].n, Table[iHPV].l, Table[iHPV].r);
          }
        }
        gDirectory->cd(".."); //Signal

        for(int iH=0; iH<nHistosPVParam; iH++)
        {
          hPVParamSignal[iH] = new TH1F((parName[iH]).Data(),(GetDirectoryPath()+parName[iH]).Data(),
                                        nBins[iH],xMin[iH],xMax[iH]);
          hPVParamSignal[iH]->GetXaxis()->SetTitle(parAxisName[iH].Data());
        }
      }      
      gDirectory->cd(".."); //L1

      gDirectory->mkdir("Pileup");
      gDirectory->cd("Pileup");
      {
        gDirectory->mkdir("FitQA");
        gDirectory->cd("FitQA");
        {
          gDirectory->mkdir("FitQAvcNMCPVTracks");
          gDirectory->cd("FitQAvcNMCPVTracks");
          {
            for(int iHPV=0; iHPV<nHistosPV-1; ++iHPV){
              hPVFitQa2D[1][0][iHPV] = new TH2F(Table[iHPV].name.Data(),(GetDirectoryPath()+Table[iHPV].title).Data(),
                                                500, 0., 5000.,
                                                Table[iHPV].n, Table[iHPV].l, Table[iHPV].r);
            }
          }
          gDirectory->cd(".."); //FitQA

          gDirectory->mkdir("FitQAvcNPVTracks");
          gDirectory->cd("FitQAvcNPVTracks");
          {
            for(int iHPV=0; iHPV<nHistosPV-1; ++iHPV){
              hPVFitQa2D[1][1][iHPV] = new TH2F(Table[iHPV].name.Data(),(GetDirectoryPath()+Table[iHPV].title).Data(),
                                                500, 0., 5000.,
                                                Table[iHPV].n, Table[iHPV].l, Table[iHPV].r);
            }
          }
          gDirectory->cd(".."); //FitQA
          
          for(int iHPV=0; iHPV<nHistosPV; ++iHPV){
            hPVFitQa[1][iHPV] = new TH1F(Table[iHPV].name.Data(),(GetDirectoryPath()+Table[iHPV].title).Data(),
                                         Table[iHPV].n, Table[iHPV].l, Table[iHPV].r);
          }
        }
        gDirectory->cd(".."); //Signal
        
        for(int iH=0; iH<nHistosPVParam; iH++)
        {
          hPVParamPileup[iH] = new TH1F((parName[iH]).Data(),(GetDirectoryPath()+parName[iH]).Data(),
                                        nBins[iH],xMin[iH],xMax[iH]);
          hPVParamPileup[iH]->GetXaxis()->SetTitle(parAxisName[iH].Data());
        }
      }      
      gDirectory->cd(".."); //L1
      
      gDirectory->mkdir("Background");
      gDirectory->cd("Background");
      {
        for(int iH=0; iH<nHistosPVParam; iH++)
        {
          hPVParamBG[iH] = new TH1F((parName[iH]).Data(),(GetDirectoryPath()+parName[iH]).Data(),
                                        nBins[iH],xMin[iH],xMax[iH]);
          hPVParamBG[iH]->GetXaxis()->SetTitle(parAxisName[iH].Data());
        }
      }      
      gDirectory->cd(".."); //L1
      
      gDirectory->mkdir("Ghost");
      gDirectory->cd("Ghost");
      {
        for(int iH=0; iH<nHistosPVParam; iH++)
        {
          hPVParamGhost[iH] = new TH1F((parName[iH]).Data(),(GetDirectoryPath()+parName[iH]).Data(),
                                        nBins[iH],xMin[iH],xMax[iH]);
          hPVParamGhost[iH]->GetXaxis()->SetTitle(parAxisName[iH].Data());
        }
      }
      gDirectory->cd(".."); //L1
    }
    gDirectory->cd(".."); //L1
    gDirectory->mkdir("TrackParameters");
    gDirectory->cd("TrackParameters");
    {
      TString chi2Name = "Chi2Prim";
      for(int iPart=0; iPart < KFPartEfficiencies::nParticles; iPart++)
      {
        TString chi2NamePart = "Chi2Prim";
        chi2NamePart += "_";
        chi2NamePart += fParteff.partName[iPart].data();
        hTrackParameters[iPart] = new TH1F(chi2NamePart.Data(), (GetDirectoryPath()+chi2NamePart).Data(), 1000, 0, 100);

      }
      hTrackParameters[KFPartEfficiencies::nParticles  ] = new TH1F("Chi2Prim_total", (GetDirectoryPath()+TString("Chi2Prim_total")), 1000, 0, 100);
      hTrackParameters[KFPartEfficiencies::nParticles+1] = new TH1F("Chi2Prim_prim", (GetDirectoryPath()+TString("Chi2Prim_prim")), 1000, 0, 100);
      hTrackParameters[KFPartEfficiencies::nParticles+2] = new TH1F("Chi2Prim_sec", (GetDirectoryPath()+TString("Chi2Prim_sec")), 1000, 0, 100);
      hTrackParameters[KFPartEfficiencies::nParticles+3] = new TH1F("Chi2Prim_ghost", (GetDirectoryPath()+TString("Chi2Prim_ghost")), 1000, 0, 100);
      
      hTrackParameters[KFPartEfficiencies::nParticles+4] = new TH1F("ProbPrim_total", (GetDirectoryPath()+TString("ProbPrim_total")), 10000, 0, 1);
      hTrackParameters[KFPartEfficiencies::nParticles+5] = new TH1F("ProbPrim_prim", (GetDirectoryPath()+TString("ProbPrim_prim")), 10000, 0, 1);
      hTrackParameters[KFPartEfficiencies::nParticles+6] = new TH1F("ProbPrim_sec", (GetDirectoryPath()+TString("ProbPrim_sec")), 10000, 0, 1);
      hTrackParameters[KFPartEfficiencies::nParticles+7] = new TH1F("ProbPrim_ghost", (GetDirectoryPath()+TString("ProbPrim_ghost")), 10000, 0, 1);
    }
    gDirectory->cd(".."); //particle directory
    curdir->cd();    
  }
}

void KFParticlePerformanceBase::CreateFitHistograms(TH1F* histo[nFitQA], int iPart)
{
  /** Creates 1D histograms with fit QA for decay with "iPart" number.
   ** \param[in,out] histo - array with pointers, for which the memory is allocated
   ** \param[in] iPart - number of the decay in the KF Particle Finder reconstruction scheme
   **/
  TString res = "res";
  TString pull = "pull";
  
  TString AxisNameResidual[nFitQA/2];
  TString AxisNamePull[nFitQA/2];

  AxisNameResidual[0] = "Residual (x^{reco} - x^{mc}) [cm]";
  AxisNameResidual[1] = "Residual (y^{reco} - y^{mc}) [cm]";
  AxisNameResidual[2] = "Residual (z^{reco} - z^{mc}) [cm]";
  AxisNameResidual[3] = "Residual (P_{x}^{reco} - P_{x}^{mc}) [GeV/c]";
  AxisNameResidual[4] = "Residual (P_{y}^{reco} - P_{y}^{mc}) [GeV/c]";
  AxisNameResidual[5] = "Residual (P_{z}^{reco} - P_{z}^{mc}) [GeV/c]";
  AxisNameResidual[6] = "Residual (E^{reco} - E^{mc}) [GeV/c^{2}]";
  AxisNameResidual[7] = "Residual (M^{reco} - M^{mc}) [GeV/c^{2}]"; 

  AxisNamePull[0] = "Pull x";
  AxisNamePull[1] = "Pull y";
  AxisNamePull[2] = "Pull z";
  AxisNamePull[3] = "Pull P_{x}";
  AxisNamePull[4] = "Pull P_{y}";
  AxisNamePull[5] = "Pull P_{z}";
  AxisNamePull[6] = "Pull E";
  AxisNamePull[7] = "Pull M";
  
  gDirectory->mkdir("FitQA");
  gDirectory->cd("FitQA");
  {
    TString parName[nFitQA/2] = {"X","Y","Z","Px","Py","Pz","E","M"};
    int nBins = 50;
//     float xMax[nFitQA/2] = {0.15,0.15,1.2,0.02,0.02,0.15,0.15,0.006};
    float xMax[nFitQA/2] = {5.f,5.f,10.f,0.1,0.1,0.5,0.5,0.006};
    float mult[nFitQA/2]={1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f};
    if(iPart>63 && iPart<75)
      for(int iMult=3; iMult<nFitQA/2; iMult++)
        mult[iMult] = 3;
    if(iPart>45 && iPart<64)
    {
#ifdef CBM
      for(int iMult=0; iMult<3; iMult++)
        mult[iMult] = 0.03;
      for(int iMult=3; iMult<nFitQA/2; iMult++)
        mult[iMult] = 3;
#else
      mult[2] = 0.1;
      for(int iMult=3; iMult<nFitQA/2; iMult++)
        mult[iMult] = 10;
      mult[5] = 2;
      mult[6] = 2;
#endif
    }
    if(iPart==44 || iPart==45)
    {
      mult[0] = 0.25;
      mult[1] = 0.5;
      mult[2] = 0.15;
      for(int iMult=3; iMult<nFitQA/2; iMult++)
        mult[iMult] = 4;
    }
    
    for( int iH=0; iH<nFitQA/2; iH++ )
    {
      histo[iH]   = new TH1F((res+parName[iH]).Data(),
                             (GetDirectoryPath()+res+parName[iH]).Data(), 
                             nBins, -mult[iH]*xMax[iH],mult[iH]*xMax[iH]);
      histo[iH]->GetXaxis()->SetTitle(AxisNameResidual[iH].Data());
      histo[iH+8] = new TH1F((pull+parName[iH]).Data(),
                             (GetDirectoryPath()+pull+parName[iH]).Data(), 
                             nBins, -6,6);
      histo[iH+8]->GetXaxis()->SetTitle(AxisNamePull[iH+8].Data());
    }
  }
  gDirectory->cd("..");
}

void KFParticlePerformanceBase::CreateEfficiencyHistograms(TProfile* histo[3][nPartEfficiency], TProfile2D* histo2[3][nPartEfficiency2D], THnSparseF* histoN[4])
{
  /** Creates efficiency plots in the current ROOT folder.
   ** \param[in,out] histo - 1D efficiency plots
   ** \param[in,out] histo2 - 2D efficiency plots
   ** \param[in,out] histoN - ND efficiency plots
   **/
  gDirectory->mkdir("Efficiency");
  gDirectory->cd("Efficiency");
  {//vs p, pt, y, z, c*tau, decay length, l, r
    TString partNameEff[nPartEfficiency] = {"EffVsP","EffVsPt","EffVsY","EffVsZ","EffVsCT","EffVsDL","EffVsL","EffVsR","EffVsMt" };
    TString partAxisNameEff[nPartEfficiency] = {"p [GeV/c]","p_{t} [GeV/c]",
                                                "y", "z [cm]", "Life time c#tau [cm]", "Decay length [cm]", 
                                                "L [cm]", "Rxy [cm]", "m_{t} [GeV/c^{2}]"};
#ifdef CBM
    int nBinsEff[nPartEfficiency]  = { 100 , 100 , 40 ,  360 ,  60 ,  60 ,  140 ,  60 , 100 };
    float xMinEff[nPartEfficiency] = {   0.,   0.,  0.,  -10., -10., -10.,    0.,   0. ,  0.};
    float xMaxEff[nPartEfficiency] = {  20.,   5.,  4.,   80.,  50.,  50.,   70.,  30. ,  4.};
#else
    int nBinsEff[nPartEfficiency]  = { 100 , 100 ,  30 ,   100 ,   60 ,   60 ,  100 ,  100 , 100  };
    float xMinEff[nPartEfficiency] = {   0.,   0.,  -3.,  -230.,  -10.,  -10.,    0.,    0.,   0. };
    float xMaxEff[nPartEfficiency] = {  10.,  10.,   0.,   230.,   50.,   50.,   50.,   50.,  10. };
#endif
    TString effTypeName[3] = {"All particles",
                              "Reconstructable daughters",
                              "Reconstructed daughters"};
    
    for(int iEff=0; iEff<3; iEff++)
    {
      gDirectory->mkdir(effTypeName[iEff].Data());
      gDirectory->cd(effTypeName[iEff].Data());
      {
        for(int iH=0; iH<nPartEfficiency; iH++)
        {
          histo[iEff][iH] = new TProfile( partNameEff[iH].Data(), (GetDirectoryPath()+partAxisNameEff[iH]).Data(), nBinsEff[iH], xMinEff[iH], xMaxEff[iH]);
          histo[iEff][iH]->GetYaxis()->SetTitle("Efficiency");                  
          histo[iEff][iH]->GetYaxis()->SetTitleOffset(1.0);  
          histo[iEff][iH]->GetXaxis()->SetTitle(partAxisNameEff[iH].Data());
        }
        
        histo2[iEff][0] = new TProfile2D( "EffVsPtVsY", (GetDirectoryPath()+partAxisNameEff[2]+partAxisNameEff[1]).Data(), 
                                          nBinsEff[2], xMinEff[2], xMaxEff[2], nBinsEff[1], xMinEff[1], xMaxEff[1]);
        histo2[iEff][0]->GetZaxis()->SetTitle("Efficiency");
        histo2[iEff][0]->GetXaxis()->SetTitle(partAxisNameEff[2].Data());
        histo2[iEff][0]->GetYaxis()->SetTitle(partAxisNameEff[1].Data());
        histo2[iEff][0]->GetYaxis()->SetTitleOffset(1.0);
        
        histo2[iEff][1] = new TProfile2D( "EffVsMtVsY", (GetDirectoryPath()+partAxisNameEff[2]+partAxisNameEff[8]).Data(), 
                                          nBinsEff[2], xMinEff[2], xMaxEff[2], nBinsEff[8], xMinEff[8], xMaxEff[8]);
        histo2[iEff][1]->GetZaxis()->SetTitle("Efficiency");
        histo2[iEff][1]->GetXaxis()->SetTitle(partAxisNameEff[2].Data());
        histo2[iEff][1]->GetYaxis()->SetTitle(partAxisNameEff[8].Data());
        histo2[iEff][1]->GetYaxis()->SetTitleOffset(1.0);
      }
      gDirectory->cd("..");// particle directory / Efficiency
    }

    if(fStore3DEfficiency)
    {
      gDirectory->mkdir("Multidimensional");
      gDirectory->cd("Multidimensional");
      {
#if 1
        constexpr int nDimensions = 3;
        //                                phi           theta    p
        int nBins[nDimensions]  {         100,             25,  25 };
        double xMin[nDimensions]{-TMath::Pi(),  TMath::Pi()/2,   0.};
        double xMax[nDimensions]{ TMath::Pi(),    TMath::Pi(),  10.};
        
        histoN[0] = new THnSparseF( "Acceptance_reco", GetDirectoryPath()+"Acceptance reco", 
                                    nDimensions, nBins, xMin, xMax);
        histoN[0] -> GetAxis(0) -> SetTitle("#phi");
        histoN[0] -> GetAxis(1) -> SetTitle("#theta");
        histoN[0] -> GetAxis(2) -> SetTitle("p [GeV/c]");
        
        gDirectory->Append(histoN[0]);
        
        histoN[1] = new THnSparseF( "Acceptance_mc", GetDirectoryPath()+"Acceptance mc", 
                                    nDimensions, nBins, xMin, xMax);
        histoN[1] -> GetAxis(0) -> SetTitle("#phi");
        histoN[1] -> GetAxis(1) -> SetTitle("#theta");
        histoN[1] -> GetAxis(2) -> SetTitle("p [GeV/c]");
        
        gDirectory->Append(histoN[1]);

        constexpr int nDimensions2 = 4;
        //                                   phi           theta    p   l/dl
        int nBins2[nDimensions2]  {          100,             25,  25 ,  20 };
        double xMin2[nDimensions2]{ -TMath::Pi(),  TMath::Pi()/2,   0.,   0.};
        double xMax2[nDimensions2]{  TMath::Pi(),    TMath::Pi(),  10., 100.};

        histoN[2] = new THnSparseF( "Cuts_reco", GetDirectoryPath()+"Cuts reco", 
                                    nDimensions2, nBins2, xMin2, xMax2);
        histoN[2] -> GetAxis(0) -> SetTitle("p_{t} [GeV/c]");
        histoN[2] -> GetAxis(1) -> SetTitle("p [GeV/c]");
        histoN[2] -> GetAxis(2) -> SetTitle("Decay length [cm]");
        
        gDirectory->Append(histoN[2]);
        
        histoN[3] = new THnSparseF( "Cuts_mc", GetDirectoryPath()+"Cuts mc", 
                                    nDimensions2, nBins2, xMin2, xMax2);
        histoN[3] -> GetAxis(0) -> SetTitle("p_{t} [GeV/c]");
        histoN[3] -> GetAxis(1) -> SetTitle("p [GeV/c]");
        histoN[3] -> GetAxis(2) -> SetTitle("Decay length [cm]");
        
        gDirectory->Append(histoN[3]);
#else
        constexpr int nDimensions = 5;
        //                                phi         theta    p      pt    eta
        int nBins[nDimensions]  {         100,         100 , 100 ,   100 ,   30 };
        double xMin[nDimensions]{-TMath::Pi(),           0.,   0.,     0.,   -3.};
        double xMax[nDimensions]{ TMath::Pi(),  TMath::Pi(),  10.,    10.,    0.};
        
        histoN[0] = new THnSparseF( "Acceptance_reco", GetDirectoryPath()+"phi-theta-p-pt-eta reco", 
                                    nDimensions, nBins, xMin, xMax);
        histoN[0] -> GetAxis(0) -> SetTitle("#phi");
        histoN[0] -> GetAxis(1) -> SetTitle("#theta");
        histoN[0] -> GetAxis(2) -> SetTitle("p [GeV/c]");
        histoN[0] -> GetAxis(3) -> SetTitle("p_{t} [GeV/c]");
        histoN[0] -> GetAxis(4) -> SetTitle("eta");
        
        gDirectory->Append(histoN[0]);
        
        histoN[1] = new THnSparseF( "Acceptance_mc", GetDirectoryPath()+"phi-theta-p-pt-eta mc", 
                                    nDimensions, nBins, xMin, xMax);
        histoN[1] -> GetAxis(0) -> SetTitle("#phi");
        histoN[1] -> GetAxis(1) -> SetTitle("#theta");
        histoN[1] -> GetAxis(2) -> SetTitle("p [GeV/c]");
        histoN[1] -> GetAxis(3) -> SetTitle("p_{t} [GeV/c]");
        histoN[1] -> GetAxis(4) -> SetTitle("eta");
        
        gDirectory->Append(histoN[1]);


        
        constexpr int nDimensions2 = 5;
        //                           pt      y    p    DL,  ldl
        int nBins2[nDimensions2]  {   50,   30, 100,  100 , 100 };
        double xMin2[nDimensions2]{    0,   -3,   0.,   0.,   0.};
        double xMax2[nDimensions2]{    5,    0,  10., 100., 100.};

        histoN[2] = new THnSparseF( "Cuts_reco", GetDirectoryPath()+"pt-y-p-DL-ldl reco", 
                                    nDimensions2, nBins2, xMin2, xMax2);
        histoN[2] -> GetAxis(0) -> SetTitle("p_t [GeV/c]");
        histoN[2] -> GetAxis(1) -> SetTitle("y");
        histoN[2] -> GetAxis(2) -> SetTitle("p [GeV/c]");
        histoN[2] -> GetAxis(3) -> SetTitle("decay length [cm]");
        histoN[2] -> GetAxis(4) -> SetTitle("l/dl");
        
        gDirectory->Append(histoN[2]);
        
        histoN[3] = new THnSparseF( "Cuts_mc", GetDirectoryPath()+"pt-y-p-DL-ldl mc", 
                                    nDimensions2, nBins2, xMin2, xMax2);
        histoN[3] -> GetAxis(0) -> SetTitle("p_t [GeV/c]");
        histoN[3] -> GetAxis(1) -> SetTitle("y");
        histoN[3] -> GetAxis(2) -> SetTitle("p [GeV/c]");
        histoN[3] -> GetAxis(3) -> SetTitle("decay length [cm]");
        histoN[3] -> GetAxis(4) -> SetTitle("l/dl");
        
        gDirectory->Append(histoN[3]);
#endif
      }
      gDirectory->cd("..");
    }
  }
  gDirectory->cd("..");// particle directory
}

void KFParticlePerformanceBase::CreateParameterHistograms(TH1F* histoParameters[KFPartEfficiencies::nParticles][nHistoPartParam],
                                                          TH2F *histoParameters2D[KFPartEfficiencies::nParticles][nHistoPartParam2D],
                                                          TH3F *histoParameters3D[KFPartEfficiencies::nParticles][nHistoPartParam3D],
                                                          int iPart, bool drawZR)
{
  /** Creates histograms with parameter distributions for decay with "iPart" number.
   ** \param[in,out] histoParameters - 1D histograms
   ** \param[in,out] histoParameters2D - 2D histograms
   ** \param[in,out] histoParameters3D - 3D histograms
   ** \param[in] iPart - number of the decay in the KF Particle Finder reconstruction scheme
   ** \param[in] drawZR - flag showing if Z-R histogram should be created
   **/
  TString parName[nHistoPartParam] = {"M","p","p_{t}","y","DecayL","c#tau","chi2ndf","prob","#theta","phi",
                                      "X","Y","Z","R", "L", "l/dl","m_{t}","Multiplicity",
                                      "dX", "dY", "dZ", "dPx", "dPy", "dPz", "dE", "dM"};
  TString parTitle[nHistoPartParam];
  TString parName2D[nHistoPartParam2D] = {"y-p_{t}", "Z-R", "Armenteros", "y-m_{t}"};
  TString parTitle2D[nHistoPartParam2D];
  TString parName3D[nHistoPartParam3D] = {"y-p_{t}-M", "y-m_{t}-M", "centrality-pt-M", "centrality-y-M", "centrality-mt-M", "ct-pt-M", "dalitz", "dalitz2","dalitz3","dalitzM2", "dalitz2M2", "dalitz3M2"};
  TString parTitle3D[nHistoPartParam3D];
  for(int iParam=0; iParam<nHistoPartParam; iParam++)
  {
    TString path = GetDirectoryPath();
    parTitle[iParam] = path + parName[iParam];
    if(iParam<nHistoPartParam2D)
      parTitle2D[iParam] = path + parName2D[iParam];
    if(iParam<nHistoPartParam3D)
      parTitle3D[iParam] = path + parName3D[iParam];
  }
  
  TString parAxisName[nHistoPartParam] = {"m [GeV/c^{2}]","p [GeV/c]","p_{t} [GeV/c]",
                                          "y","Decay length [cm]","Life time c#tau [cm]",
                                          "chi2/ndf","prob","#theta [rad]",
                                          "phi [rad]","x [cm]","y [cm]","z [cm]","Rxy [cm]", "L [cm]", "L/dL","m_{t} [GeV/c^{2}]","Multiplicity",
                                          "#sigma_{X}, [cm]", "#sigma_{Y}, [cm]", "#sigma_{Z}, [cm]",
                                          "#sigma_{p_{x}}/|p_{x}|", "#sigma_{p_{y}}/|p_{y}|", "#sigma_{p_{z}}/|p_{z}|",
                                          "#sigma_{E}/E", "#sigma_{M}, [GeV/c^{2}]"};
#ifdef CBM
  int nBins[nHistoPartParam] =  {1000, // M
                                  100, // p
                                  100, // pt
                                   40, // y
                                   60, // DecayL
                                   60, // ctau
                                  100, // chi2/ndf
                                  100, // prob
                                  100, // theta
                                  100, // phi
                                  200, // X
                                  200, // Y
                                  360, // Z
                                   60, // R
                                  140, // L
                                  200, // L/dL
                                  100, // Mt
                                  fParteff.partMaxMult[iPart]+1};
  float xMin[nHistoPartParam] = { fParteff.partMHistoMin[iPart], // M
                                  0.f, // p
                                  0.f, // pt
                                  0.f, // y
                                -10.f, // DecayL
                                -10.f, // ctau
                                  0.f, // chi2/ndf
                                  0.f, // prob
                                  0.f, // theta
                             -3.1416f, // phi
                                -50.f, // X
                                -50.f, // Y
                                -10.f, // Z
                                  0.f, // R
                                  0.f, // L
                                 -1.f, // L/dL
                                  0.f, // Mt
                                 -0.5f };
  float xMax[nHistoPartParam] = { fParteff.partMHistoMax[iPart], // M
                                  20.f, // p
                                   5.f, // pt
                                   4.f, // y
                                  50.f, // DecayL
                                  50.f, // ctau
                                  20.f, // chi2/ndf
                                   1.f, // prob
                               3.1416f, // theta
                               3.1416f, // phi
                                  50.f, // X
                                  50.f, // Y
                                  80.f, // Z
                                  30.f, // R
                                  70.f, // L
                                  35.f, // L/dL
                                  4.f, // Mt
                                  float(fParteff.partMaxMult[iPart])+0.5f};
#else
  int nBins[nHistoPartParam] = {1000, // M
                                 100, // p
                                 100, // pt
                                  30, // y
                                  60, // DecayL
                                  60, // ctau
                                 100, // chi2/ndf
                                 100, // prob
                                 100, // theta
                                 100, // phi
                                 100, // X
                                 100, // Y
                                 100, // Z
                                 500, // R
                                 500, // L
                                1000, // L/dL
                                 100, // Mt
                                 fParteff.partMaxMult[iPart]+1,
                                 100, // dX
                                 100, // dY
                                 100, // dZ
                                 100, // dPx
                                 100, // dPy
                                 100, // dPz
                                 100, // dE
                                 100  // dM
  };
  float xMin[nHistoPartParam] = { fParteff.partMHistoMin[iPart], // M
                                  0.f, // p
                                  0.f, // pt
                                -3.0f, // y
                                -10.f, // DecayL
                                -10.f, // ctau
                                  0.f, // chi2/ndf
                                  0.f, // prob
                                  0.f, // theta
                         -TMath::Pi(), // phi
                                -10.f, // X
                                -10.f, // Y
                               -230.f, // Z
                                  0.f, // R
                                  0.f, // L
                                 -1.f, // L/dL
                                  0.f, // Mt
                                 -0.5f,
                                  0.f, // dX
                                  0.f, // dY
                                  0.f, // dZ
                                  0.f, // dPx
                                  0.f, // dPy
                                  0.f, // dPz
                                  0.f, // dE
                                  0.f  // dM
  };
  float xMax[nHistoPartParam] = { fParteff.partMHistoMax[iPart], // M
                                  10.f, // p
                                  10.f, // pt
                                  0.f, // y
                                  50.f, // DecayL
                                  50.f, // ctau
                                  20.f, // chi2/ndf
                                   1.f, // prob
                           TMath::Pi(), // theta
                           TMath::Pi(), // phi
                                  10.f, // X
                                  10.f, // Y
                                 230.f, // Z
                                 200.f, // R
                                 200.f, // L
                                 200.f, // L/dL
                                  10.f, // Mt
                                  float(fParteff.partMaxMult[iPart])+0.5f,
                                  10.f, // dX
                                  10.f, // dY
                                  10.f, // dZ
                                  0.2f, // dPx
                                  0.2f, // dPy
                                  0.2f, // dPz
                                  0.2f, // dE
                                 0.02f  // dM
  };
  if(iPart < 9)
  {
    xMin[10] =-50; xMin[11] =-50; xMin[12] =-100;
    xMax[10] = 50; xMax[11] = 50; xMax[12] = 100; xMax[13] = 50; xMax[14] = 50; 
  }
#endif
  for(int iH=0; iH<nHistoPartParam; iH++)
  {
    histoParameters[iPart][iH] = new TH1F(parName[iH].Data(),parTitle[iH].Data(),
                                          nBins[iH],xMin[iH],xMax[iH]);
    histoParameters[iPart][iH]->GetXaxis()->SetTitle(parAxisName[iH].Data());
  }

  histoParameters2D[iPart][0] = new TH2F(parName2D[0].Data(),parTitle2D[0].Data(),
                                    nBins[3],xMin[3],xMax[3],
                                    nBins[2],xMin[2],xMax[2]);
  histoParameters2D[iPart][0]->GetXaxis()->SetTitle("y");
  histoParameters2D[iPart][0]->GetYaxis()->SetTitle("p_{t} [GeV/c]");
  histoParameters2D[iPart][0]->GetYaxis()->SetTitleOffset(1.0);

  if(drawZR)
  {
    histoParameters2D[iPart][1] = new TH2F(parName2D[1].Data(),parTitle2D[1].Data(),
                                      nBins[12],xMin[12],xMax[12],
                                      nBins[13],xMin[13],xMax[13]);
    histoParameters2D[iPart][1]->GetXaxis()->SetTitle("Z [cm]");
    histoParameters2D[iPart][1]->GetYaxis()->SetTitle("R [cm]");
    histoParameters2D[iPart][1]->GetYaxis()->SetTitleOffset(1.0);
  }
  else
    histoParameters2D[iPart][1] = NULL;
  
  //create armenteros plot
  if(IsCollectArmenteros(iPart))
  {
    histoParameters2D[iPart][2] = new TH2F(parName2D[2].Data(),parTitle2D[2].Data(),
                                           50, -1.f, 1.f,
                                          150,  0.f, 1.f);
    histoParameters2D[iPart][2]->GetXaxis()->SetTitle("#alpha (p_{L}^{+}-p_{L}^{-})/(p_{L}^{+}+p_{L}^{-})");
    histoParameters2D[iPart][2]->GetYaxis()->SetTitle("q_{t} [GeV/c]");
    histoParameters2D[iPart][2]->GetYaxis()->SetTitleOffset(1.0);
  }
  else
    histoParameters2D[iPart][2] = NULL;
  //create y-mt plot
  histoParameters2D[iPart][3] = new TH2F(parName2D[3].Data(),parTitle2D[3].Data(),
                                         nBins[3],xMin[3], xMax[3],     //y
                                         nBins[16],xMin[16],xMax[16]); //Mt
  histoParameters2D[iPart][3]->GetXaxis()->SetTitle("y");
  histoParameters2D[iPart][3]->GetYaxis()->SetTitle("m_{t} [GeV/c]");
  histoParameters2D[iPart][3]->GetYaxis()->SetTitleOffset(1.0);
  
  
  if( histoParameters3D && IsCollect3DHistogram(iPart) )
  {
    histoParameters3D[iPart][0] = new TH3F(parName3D[0].Data(),parTitle3D[0].Data(),
                                      nBins[3],xMin[3],xMax[3],
                                      nBins[2],xMin[2],xMax[2],
                                      nBins[0],xMin[0],xMax[0]);
    histoParameters3D[iPart][0]->GetXaxis()->SetTitle("y");
    histoParameters3D[iPart][0]->GetYaxis()->SetTitle("p_{t} [GeV/c]");
    histoParameters3D[iPart][0]->GetYaxis()->SetTitleOffset(1.0);
    histoParameters3D[iPart][0]->GetZaxis()->SetTitle("M");
    
    histoParameters3D[iPart][1] = new TH3F(parName3D[1].Data(),parTitle3D[1].Data(),
                                           nBins[3],xMin[3],xMax[3],
                                           nBins[16],xMin[16],xMax[16],
                                           nBins[0],xMin[0],xMax[0]);
    histoParameters3D[iPart][1]->GetXaxis()->SetTitle("y");
    histoParameters3D[iPart][1]->GetYaxis()->SetTitle("m_{t} [GeV/c]");
    histoParameters3D[iPart][1]->GetYaxis()->SetTitleOffset(1.0);
    histoParameters3D[iPart][1]->GetZaxis()->SetTitle("M");
    
    int centralityHisto[3] = {2,3,16};
    for(int iCH = 0; iCH<3; iCH++)
    {
      histoParameters3D[iPart][2+iCH] = new TH3F(parName3D[2+iCH].Data(),parTitle3D[2+iCH].Data(),
                                                 10,0.,10.,
                                                 nBins[centralityHisto[iCH]],xMin[centralityHisto[iCH]],xMax[centralityHisto[iCH]],
                                                 nBins[0],xMin[0],xMax[0]);
      histoParameters3D[iPart][2+iCH]->GetXaxis()->SetTitle("centrality bin");
      histoParameters3D[iPart][2+iCH]->GetYaxis()->SetTitle(parAxisName[centralityHisto[iCH]]);
      histoParameters3D[iPart][2+iCH]->GetYaxis()->SetTitleOffset(1.0);
      histoParameters3D[iPart][2+iCH]->GetZaxis()->SetTitle("M");
    }
    
    histoParameters3D[iPart][5] = new TH3F(parName3D[5].Data(),parTitle3D[5].Data(),
                                           nBins[5],xMin[5],xMax[5],
                                           nBins[2],xMin[2],xMax[2],
                                           nBins[0],xMin[0],xMax[0]);
    histoParameters3D[iPart][5]->GetXaxis()->SetTitle("c#tau [cm]");
    histoParameters3D[iPart][5]->GetYaxis()->SetTitle("p_{t} [GeV/c]");
    histoParameters3D[iPart][5]->GetYaxis()->SetTitleOffset(1.0);
    histoParameters3D[iPart][5]->GetZaxis()->SetTitle("M");
    
    if(IsCollectDalitz(iPart))
    {
      int nBinsM12  = 100;
      int nBinsM23  = 100;
      int nBinsM13  = 100;
      int nBinsMass = 100;
      
      const int pdg1 = fParteff.partDaughterPdg[iPart][0];
      const int pdg2 = fParteff.partDaughterPdg[iPart][1];
      const int pdg3 = fParteff.partDaughterPdg[iPart][2];
      
      int index1 = fParteff.GetParticleIndex(pdg1);
      int index2 = fParteff.GetParticleIndex(pdg2);
      int index3 = fParteff.GetParticleIndex(pdg3);
      
      float m1 = fParteff.partMass[index1];
      float m2 = fParteff.partMass[index2];
      float m3 = fParteff.partMass[index3];
      
      const float MMin = m1+m2+m3;
      const float MMax = fParteff.partMass[iPart] + 0.05;
      
      const float m12Min = (m1+m2)*0.998;
      const float m12Max = (MMax-m3);

      const float m23Min = (m2+m3)*0.998;
      const float m23Max = (MMax-m1);      
      
      const float m13Min = (m1+m3)*0.998;
      const float m13Max = (MMax-m2);  
      
      TString axis12 = "m_{"; 
      axis12 += fParteff.partName[index1];
      axis12 += fParteff.partName[index2];
      axis12 += "}";

      TString axis23 = "m_{"; 
      axis23 += fParteff.partName[index2];
      axis23 += fParteff.partName[index3];
      axis23 += "}";

      TString axis13 = "m_{"; 
      axis13 += fParteff.partName[index1];
      axis13 += fParteff.partName[index3];
      axis13 += "}";

      histoParameters3D[iPart][6] = new TH3F(parName3D[6].Data(),parTitle3D[6].Data(),
                                             nBinsM12, m12Min, m12Max,
                                             nBinsM23, m23Min, m23Max,
                                             nBinsMass,MMin,   MMax);      
      histoParameters3D[iPart][6]->GetXaxis()->SetTitle(axis12 + " [GeV/c]");
      histoParameters3D[iPart][6]->GetYaxis()->SetTitle(axis23 + " [GeV/c]");
      histoParameters3D[iPart][6]->GetYaxis()->SetTitleOffset(1.0);
      histoParameters3D[iPart][6]->GetZaxis()->SetTitle("M");
      
      histoParameters3D[iPart][7] = new TH3F(parName3D[7].Data(),parTitle3D[7].Data(),
                                             nBinsM13, m13Min, m13Max,
                                             nBinsM23, m23Min, m23Max,
                                             nBinsMass,MMin,   MMax);
      histoParameters3D[iPart][7]->GetXaxis()->SetTitle(axis13 + " [GeV/c]");
      histoParameters3D[iPart][7]->GetYaxis()->SetTitle(axis23 + " [GeV/c]");
      histoParameters3D[iPart][7]->GetYaxis()->SetTitleOffset(1.0);
      histoParameters3D[iPart][7]->GetZaxis()->SetTitle("M");

      histoParameters3D[iPart][8] = new TH3F(parName3D[8].Data(),parTitle3D[8].Data(),
                                             nBinsM12, m12Min, m12Max,
                                             nBinsM13, m13Min, m13Max,
                                             nBinsMass,MMin,   MMax);
      histoParameters3D[iPart][8]->GetXaxis()->SetTitle(axis12 + " [GeV/c]");
      histoParameters3D[iPart][8]->GetYaxis()->SetTitle(axis13 + " [GeV/c]");
      histoParameters3D[iPart][8]->GetYaxis()->SetTitleOffset(1.0);
      histoParameters3D[iPart][8]->GetZaxis()->SetTitle("M");
      
      histoParameters3D[iPart][9] = new TH3F(parName3D[9].Data(),parTitle3D[9].Data(),
                                             nBinsM12, m12Min*m12Min, m12Max*m12Max,
                                             nBinsM23, m23Min*m23Min, m23Max*m23Max,
                                             nBinsMass,MMin,   MMax);      
      histoParameters3D[iPart][9]->GetXaxis()->SetTitle(axis12 + "^{2} [GeV^{2}/c^{2}]");
      histoParameters3D[iPart][9]->GetYaxis()->SetTitle(axis23 + "^{2} [GeV^{2}/c^{2}]");
      histoParameters3D[iPart][9]->GetYaxis()->SetTitleOffset(1.0);
      histoParameters3D[iPart][9]->GetZaxis()->SetTitle("M");
      
      histoParameters3D[iPart][10] = new TH3F(parName3D[10].Data(),parTitle3D[10].Data(),
                                              nBinsM13, m13Min*m13Min, m13Max*m13Max,
                                              nBinsM23, m23Min*m23Min, m23Max*m23Max,
                                              nBinsMass,MMin,   MMax);
      histoParameters3D[iPart][10]->GetXaxis()->SetTitle(axis13 + "^{2} [GeV^{2}/c^{2}]");
      histoParameters3D[iPart][10]->GetYaxis()->SetTitle(axis23 + "^{2} [GeV^{2}/c^{2}]");
      histoParameters3D[iPart][10]->GetYaxis()->SetTitleOffset(1.0);
      histoParameters3D[iPart][10]->GetZaxis()->SetTitle("M");
      
      histoParameters3D[iPart][11] = new TH3F(parName3D[11].Data(),parTitle3D[11].Data(),
                                              nBinsM12, m12Min*m12Min, m12Max*m12Max,
                                              nBinsM13, m13Min*m13Min, m13Max*m13Max,
                                              nBinsMass,MMin,   MMax);
      histoParameters3D[iPart][11]->GetXaxis()->SetTitle(axis12 + "^{2} [GeV^{2}/c^{2}]");
      histoParameters3D[iPart][11]->GetYaxis()->SetTitle(axis13 + "^{2} [GeV^{2}/c^{2}]");
      histoParameters3D[iPart][11]->GetYaxis()->SetTitleOffset(1.0);
      histoParameters3D[iPart][11]->GetZaxis()->SetTitle("M");
#if 0
      TString nameXAxis = "#sqrt{3}*(T_{";
      nameXAxis += fParteff.partName[index1];
      nameXAxis += "} - T_{";
      nameXAxis += fParteff.partName[index2];
      nameXAxis += "})/Q";
      TString nameYAxis = "(2T_{";
      nameYAxis += fParteff.partName[index3];
      nameYAxis += "} - T_{";
      nameYAxis += fParteff.partName[index1];
      nameYAxis += "} - T_{";
      nameYAxis += fParteff.partName[index2];
      nameYAxis += "})/Q";
      histoParameters3D[iPart][10] = new TH3F(parName3D[10].Data(),parTitle3D[10].Data(),
                                              100, -1.2, 1.2,
                                              100, -1.2, 1.2,
                                              nBinsMass,MMin,   MMax);
      histoParameters3D[iPart][10]->GetXaxis()->SetTitle(nameXAxis);
      histoParameters3D[iPart][10]->GetYaxis()->SetTitle(nameYAxis);
      histoParameters3D[iPart][10]->GetYaxis()->SetTitleOffset(1.0);
      histoParameters3D[iPart][10]->GetZaxis()->SetTitle("M");
#endif
    }
  }
  else if(histoParameters3D)
  {
    histoParameters3D[iPart][0] = NULL;
    histoParameters3D[iPart][1] = NULL;
    for(int iCH = 0; iCH<3; iCH++)
      histoParameters3D[iPart][2+iCH] = NULL;
    histoParameters3D[iPart][5] = NULL;
  }
}

bool KFParticlePerformanceBase::IsCollectZRHistogram(int iParticle) const
{
  /** Checks if Z-R histogram for decay "iParticle" should be created. */
  return (abs(fParteff.partPDG[iParticle]) == 310 ||
          abs(fParteff.partPDG[iParticle]) == 3122 ||
          abs(fParteff.partPDG[iParticle]) == 3312 ||
          abs(fParteff.partPDG[iParticle]) == 3334 ||
          abs(fParteff.partPDG[iParticle]) == 22) && 
          fStoreMCHistograms && fStoreZRHistograms && (!fStore3DEfficiency);
}

bool KFParticlePerformanceBase::IsCollect3DHistogram(int iParticle) const
{
  /** Checks if 3D histograms for decay "iParticle" should be created. */
  return 0; //TODO
//   return (!fStore3DEfficiency) && (abs(fParteff.partPDG[iParticle]) == 310 ||
//          abs(fParteff.partPDG[iParticle]) == 3122 ||
//          abs(fParteff.partPDG[iParticle]) == 3312 ||
//          abs(fParteff.partPDG[iParticle]) == 3334 ||
//          (abs(fParteff.partPDG[iParticle]) >= 3000 && 
//           abs(fParteff.partPDG[iParticle]) <= 3029) ||
//          abs(fParteff.partPDG[iParticle]) == 3103 ||
//          abs(fParteff.partPDG[iParticle]) == 3203 ||
// #ifdef CBM
//          abs(fParteff.partPDG[iParticle]) == 7003112 ||
//          abs(fParteff.partPDG[iParticle]) == 7003222 ||
//          abs(fParteff.partPDG[iParticle]) == 7003312 ||
//          abs(fParteff.partPDG[iParticle]) == 8003222 ||
//          abs(fParteff.partPDG[iParticle]) == 9000321);
// #else
//          abs(fParteff.partPDG[iParticle]) == 421 ||
//          abs(fParteff.partPDG[iParticle]) == 429 ||
//          abs(fParteff.partPDG[iParticle]) == 426 ||
//          abs(fParteff.partPDG[iParticle]) == 411 ||
//          abs(fParteff.partPDG[iParticle]) == 431 ||
//          abs(fParteff.partPDG[iParticle]) == 4122 ||
//          abs(fParteff.partPDG[iParticle]) == 521 ||
//          abs(fParteff.partPDG[iParticle]) == 511);
// #endif
}

bool KFParticlePerformanceBase::IsCollectArmenteros(int iParticle) const
{
  /** Checks if Armenteros-Podoliansky plot for decay "iParticle" should be created. */
  return 0; //TODO
//   return (!fStore3DEfficiency) && (abs(fParteff.partPDG[iParticle]) == 310 ||
//          abs(fParteff.partPDG[iParticle]) == 3122 ||
//          abs(fParteff.partPDG[iParticle]) == 3312 ||
//          abs(fParteff.partPDG[iParticle]) == 3334 ||
//          abs(fParteff.partPDG[iParticle]) == 22 ||
//          abs(fParteff.partPDG[iParticle]) == 111 ||
//          abs(fParteff.partPDG[iParticle]) == 3003 ||
//          abs(fParteff.partPDG[iParticle]) == 3103 ||
//          abs(fParteff.partPDG[iParticle]) == 3004 ||
//          abs(fParteff.partPDG[iParticle]) == 3005 ||
//          abs(fParteff.partPDG[iParticle]) == 3203 ||
//          abs(fParteff.partPDG[iParticle]) == 3008 ||
//          abs(fParteff.partPDG[iParticle]) == 3000 ||
//          abs(fParteff.partPDG[iParticle]) == 333 ||
// #ifdef CBM
//          abs(fParteff.partPDG[iParticle]) == 7003112 ||
//          abs(fParteff.partPDG[iParticle]) == 7003222 ||
//          abs(fParteff.partPDG[iParticle]) == 7003312 ||
//          abs(fParteff.partPDG[iParticle]) == 8003222 ||
//          abs(fParteff.partPDG[iParticle]) == 9000321);
// #else
//          abs(fParteff.partPDG[iParticle]) == 421 ||
//          abs(fParteff.partPDG[iParticle]) == 420 ||
//          abs(fParteff.partPDG[iParticle]) == 426 ||
//          abs(fParteff.partPDG[iParticle]) == 521 ||
//          abs(fParteff.partPDG[iParticle]) == 511);        
// #endif
}

bool KFParticlePerformanceBase::IsCollectDalitz(int iParticle) const
{
  /** Checks if Armenteros-Podoliansky plot for decay "iParticle" should be created. */
  return 0; // TODO
//   return fParteff.partPDG[iParticle] == 3006 ||
//          fParteff.partPDG[iParticle] == 3007 ||
//          fParteff.partPDG[iParticle] == 3012 ||
//          fParteff.partPDG[iParticle] == 3013 ||
//          
//          fParteff.partPDG[iParticle] == 3014 ||
//          fParteff.partPDG[iParticle] == 3015 ||
//          fParteff.partPDG[iParticle] == 3017 ||
//          fParteff.partPDG[iParticle] == 3018 ||
//          fParteff.partPDG[iParticle] == 3020 ||
//          fParteff.partPDG[iParticle] == 3021 ||
//          fParteff.partPDG[iParticle] == 3023 ||
//          fParteff.partPDG[iParticle] == 3024 ||
//          fParteff.partPDG[iParticle] == 3026 ||
//          fParteff.partPDG[iParticle] == 3027 ||
//          
//          fParteff.partPDG[iParticle] == 3028 || 
//          fParteff.partPDG[iParticle] == 3029 ;
}

void KFParticlePerformanceBase::CreateParameterSubfolder(TString folderName, 
                                                         TH1F* histoParameters[nParametersSet][KFPartEfficiencies::nParticles][nHistoPartParam],
                                                         TH2F* histoParameters2D[nParametersSet][KFPartEfficiencies::nParticles][nHistoPartParam2D],
                                                         TH1F* histoFit[KFPartEfficiencies::nParticles][nFitQA], int iPart, bool withWrongPVHypothesis)
{
  /** Creates all subfolders in the current ROOT directory for the current decay channel. */
  gDirectory->mkdir(folderName.Data());
  gDirectory->cd(folderName.Data());
  {
    gDirectory->mkdir("Signal");
    gDirectory->cd("Signal");
    {
      CreateParameterHistograms(histoParameters[1], histoParameters2D[1], 0, iPart);
    }
    gDirectory->cd("..");
    if(withWrongPVHypothesis)
    {
      gDirectory->mkdir("WrongPVHypothesis");
      gDirectory->cd("WrongPVHypothesis");
      {
        CreateParameterHistograms(histoParameters[4], histoParameters2D[4], 0, iPart);
      }
      gDirectory->cd("..");
    }
    gDirectory->mkdir("Background");
    gDirectory->cd("Background");
    {
      CreateParameterHistograms(histoParameters[2], histoParameters2D[2], 0, iPart);
    }
    gDirectory->cd("..");
    gDirectory->mkdir("Ghost");
    gDirectory->cd("Ghost");
    {
      CreateParameterHistograms(histoParameters[3], histoParameters2D[3], 0, iPart);
    }
    gDirectory->cd("..");
    
    CreateParameterHistograms(histoParameters[0], histoParameters2D[0], 0, iPart);
    if(histoFit!=0)
      CreateFitHistograms(histoFit[iPart], iPart);
  }
  gDirectory->cd("..");
}

TString KFParticlePerformanceBase::GetDirectoryPath()
{
  /** Returns the path to the current folder. It is used as an addition to the histogram name. */
  TString path = gDirectory->GetPath();
  int fileNamePosition = path.Index("Finder/");
  path.Remove(0, fileNamePosition+7);
  path.ReplaceAll("Particles/", "");
  path.ReplaceAll("/Parameters", "");
  path+=" ";
  return path;
}

#endif //DO_TPCCATRACKER_EFF_PERFORMANCE

