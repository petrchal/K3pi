


 //the declarations must be here, not global 0,0 -automatic range and #bins, otherwise 100 bins
  TPlotDefinitions Event_plots_FXT{
    TPlotDef{"Evt.Vz","Prim. vtx position","Vz[cm]",0,500,"Vz"},
    TPlotDef{"Evt.eventId","eventId","eventId",0,10000000},
    TPlotDef{"Evt.vzVpd","Vpd vtx position","VZ_vpd[cm]",-0,500,"VPDdif,Vz"},
    TPlotDef{"Evt.Vz-Evt.vzVpd","VpdVz-Vz","diff Vz[cm]",-8,8,"VPDdif"},
    TPlotDef{"Evt.ZDCx","ZDC coincidence rate","f[hz]",0,800,"ZDC"}, 
    TPlotDef{"Evt.BBCx","BBC coincidence rate","f[Hz]",0,300000,"BBC"},
    TPlotDef{"Evt.nBTOFMatch","Num of TOF matched track","#TOF matched",0,500,"TOFmatch"},
    TPlotDef{"Evt.refMult","RefMult ","RefMult",0,1000,"refMult"},
    TPlotDef{"Evt.gRefMult","gRefMult","gRefMult",0,1000,"gRefMult"},
    TPlotDef{"Evt.runId","runId","runId",0,0},
    TPlotDef{"Evt.triggerIds","triggerIds","Trigger Ids",0,1000000},
    TPlotDef{"Evt.nK3piP","number of 3pi+ vertexes","nK3piP",-2,10},
    TPlotDef{"Evt.nK3piN","number of 3pi- vertexes","nK3piN",-2,10},
    };

    TPlotDefinitions_2D Event_plots_2D_FXT{
    TPlotDef_2D{"Evt.Vz","Evt.vzVpd", "Prim. vtx Vz TPC vs VPD ","Vz TPC[cm]","Vz VPD[cm]",0,300,0,500},
     //NOTE: for 2D histograms one cannot use automatic range finding
     TPlotDef_2D{"Evt.runId","Evt.eventId", "run vs event Id ","runId","eventid",0,3.*1e7,0,10000000},
     //TPlotDef_2D{"Evt.ZDCx","Evt.BBCx", "ZDC vs BBC rates","ZDC coinc. rate [Hz]","BBC coinc. rate [Hz]",0,200,0,200000},
     //TPlotDef_2D{"Evt.Vz","Evt.BBCx", "event Vz vs BBC rates","Vz TPC [cm]","BBC coinc. rate [Hz]",-100,100,0,200000},
     //TPlotDef_2D{"Evt.Vz","Evt.ZDCx", "event Vz vs ZDZ rates","Vz TPC [cm]","ZDC coinc. rate [Hz]",-100,100,0,200},
     TPlotDef_2D{"Evt.refMult","Evt.gRefMult", "refMult vs gRefMult","refMult","gRefMult",0,100,0,100},
     TPlotDef_2D{"Evt.refMult","Evt.nBTOFMatch", "refMult vs nBTOFMatch","refMult","nBTOFMatch",0,100,0,200},
     TPlotDef_2D{"Evt.gRefMult","Evt.nBTOFMatch", "gRefMult vs nBTOFMatch","gRefMult","nBTOFMatch",0,100,0,200},
     TPlotDef_2D{"Evt.nBTOFMatch","Evt.nK3piP","nBTOFMatch vs number of 3pi+ vertexes","nBTOFMatch","nK3piP",0,200,-2,10},
     TPlotDef_2D{"Evt.nBTOFMatch","Evt.nK3piN","nBTOFMatch vs number of 3pi- vertexes","nBTOFMatch","nK3piN",0,200,-2,10},
     TPlotDef_2D{"Evt.nK3piP","Evt.nK3piN","number of 3pi+ vs 3p- vertexes","nK3piP","nK3piN",-2,10,-2,10},
};

 TPlotDefinitions RecoVtx_plots{
    //mother plots 
    //variable to plot , name , ranges,  comma-separated list of cuts which are disabled before plotting
    
    TPlotDef{"mother_m","mass of 3pi vertex","m[GeV/c]",0.4,.6,"minv"}, 
    TPlotDef{"mother_PID","3pi vertex PID","PID",0,0},
    TPlotDef{"mother_isMc","Is 3pi MC","",-5,5},
    TPlotDef{"matchedKF","matched by KF","",-5,5},
    TPlotDef{"matchedGeom","matched by dp","",-5,5},
    TPlotDef{"mother_chi2ndf","3piVtx_chi2ndf","xi/ndf",0,0,"cut_mom_ch2ndf"},
    TPlotDef{"mother_PV_chi2","3piVtx_PV_chi2","xi/ndf",0,150}, //note in data are values ..-600
    TPlotDef{"mother_PV_l","3piVtx_PV_l","l[cm]",0,15},
    TPlotDef{"mother_PV_dl","3piVtx_PV_dl","\sigma l[cm]",0,10},
    TPlotDef{"(mother_PV_dl>0)?mother_PV_l/mother_PV_dl:0","3piVtx l per dl","\sigma l/dl",0,15},
  
    TPlotDef{"mother_eta_PVX","eta from 3pi vertex at PVX","eta",-2.,2.,"eta"},
    TPlotDef{"decay_Vr","radial position of decay","Vr[cm]",50,200,"decay_Vr"}, 
    TPlotDef{"decay_Vz","z position of decay","Vz[cm]",-100,200}, 
    TPlotDef{"mother_pt_PVX","pt from 3pi vertex at PVX","pt[GeV/c^2]",0,1.5,"pt"},
    TPlotDef{"mother_phi_PVX","phi from 3pi vertex at PVX","phi",-7,7.},
    TPlotDef{"MaxHits(decay_Vr)","maximum number of hits of matched kaon","MaxHits(decay_Vr)",-0.5,99.5,"decay_Vr"},
 

    TPlotDef{"d.qaTruth[0]+d.qaTruth[1]+d.qaTruth[2]","sum of daughter qaTruths","",0.,300.},
  
    //first decay daugter 
    TPlotDef{"d.index[0]","daughter[0] index","index",-0.5,6.5},
    TPlotDef{"d.nhits[0]","daughter[0] nHitsFit","nHitsFit",-0.5,99.5,"nhits_daughters"},
    TPlotDef{"d.nhits_pos[0]","daughter[0] NHits possible from DST","nHitsPos",-0.5,99.5,"nhits_possible"},
   // TPlotDef{"MaxHitsDaughter","daughter[0] NHitsf possible from Vr","lastPoint - MaxHits(decay_V)",-0.5,99.5,"nhits_possible"},
 

    //second decay daugter 
    //TPlotDef{"d.index[1]","daughter[1] index","index",-0.5,6.5},
    TPlotDef{"d.nhits[1]","daughter[1] nHitsFit","nHitsFit",-0.5,99.5,"nhits_daughters"},
    TPlotDef{"d.nhits_pos[1]","daughter[1] NHits possible from DST","nHitsPos",-0.5,99.5,"nhits_possible"},
    //TPlotDef{"MaxHitsDaughter[1]","daughter[1] NHitsf possible from Vr","MaxHits(Vr)",-0.5,99.5,"nhits_possible"},
 
    //thisrd decay daugter 
    //TPlotDef{"d.index[1]","daughter[1] index","index",-0.5,6.5},
    TPlotDef{"d.nhits[2]","daughter[2] nHitsFit","nHitsFit",-0.5,99.5,"nhits_daughters"},
    TPlotDef{"d.nhits_pos[2]","daughter[2] NHits possible from DST","nHitsPos",-0.5,99.5,"nhits_possible"},
    //TPlotDef{"MaxHitsDaughter[2]","daughter[2] NHitsf possible from Vr","MaxHits(Vr)",-0.5,99.5,"nhits_possible"},
 
 

/*
    TPlotDef{"(float)d.nhits[0]/(float)d.nhits_pos[0]","daughter NHits/possible","nHitsFit/nHitsPos",-0.5,2,"nhits_posrat"},
    TPlotDef{"(float)d.nhits[0]/(float)MaxHitsDaughter","daughter NHits/possible recalc","nHitsFit/nHitsPos recalculated",-0.5,2,"nhits_posrat"},
    TPlotDef{"d.PvtxDca_official[0]","daughter Prim. vtx DCA","dca[cm]",-1,100},
    TPlotDef{"d.PvtxDcaXY_official[0]","daughter Prim. vtx DCA_XY","dca_xy[cm]",-1,10},
    TPlotDef{"d.PvtxDcaZ_official[0]","daughter Prim. vtx DCA_Z","dca_z[cm]",-10,10},
    TPlotDef{"d.lastPointR[0]","daughter last hit position","r[cm]",80,200},
    TPlotDef{"d.match_chi2[0]","daughter match_chi2","chi2",0,5000},
    //TPlotDef{"d.idTruth[0]","daughter idTruth","idTruth",0,0,},
    TPlotDef{"d.qaTruth[0]","daughter qaTruth","qaTruth",0,150},
    */
 }; //3piVtx_plot

TPlotDefinitions_2D RecoVtx_plots_2D{
    //vertex 
    TPlotDef_2D{"mother_chi2ndf", "mother_m","3pi vertex chi2 vs mass","xi/ndf of found 3pi","m[GeV/c]",-1,11,0.4,0.6,"minv,3piVtx_chi"},
    TPlotDef_2D{"mother_PV_chi2", "mother_m","3pi chi2 to PV vs mass","xi to PV of found  3pi","m[GeV/c]",-1,11,0.4,0.6,"minv,3piVtx_chi"},
   
    TPlotDef_2D{"decay_Vz", "decay_Vr","3pi - reconstructed position Z vs radius","z[cm]","r[cm]",-500,500.,0.,200,"decay_Vr,eta,pt"},
    TPlotDef_2D{"decay_Vx", "decay_Vy","3pi - reconstructed position X-Y","x[cm]","y[cm]",-200,200.,-200,200,"decay_Vr,eta,pt"},
    TPlotDef_2D{"mother_eta_PVX", "mother_pt_PVX","3pi - pt vs eta","eta","pt[GeV/c^2]",-3.,3.,0,2.,"eta,pt"},
    TPlotDef_2D{"decay_Vr", "mother_pt_PVX","3pi - pt vs r","r[cm]","pt[GeV/c^2]",0.,200.,0,2.,"eta,pt,decay_Vr"},
    TPlotDef_2D{"mother_PV_l", "mother_PV_dl","3pi - distance to PV vs error","l to pV[cm]","sigma l to PV[cm]",0,15.,-0,10},
    TPlotDef_2D{"mother_PV_l", "(mother_PV_dl>0)?mother_PV_l/mother_PV_dl:0","3pi - distance to PV vs nsigma","l to pV[cm]","nsigma l to PV",0,15.,-0,10},
     
    //daughters[0]
     TPlotDef_2D{"d.nhits_pos[0]", "d.nhits[0]","daughter[0] - nhits possible vs measured","nhits possible","hnits",-0.5,99.5,-0.5,99.5,"decay_Vr"},
     TPlotDef_2D{"d.nhits_pos[0]", "decay_Vr","daughter[0] - nhits possible vs r","nhits possible","r[cm]",-0.5,99.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits[0]",     "decay_Vr","daughter[0] - nhits measured vs r","nhits","r[cm]",-0.5,99.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits_pos[0]", "d.lastPointR[0]","daughter[0] - nhits possible vs lastPointR","nhits possible","lastPoint[cm]",0,200.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits[0]",     "d.lastPointR[0]","daughter[0] - nhits measured vs lastPointR","nhits","lastPoint[cm]",0,200.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.lastPointR[0]", "decay_Vr","daughter[0] - lastPointR vs r","lastPoint[cm]","r[cm]",0.,200.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits_pos[0]", "MaxHits(decay_Vr)","daughter[0] - nhits possible vs calculated max","nhits possible","MaxHits(Vr)",-0.5,99.5,-0.5,99.5,"decay_Vr"},
   
     TPlotDef_2D{"d.nhits_pos[1]", "d.nhits[1]","daughter[1] - nhits possible vs measured","nhits possible","hnits",-0.5,99.5,-0.5,99.5,"decay_Vr"},
     TPlotDef_2D{"d.nhits_pos[1]", "decay_Vr","daughter[1] - nhits possible vs r","nhits possible","r[cm]",-0.5,99.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits[1]", "decay_Vr","daughter[1] - nhits measured vs r","nhits","r[cm]",-0.5,99.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits_pos[1]", "d.lastPointR[1]","daughter[1] - nhits possible vs lastPointR","nhits possible","lastPoint[cm]",0,200.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits[1]", "d.lastPointR[1]","daughter[1] - nhits measured vs lastPointR","nhits","lastPoint[cm]",0,200.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.lastPointR[1]", "decay_Vr","daughter[1] - lastPointR vs r","lastPoint[cm]","r[cm]",0.,200.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits_pos[1]", "MaxHits(decay_Vr)","daughter[1] - nhits possible vs calculated max","nhits possible","MaxHits(Vr)",-0.5,99.5,-0.5,99.5,"decay_Vr"},
    
     TPlotDef_2D{"d.nhits_pos[2]", "d.nhits[2]","daughter[2] - nhits possible vs measured","nhits possible","hnits",-0.5,99.5,-0.5,99.5,"decay_Vr"},
     TPlotDef_2D{"d.nhits_pos[2]", "decay_Vr","daughter[2] - nhits possible vs r","nhits possible","r[cm]",-0.5,99.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits[2]", "decay_Vr","daughter[2] - nhits measured vs r","nhits","r[cm]",-0.5,99.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits_pos[2]", "d.lastPointR[2]","daughter[2] - nhits possible vs lastPointR","nhits possible","lastPoint[cm]",0,200.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits[2]", "d.lastPointR[2]","daughter[2] - nhits measured vs lastPointR","nhits","lastPoint[cm]",0,200.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.lastPointR[2]", "decay_Vr","daughter[2] - lastPointR vs r","lastPoint[cm]","r[cm]",0.,200.5,0.,200.,"decay_Vr"},
     TPlotDef_2D{"d.nhits_pos[2]", "MaxHits(decay_Vr)","daughter[2] - nhits possible vs calculated max","nhits possible","MaxHits(Vr)",-0.5,99.5,-0.5,99.5,"decay_Vr"},
   

    /*
    TPlotDef{"mother_eta_PVX","eta from 3pi vertex at PVX","eta",-2.,2.,"eta"},
    TPlotDef{"mother_pt_PVX","pt from 3pi vertex at PVX","pt[GeV/c^2]",0,1.5,"pt"},
    TPlotDef{"mother_phi_PVX","phi from 3pi vertex at PVX","phi",-7,7.},
    TPlotDef{"MaxHits(decay_Vr)","maximum number of hits of matched kaon","hnits",-0.5,99.5,"decay_Vr"},
 

    TPlotDef{"d.qaTruth[0]+d.qaTruth[1]+d.qaTruth[2]","sum of daughter qaTruths","",0.,300.},
  
    //first decay daugter 
  
    TPlotDef{"(float)d.nhits[0]/(float)d.nhits_pos[0]","daughter NHits/possible","nHitsFit/nHitsPos",-0.5,2,"nhits_posrat"},
    TPlotDef{"(float)d.nhits[0]/(float)MaxHitsDaughter","daughter NHits/possible recalc","nHitsFit/nHitsPos recalculated",-0.5,2,"nhits_posrat"},
    TPlotDef{"d.PvtxDca_official[0]","daughter Prim. vtx DCA","dca[cm]",-1,100},
    TPlotDef{"d.PvtxDcaXY_official[0]","daughter Prim. vtx DCA_XY","dca_xy[cm]",-1,10},
    TPlotDef{"d.PvtxDcaZ_official[0]","daughter Prim. vtx DCA_Z","dca_z[cm]",-10,10},
    TPlotDef{"d.lastPointR[0]","daughter last hit position","r[cm]",80,200},
    TPlotDef{"d.match_chi2[0]","daughter match_chi2","chi2",0,5000},
    //TPlotDef{"d.idTruth[0]","daughter idTruth","idTruth",0,0,},
    TPlotDef{"d.qaTruth[0]","daughter qaTruth","qaTruth",0,150},
    */
 }; //3piVtx_plot



//const int K_match=4;//[4]..is matched (via dp from Dst) ...d[3] via KFP
#define K_match 4

TPlotDefinitions MatchedKaon_plots{
    TPlotDef{"d.nhits[K_match]","kaon nHitsFit","nHitsFit",-0.5,99.5,"kaon_nhits"},
    TPlotDef{"d.nhits_pos[K_match]","kaon nHitsPos (from MuDst)","nHitsPos",-0.5,99.5,},
    TPlotDef{"d.nhits_dEdx[K_match]","kaon dEdx hist","nhit dEdx",-0.5,99.5},
    TPlotDef{"d.PvtxDca_official[K_match]","kaon: Prim. vtx DCA","dca[cm]",-1,10,"kaon_DCA"},
    //TPlotDef{"PvtxDca_corrected","corrected kaon Prim. vtx DCA","dca[cm]",-1,10,"kaon_DCA"},
    TPlotDef{"d.PvtxDcaXY_official[K_match]","kaon Prim. vtx DCA_XY","dca_xy[cm]",-3,3,"kaon_DCA"},
    //TPlotDef{"PvtxDcaXY_corrected","corrected kaon Prim. vtx DCA_XY","dca_xy[cm]",-3,3,"kaon_DCA"},
    //TPlotDef{"PvtxDcaXY_corrected-d.PvtxDcaXY_official[K_match]","difference corrected - kaon Prim. vtx DCA_XY","dca_xy[cm]",-3,3,"kaon_DCA"},
    TPlotDef{"d.PvtxDcaZ_official[K_match]","kaon Prim. vtx DCA_Z","dca_z[cm]",-3,3,"kaon_DCA"},
    TPlotDef{"d.PvtxDca_mu[K_match]","recalc from helix Prim. vtx DCA ","dca[cm]",-1,10,"kaon_DCA"},
    TPlotDef{"d.PvtxDcaXY_mu[K_match]","recalc from helix Prim. vtx DCA XY ","dca_xy[cm]",-3,3,"kaon_DCA"},
    TPlotDef{"d.lastPointR[K_match]","kaon: last hit position","r[cm]",0,200},
    TPlotDef{"d.match_chi2[K_match]","kaon: match_chi2","chi2",0,0},
    TPlotDef{"d.idTruth[K_match]","kaon: idTruth","idTruth",0,0},
    TPlotDef{"d.qaTruth[K_match]","kaon: qaTruth","qaTruth",0,150},
    TPlotDef{"d.isBest[K_match]","kaon: isBest","",-3,3,"K_best"},
    TPlotDef{"d.lastPointR[K_match]-decay_Vr","kaon: dR -radial difference betwen last hit and decay vtx","dR[cm]",-100,200,"K_lastR"},
    TPlotDef{"d.dp_Decay[K_match]","kaon: dP at decay","dp[Gev/c^2]",0,0.4,"K_dp"},
  //  TPlotDef{"d.dp_decay_KF[K_match]","kaon: dP at decay from KF","dp[Gev/c^2]",0,0.2},
    TPlotDef{"d.DecayDca_mu[K_match]","kaon: DCA at decay from StHelix","dca[cm]",0,3,"K_dp"},
    //TPlotDef{"d.DecayDcaXY_mu[K_match]","kaon: DCA_xy at decay from StHelix","dca_xy[cm]",0,3,"K_dp"},
 //   TPlotDef{"d.DecayDca_KF[K_match]","kaon: DCA at decay from KF track","dca",0,3},   
    TPlotDef{"(float)d.nhits[K_match]/(float)MaxHits(decay_Vr)","nhits/nposhits (calculated).","nhits/nposhits",-1,2,"kaon_hits_ratio"},
    TPlotDef{"(float)d.nhits[K_match]/(float)d.nhits_pos[K_match]","nhits/nposhits (from MUDst).","nhits/nposhits",-1,2,"kaon_hits_ratio"}

 };

 //TPlotDefinitions2D MatchedKaon_plots_2D{

  TPlotDefinitions FXT_Efficiency_plots{
    //3piVtx plots 
    //variable to plot , name , ranges,  comma-separated list of cuts which are disabled before plotting
    TPlotDef{"decay_Vr","radial decay position","r[cm]",50,200},
    TPlotDef{"mother_pt_PVX","matched kaon pt","pt[GeV/c^2]",0,1.2},
    TPlotDef{"mother_m","reconstructed mass","M_inv[GeV/c]",0.45,0.55},
    TPlotDef{"mother_eta_PVX","pseudorapidity","eta",-2,2},
    TPlotDef{"mother_phi_PVX","phi","phi",-4,4},
       //TPlotDef{"3piVtx_PV_l/3piVtx_PV_dl","nsig l","nsig l",0,15}, //bad idea
       //TPlotDef{"d.nhits[4]","number of hits","nhits",-0.5,99.5},   //bad idea
    TPlotDef{"Evt.Vz","event vertex z-position","Vz[cm]",-180,220},
    TPlotDef{"decay_Vz","decay z-position","Z_decay[cm]",-200,200},
    TPlotDef{"MaxHits(decay_Vr)","max number of fitted hist. (calculated)","max #hits",-0.5,99.5},
 };
