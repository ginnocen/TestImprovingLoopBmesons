//Working tag v0
//16 July 2014
//13PM 
//GMI

#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <iostream>
#include <TNtuple.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <cmath>
#include "loop.h"

#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677
#define KSHORT_MASS 0.497614
#define KSTAR_MASS  0.89594
#define PHI_MASS    1.019455
#define JPSI_MASS   3.096916

Double_t massmincut=5.;
Double_t massmaxcut=6.;
Double_t ptmincut=10.;
Float_t mumumassCut=0.15;

Float_t trkPtCutBplus=0.9;
Float_t chi2clCutBplus=1.32e-02;
Float_t d0d0ErrCutBplus=3.41;
Float_t cosdthetaCutBplus=-3.46e-01;

Float_t trkPtCutBzero=0.7;
Float_t chi2clCutBzero=1.65e-01;
Float_t d0d0ErrCutBzero=4.16;
Float_t cosdthetaCutBzero=7.50e-01;
Float_t invmasstktkBzero=2.33e-01;


void fillTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, int typesize, float track_mass1, float track_mass2, int REAL, int PbpMC)
{ 
  void TrackDirectGenLabel(int,int,int,int,int,int,int&,int&,int&);
  void ParticleResonanceGenLabel(int,int,int,int,int,int,int,int,int,int&,int&,int&);

  //Event Info
  Event = EvtInfo_EvtNo;
  Run = EvtInfo_RunNo+10*PbpMC;
  size = typesize+1;
  besttktkmass = 0;
  bestchi2 = 0;

  // Trigger Info
  HLT_PAL1DoubleMu0_v1 = Bfr_HLT_PAL1DoubleMu0_v1;
  HLT_PAL1DoubleMu0_v1_Prescl = Bfr_HLT_PAL1DoubleMu0_v1_Prescl;
  HLT_PADimuon0_NoVertexing_v1 = Bfr_HLT_PADimuon0_NoVertexing_v1;
  HLT_PADimuon0_NoVertexing_v1_Prescl = Bfr_HLT_PADimuon0_NoVertexing_v1_Prescl;
  HLT_PAL1DoubleMu0_HighQ_v1 = Bfr_HLT_PAL1DoubleMu0_HighQ_v1;
  HLT_PAL1DoubleMu0_HighQ_v1_Prescl = Bfr_HLT_PAL1DoubleMu0_HighQ_v1_Prescl;
  HLT_PAL1DoubleMuOpen_v1 = Bfr_HLT_PAL1DoubleMuOpen_v1;
  HLT_PAL1DoubleMuOpen_v1_Prescl = Bfr_HLT_PAL1DoubleMuOpen_v1_Prescl;
  HLT_PAL2DoubleMu3_v1 = Bfr_HLT_PAL2DoubleMu3_v1;
  HLT_PAL2DoubleMu3_v1_Prescl = Bfr_HLT_PAL2DoubleMu3_v1_Prescl;
  HLT_PAMu3_v1 = Bfr_HLT_PAMu3_v1;
  HLT_PAMu3_v1_Prescl = Bfr_HLT_PAMu3_v1_Prescl;
  HLT_PAMu7_v1 = Bfr_HLT_PAMu7_v1;
  HLT_PAMu7_v1_Prescl = Bfr_HLT_PAMu7_v1_Prescl;
  HLT_PAMu12_v1 = Bfr_HLT_PAMu12_v1;
  HLT_PAMu12_v1_Prescl = Bfr_HLT_PAMu12_v1_Prescl;

  //B Info
  bP->SetXYZ(BInfo_px[j],BInfo_py[j],BInfo_pz[j]*0);
  bVtx->SetXYZ(BInfo_vtxX[j]-EvtInfo_PVx,
	       BInfo_vtxY[j]-EvtInfo_PVy,
	       BInfo_vtxZ[j]*0-EvtInfo_PVz*0);
  b4P->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);

  bindex[typesize] = typesize;
  y[typesize] = b4P->Rapidity();
  dtheta[typesize] = bP->Angle(*bVtx);
  pt[typesize] = BInfo_pt[j];
  eta[typesize] = BInfo_eta[j];
  phi[typesize] = BInfo_phi[j];
  chi2cl[typesize] = TMath::Prob(BInfo_vtxchi2[j],BInfo_vtxdof[j]);
  d0[typesize] = sqrt((BInfo_vtxX[j]-EvtInfo_PVx)*(BInfo_vtxX[j]-EvtInfo_PVx)+(BInfo_vtxY[j]-EvtInfo_PVy)*(BInfo_vtxY[j]-EvtInfo_PVy));
  vx[typesize] = BInfo_vtxX[j] - EvtInfo_PVx;
  vy[typesize] = BInfo_vtxY[j] - EvtInfo_PVy;
  d0Err[typesize] = sqrt(BInfo_vtxXE[j]*BInfo_vtxXE[j]+BInfo_vtxYE[j]*BInfo_vtxYE[j]);
  mass[typesize] = BInfo_mass[j];
  tktkmass[typesize] = BInfo_tktk_mass[j];
  chi2ndf[typesize] = BInfo_vtxchi2[j]/BInfo_vtxdof[j];
  lxy[typesize] = ((BInfo_vtxX[j]-EvtInfo_PVx)*BInfo_px[j] + (BInfo_vtxY[j]-EvtInfo_PVy)*BInfo_py[j])/BInfo_pt[j];
  isbestchi2[typesize] = 0;
  isbesttktkmass[typesize] = 0;
  kstar[typesize] = 0;
  if(BInfo_type[j]==4) kstar[typesize]=1;
  else if(BInfo_type[j]==5) kstar[typesize]=2;
  
  float mu1px,mu1py,mu1pz,mu1E;
  float mu2px,mu2py,mu2pz,mu2E;

  mu1px = MuonInfo_pt[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]*cos(MuonInfo_phi[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]);
  mu1py = MuonInfo_pt[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]*sin(MuonInfo_phi[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]);
  mu1pz = MuonInfo_pt[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]*sinh(MuonInfo_eta[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]);
  b4P->SetXYZM(mu1px,mu1py,mu1pz,MUON_MASS);
  mu1E = b4P->E();
  mu1eta[typesize] = b4P->Eta();
  mu1phi[typesize] = b4P->Phi();
  mu1y[typesize] = b4P->Rapidity();
  mu1pt[typesize] = b4P->Pt();
  mu1p[typesize] = b4P->P();

  mu2px = MuonInfo_pt[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]*cos(MuonInfo_phi[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]);
  mu2py = MuonInfo_pt[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]*sin(MuonInfo_phi[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]);
  mu2pz = MuonInfo_pt[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]*sinh(MuonInfo_eta[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]);
  b4P->SetXYZM(mu2px,mu2py,mu2pz,MUON_MASS);
  mu2E = b4P->E();
  mu2eta[typesize] = b4P->Eta();
  mu2phi[typesize] = b4P->Phi();
  mu2y[typesize] = b4P->Rapidity();
  mu2pt[typesize] = b4P->Pt();
  mu2p[typesize] = b4P->P();

  b4P->SetPxPyPzE(mu1px+mu2px,
		  mu1py+mu2py,
		  mu1pz+mu2pz,
		  mu1E+mu2E);
  mumumass[typesize] = b4P->Mag();
  mumueta[typesize] = b4P->Eta();
  mumuphi[typesize] = b4P->Phi();
  mumuy[typesize] = b4P->Rapidity();
  mumupt[typesize] = b4P->Pt();


  //jpsi section
  ujmass[typesize] = BInfo_uj_mass[BInfo_rfuj_index[j]];
  ujvProb[typesize] = TMath::Prob(BInfo_uj_vtxchi2[BInfo_rfuj_index[j]],BInfo_uj_vtxdof[BInfo_rfuj_index[j]]);
  b4P->SetXYZM(BInfo_uj_px[BInfo_rfuj_index[j]],
	       BInfo_uj_py[BInfo_rfuj_index[j]],
	       BInfo_uj_pz[BInfo_rfuj_index[j]],
	       BInfo_uj_mass[BInfo_rfuj_index[j]]);
  ujpt[typesize] = b4P->Pt();
  ujeta[typesize] = b4P->PseudoRapidity();
  ujphi[typesize] = b4P->Phi();
  ujy[typesize] = b4P->Rapidity();
  ujlxy[typesize] = ((BInfo_uj_vtxX[BInfo_rfuj_index[j]]-EvtInfo_PVx)*BInfo_uj_px[BInfo_rfuj_index[j]] + (BInfo_uj_vtxY[BInfo_rfuj_index[j]]-EvtInfo_PVy)*BInfo_uj_py[BInfo_rfuj_index[j]])/ujpt[typesize];

  //track section
  float tk1px,tk1py,tk1pz,tk1E;
  float tk2px,tk2py,tk2pz,tk2E;

  if(BInfo_type[j]==1 || BInfo_type[j]==2)
    {
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],track_mass1);
      trk1Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
      trk1D0Err[typesize] = TrackInfo_d0error[BInfo_rftk1_index[j]];
      trk1PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
      trk1StripHit[typesize] = TrackInfo_striphit[BInfo_rftk1_index[j]];
      trk1Pt[typesize] = TrackInfo_pt[BInfo_rftk1_index[j]];
      trk1Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
      trk1Eta[typesize] = TrackInfo_eta[BInfo_rftk1_index[j]];
      trk1Phi[typesize] = TrackInfo_phi[BInfo_rftk1_index[j]];
      trk1Y[typesize] = b4P->Rapidity();

      trk2Dxy[typesize] = -1;
      trk2D0Err[typesize] = -1;
      trk2PixelHit[typesize] = -1;
      trk2StripHit[typesize] = -1;
      trk2Pt[typesize] = -1;
      trk2Chi2ndf[typesize] = -1;
      trk2Eta[typesize] = -20;
      trk2Phi[typesize] = -20;
      trk2Y[typesize] = -1;

      tktkmass[typesize] = -1;
      tktkvProb[typesize] = -1;
      tktkpt[typesize] = -1;
      tktketa[typesize] = -20;
      tktkphi[typesize] = -20;
      tktky[typesize] = -1;
      doubletmass[typesize] = -1;
      doubletpt[typesize] = -1;
      doubleteta[typesize] = -20;
      doubletphi[typesize] = -20;
      doublety[typesize] = -1;
    }  
  else if(BInfo_type[j]==5)
    {
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk2_index[j]],TrackInfo_eta[BInfo_rftk2_index[j]],TrackInfo_phi[BInfo_rftk2_index[j]],track_mass1);
      trk1Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk2_index[j]];
      trk1D0Err[typesize] = TrackInfo_d0error[BInfo_rftk2_index[j]];
      trk1PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk2_index[j]];
      trk1StripHit[typesize] = TrackInfo_striphit[BInfo_rftk2_index[j]];
      trk1Pt[typesize] = TrackInfo_pt[BInfo_rftk2_index[j]];
      trk1Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk2_index[j]]/TrackInfo_ndf[BInfo_rftk2_index[j]];
      trk1Eta[typesize] = TrackInfo_eta[BInfo_rftk2_index[j]];
      trk1Phi[typesize] = TrackInfo_phi[BInfo_rftk2_index[j]];
      trk1Y[typesize] = b4P->Rapidity();
      tk1px = b4P->Px();
      tk1py = b4P->Py();
      tk1pz = b4P->Pz();
      tk1E = b4P->E();

      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],track_mass2);
      trk2Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
      trk2D0Err[typesize] = TrackInfo_d0error[BInfo_rftk1_index[j]];
      trk2PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
      trk2StripHit[typesize] = TrackInfo_striphit[BInfo_rftk1_index[j]];
      trk2Pt[typesize] = TrackInfo_pt[BInfo_rftk1_index[j]];
      trk2Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
      trk2Eta[typesize] = TrackInfo_eta[BInfo_rftk1_index[j]];
      trk2Phi[typesize] = TrackInfo_phi[BInfo_rftk1_index[j]];
      trk2Y[typesize] = b4P->Rapidity();
      tk2px = b4P->Px();
      tk2py = b4P->Py();
      tk2pz = b4P->Pz();
      tk2E = b4P->E();

      b4P->SetPxPyPzE(tk1px+tk2px,
		      tk1py+tk2py,
		      tk1pz+tk2pz,
		      tk1E+tk2E);
      tktkmass[typesize] = b4P->Mag();
      tktketa[typesize] = b4P->Eta();
      tktkphi[typesize] = b4P->Phi();
      tktky[typesize] = b4P->Rapidity();
      tktkpt[typesize] = b4P->Pt();
      tktkvProb[typesize] = TMath::Prob(BInfo_tktk_vtxchi2[j],BInfo_tktk_vtxdof[j]);
      doubletmass[typesize] = BInfo_tktk_mass[j];
      b4P->SetXYZM(BInfo_tktk_px[j],BInfo_tktk_py[j],BInfo_tktk_pz[j],BInfo_tktk_mass[j]);
      doubletpt[typesize] = b4P->Pt();
      doubleteta[typesize] = b4P->PseudoRapidity();
      doubletphi[typesize] = b4P->Phi();
      doublety[typesize] = b4P->Rapidity();
    }
  else
    {
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],track_mass1);
      trk1Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
      trk1D0Err[typesize] = TrackInfo_d0error[BInfo_rftk1_index[j]];
      trk1PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
      trk1StripHit[typesize] = TrackInfo_striphit[BInfo_rftk1_index[j]];
      trk1Pt[typesize] = TrackInfo_pt[BInfo_rftk1_index[j]];
      trk1Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
      trk1Eta[typesize] = TrackInfo_eta[BInfo_rftk1_index[j]];
      trk1Phi[typesize] = TrackInfo_phi[BInfo_rftk1_index[j]];
      trk1Y[typesize] = b4P->Rapidity();
      tk1px = b4P->Px();
      tk1py = b4P->Py();
      tk1pz = b4P->Pz();
      tk1E = b4P->E();

      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk2_index[j]],TrackInfo_eta[BInfo_rftk2_index[j]],TrackInfo_phi[BInfo_rftk2_index[j]],track_mass2);
      trk2Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk2_index[j]];
      trk2D0Err[typesize] = TrackInfo_d0error[BInfo_rftk2_index[j]];
      trk2PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk2_index[j]];
      trk2StripHit[typesize] = TrackInfo_striphit[BInfo_rftk2_index[j]];
      trk2Pt[typesize] = TrackInfo_pt[BInfo_rftk2_index[j]];
      trk2Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk2_index[j]]/TrackInfo_ndf[BInfo_rftk2_index[j]];
      trk2Eta[typesize] = TrackInfo_eta[BInfo_rftk2_index[j]];
      trk2Phi[typesize] = TrackInfo_phi[BInfo_rftk2_index[j]];
      trk2Y[typesize] = b4P->Rapidity();
      tk2px = b4P->Px();
      tk2py = b4P->Py();
      tk2pz = b4P->Pz();
      tk2E = b4P->E();

      b4P->SetPxPyPzE(tk1px+tk2px,
		      tk1py+tk2py,
		      tk1pz+tk2pz,
		      tk1E+tk2E);
      tktkmass[typesize] = b4P->Mag();
      tktketa[typesize] = b4P->Eta();
      tktkphi[typesize] = b4P->Phi();
      tktky[typesize] = b4P->Rapidity();
      tktkpt[typesize] = b4P->Pt();
      tktkvProb[typesize] = TMath::Prob(BInfo_tktk_vtxchi2[j],BInfo_tktk_vtxdof[j]);
      doubletmass[typesize] = BInfo_tktk_mass[j];
      b4P->SetXYZM(BInfo_tktk_px[j],BInfo_tktk_py[j],BInfo_tktk_pz[j],BInfo_tktk_mass[j]);
      doubletpt[typesize] = b4P->Pt();
      doubleteta[typesize] = b4P->PseudoRapidity();
      doubletphi[typesize] = b4P->Phi();
      doublety[typesize] = b4P->Rapidity();
    }

  //gen info judgement

  if(!REAL){
    gen[typesize] = 0;//gen init
    genIndex[typesize] = -1;//gen init
    genpt[typesize] = -1;
    geneta[typesize] = -20;
    genphi[typesize] = -20;
    geny[typesize] = -1;
      
    int mGenIdxTk1=-1;
    int mGenIdxTk2=-1;
    int bGenIdxTk1=-1;
    int bGenIdxTk2=-1;
    int bGenIdxMu1=-1;
    int bGenIdxMu2=-1;
    int ujGenIdxMu1=-1;
    int ujGenIdxMu2=-1;
      
    float BId,MId,tk1Id,tk2Id;
    //tk1:positive, tk2:negtive
    if(BInfo_type[j]==1){
	  BId = 521;//B+-
	  MId = -1;
	  tk1Id = 321;//K+-
	  tk2Id = -1;
	}
    if(BInfo_type[j]==2){
	  BId = 521;//B+-
	  MId = -1;
	  tk1Id = 211;//pi+-
	  tk2Id = -1;
	}
    if(BInfo_type[j]==3){
	  BId = 511;//B0
	  MId = 310;//Ks
	  tk1Id = 211;//pi+
	  tk2Id = 211;//pi-
	}
    if(BInfo_type[j]==4){
	  BId = 511;//B0
	  MId = 313;//K*0
	  tk1Id = 321;//K+
	  tk2Id = 211;//pi-
	}
    if(BInfo_type[j]==5){
	  BId = 511;//B0
	  MId = 313;//K*0
	  tk1Id = 211;//pi+
	  tk2Id = 321;//K-
	}
    if(BInfo_type[j]==6){
	  BId = 531;//Bs
	  MId = 333;//phi
	  tk1Id = 321;//K+
	  tk2Id = 321;//K-
	}
      
    int twoTks,kStar,flagkstar=0;
    if(BInfo_type[j]==1 || BInfo_type[j]==2) twoTks=0;
    else twoTks=1;
    if(BInfo_type[j]==4 || BInfo_type[j]==5) kStar=1;
    else kStar=0;
    //int nonprompt=0
    int prompt=0;
    
    int leveltk1=0;
    int leveltk2=0;
    int levelmuon1=0;
    int levelmuon2=0;
    
    int tk1geninfo=TrackInfo_geninfo_index[BInfo_rftk1_index[j]];
    int tk1pdg=abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]);
    int mothertk1geninfo=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]];
    int mothertk1pdg=abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]);
    int grandmothertk1geninfo=GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]];
    int grandmothertk1pdg=abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]]);
    
    int tk2geninfo=TrackInfo_geninfo_index[BInfo_rftk2_index[j]];
    int tk2pdg=abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]);
    int mothertk2geninfo=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]];
    int mothertk2pdg=abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]);
    int grandmothertk2geninfo=GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]];
    int grandmothertk2pdg=abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]]);
    
    int muon1geninfo=MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
    int muon1pdg=abs(GenInfo_pdgId[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]);
    int mothermuon1geninfo=GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]];
    int mothermuon1pdg=GenInfo_pdgId[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]];
    int grandmothermuon1geninfo=GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]];
    int grandmothermuon1pdg=abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]]]);

    int muon2geninfo=MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
    int muon2pdg=abs(GenInfo_pdgId[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]);
    int mothermuon2geninfo=GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]];
    int mothermuon2pdg=GenInfo_pdgId[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]];
    int grandmothermuon2geninfo=GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]];
    int grandmothermuon2pdg=abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]]]);

    if(!twoTks) TrackDirectGenLabel(tk1geninfo,tk1pdg,mothertk1geninfo,mothertk1pdg,tk1Id,BId,mGenIdxTk1,bGenIdxTk1,leveltk1);      
    else ParticleResonanceGenLabel(tk1geninfo,tk1pdg,mothertk1geninfo,mothertk1pdg,grandmothertk1geninfo,grandmothertk1pdg,tk1Id,MId,BId,mGenIdxTk1,bGenIdxTk1,leveltk1);  
    gen[typesize]=leveltk1;

      //tk2
    if(!twoTks){  //one trk channel
	  gen[typesize]+=30;
	  mGenIdxTk2=0;
	  bGenIdxTk2=0;
	}
    else{
      ParticleResonanceGenLabel(tk2geninfo,tk2pdg,mothertk2geninfo,mothertk2pdg,grandmothertk2geninfo,grandmothertk2pdg,tk2Id,MId,BId,mGenIdxTk2,bGenIdxTk2,leveltk2);  
      gen[typesize]+=(leveltk2*10);
    }
    
    ParticleResonanceGenLabel(muon1geninfo,muon1pdg,mothermuon1geninfo,mothermuon1pdg,grandmothermuon1geninfo,grandmothermuon1pdg,13,443,BId,ujGenIdxMu1,bGenIdxMu1,levelmuon1);
    gen[typesize]+=(levelmuon1*100);
    if(levelmuon1==3) flagkstar++;
    
    ParticleResonanceGenLabel(muon2geninfo,muon2pdg,mothermuon2geninfo,mothermuon2pdg,grandmothermuon2geninfo,grandmothermuon2pdg,13,443,BId,ujGenIdxMu2,bGenIdxMu2,levelmuon2);
    gen[typesize]+=(levelmuon2*1000);
    if(levelmuon2==3) flagkstar++;

    int levelBcand=0;
    if(mGenIdxTk1!=-1 && mGenIdxTk2!=-1){
	  if(!twoTks) levelBcand=1;
	  else{
	    if(mGenIdxTk1==mGenIdxTk2) levelBcand=1;
	  }
	}
    if(bGenIdxMu1!=-1 && bGenIdxMu1==bGenIdxMu2 && bGenIdxMu1==bGenIdxTk1){
	  if(!twoTks){
	    levelBcand=2;
	    genIndex[typesize] = bGenIdxMu1;
	  }
	  else if(bGenIdxMu1==bGenIdxTk2){
	    levelBcand=2;
	    genIndex[typesize] = bGenIdxMu1;
	  }
	}
    
    gen[typesize]+=(levelBcand*10000);

    if(kStar){ //here we recover cases in which we did the wrong assumption on the K and pi mass.
    
      int leveltk1swapped=0;
      int leveltk2swapped=0;
      
      ParticleResonanceGenLabel(tk1geninfo,tk1pdg,mothertk1geninfo,mothertk1pdg,grandmothertk1geninfo,grandmothertk1pdg,tk2Id,MId,BId,mGenIdxTk1,bGenIdxTk1,leveltk1swapped);  
      ParticleResonanceGenLabel(tk2geninfo,tk2pdg,mothertk2geninfo,mothertk2pdg,grandmothertk2geninfo,grandmothertk2pdg,tk1Id,MId,BId,mGenIdxTk2,bGenIdxTk2,leveltk2swapped);  
      if(leveltk1swapped==3) flagkstar++;
      if(leveltk2swapped==3) flagkstar++;

	  if(flagkstar==4){
	    if((bGenIdxMu1!=-1)&&(bGenIdxMu1==bGenIdxMu2)&&(bGenIdxMu1==bGenIdxTk1)&&(bGenIdxMu1==bGenIdxTk2)){
		  gen[typesize]=41000; 
		}
	  }
	}

    int tgenIndex=genIndex[typesize];
    if(gen[typesize]==23333 || gen[typesize]==41000){
    
	  genpt[typesize] = GenInfo_pt[tgenIndex];
	  geneta[typesize] = GenInfo_eta[tgenIndex];
	  genphi[typesize] = GenInfo_phi[tgenIndex];
	  b4P->SetXYZM(GenInfo_pt[tgenIndex]*cos(GenInfo_phi[tgenIndex]),GenInfo_pt[tgenIndex]*sin(GenInfo_phi[tgenIndex]),GenInfo_pt[tgenIndex]*sinh(GenInfo_eta[tgenIndex]),GenInfo_mass[tgenIndex]);
	  geny[typesize] = b4P->Rapidity();
	}
  }
}

int signalGen(int Btype, int j)
{
  float BId,MId,tk1Id,tk2Id;
  int twoTks;
  //tk1:positive, tk2:negtive
  if(Btype==1)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 321;//K+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==2)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 211;//pi+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==3)
    {
      BId = 511;//B0
      MId = 310;//Ks
      tk1Id = 211;//pi+
      tk2Id = -211;//pi-
      twoTks = 1;
    }
  if(Btype==4)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = 321;//K+
      tk2Id = -211;//pi-
      twoTks = 1;
    }
  if(Btype==5)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = -321;//pi+
      tk2Id = 211;//K-
      twoTks = 1;
    }
  if(Btype==6)
    {
      BId = 531;//Bs
      MId = 333;//phi
      tk1Id = 321;//K+
      tk2Id = -321;//K-
      twoTks = 1;
    }

  int flag=0;
  if (abs(GenInfo_pdgId[j])==BId&&GenInfo_nDa[j]==2&&GenInfo_da1[j]!=-1&&GenInfo_da2[j]!=-1){
    if (abs(GenInfo_pdgId[GenInfo_da1[j]]==443)){//jpsi
	  if(GenInfo_da1[GenInfo_da1[j]]!=-1&&GenInfo_da2[GenInfo_da1[j]]!=-1){
	    if(abs(GenInfo_pdgId[GenInfo_da1[GenInfo_da1[j]]])==13&&abs(GenInfo_pdgId[GenInfo_da2[GenInfo_da1[j]]])==13){
		  if(!twoTks){if(abs(GenInfo_pdgId[GenInfo_da2[j]])==tk1Id) flag++;}
		  else{
		    if (abs(GenInfo_pdgId[GenInfo_da2[j]])==MId){
			  if(GenInfo_da1[GenInfo_da2[j]]!=-1 && GenInfo_da2[GenInfo_da2[j]]!=-1){
			    if(GenInfo_pdgId[GenInfo_da1[GenInfo_da2[j]]]==tk1Id && GenInfo_pdgId[GenInfo_da2[GenInfo_da2[j]]]==tk2Id) flag++;
			  }
			}
		  }
		}
	  }
	}
  }
  return flag;
}




void loopGMI(string infile="../Input/Bfinder_all_151_1_Y7s.root", 
             string outfile="Modified/testModified.root", bool REAL=0,bool PbpMC=0,int nEntries=0){
	  
  Bool_t IsEventSelected_pPb5TeV(bool,bool,int,double);
  Bool_t IsBplusCandidateSelected(double,double,double,double,double);
  Bool_t IsBzeroCandidateSelected(double,double,double,double,double,double,double);
    
  const char* infname;
  const char* outfname;

  if(REAL) cout<<"--- REAL DATA ---"<<endl;
  else{
    cout<<"--- MC ---"<<endl;
    if(PbpMC) cout<<"--- Pbp ---"<<endl;
    else cout<<"--- pPb ---"<<endl;
  }

  infname = infile.c_str();
  outfname = outfile.c_str();

  TFile *f = new TFile(infname);
  TTree *root = (TTree*)f->Get("demo/root");
  TTree *hlt = (TTree*)f->Get("hltanalysis/HltTree");
  if (root->GetEntries()!=hlt->GetEntries()) {
    cout <<"Inconsistent number of entries!!!"<<endl;
    cout <<"HLT tree: "<<hlt->GetEntries()<<endl;
    cout <<"Bfinder tree: "<<root->GetEntries()<<endl;
  }
    
  TFile *outf = new TFile(outfname,"recreate");

  setBranch(root);
  setHltBranch(hlt);
    
  int ifchannel[7];
  ifchannel[0] = 1; //jpsi+Kp
  ifchannel[3] = 1; //jpsi+K*(K+,pi-)
  ifchannel[4] = 1; //jpsi+K*(K-,pi+)

  
  TTree* nt0 = new TTree("ntKp","");
  buildBranch(nt0);
  TTree* nt3 = new TTree("ntKstar","");
  buildBranch(nt3);
  TTree* ntGen = new TTree("ntGen","");
  buildGenBranch(ntGen);

  cout<<"--- Tree building finished ---"<<endl;
  
  Long64_t nentries = root->GetEntries();
  Long64_t nbytes = 0;
  TVector3* bP = new TVector3;
  TVector3* bVtx = new TVector3;
  TLorentzVector* b4P = new TLorentzVector;
  TLorentzVector* b4Pout = new TLorentzVector;
  TLorentzVector bGen;
  
  int type,flag;
  int flagEvt=0;  
  int offsetHltTree=0;
  
  if (nEntries!=0) nentries=nEntries;

  for (Long64_t i=0; i<nentries;i++) {
    nbytes += root->GetEntry(i);
    flagEvt=0;
    while (flagEvt==0){
      hlt->GetEntry(i+offsetHltTree);
      if (Bfr_HLT_Event==EvtInfo_EvtNo && Bfr_HLT_Run==EvtInfo_RunNo) flagEvt=1; else offsetHltTree++;
    } 

    if (i%10000==0) cout <<i<<" / "<<nentries<<"   offset HLT:"<<offsetHltTree<<endl;

    int type1size=0;
    int type4size=0;
    
    float best,best2,temy;
    int bestindex,best2index;

    // Bplus section

    size=0;
    best=-1;
    best2=10000;
    bestindex=-1;
    best2index=-1;
    
    bool iscandselected=false;
    
    for (int j=0;j<BInfo_size;j++){
	  if(BInfo_type[j]>7) continue;
	  if(ifchannel[BInfo_type[j]-1]!=1) continue;

	  if (!MuonInfo_passMuID[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]) continue;
	  if (!MuonInfo_passMuID[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]) continue;

	  b4Pout->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);
	  temy = b4Pout->Rapidity();
	
	  if(!IsEventSelected_pPb5TeV(REAL,PbpMC,EvtInfo_RunNo,temy)) continue;
		
	  if(BInfo_mass[j]<massmincut || BInfo_mass[j]>massmaxcut) continue;
	  if(BInfo_pt[j]<ptmincut) continue;
	  iscandselected=false;
	
	  if(BInfo_type[j]==1){
	    fillTree(bP,bVtx,b4P,j,type1size,KAON_MASS,0,REAL,PbpMC);
	    iscandselected=IsBplusCandidateSelected(chi2cl[type1size],trk1Pt[type1size],mumumass[type1size],(d0[type1size]/d0Err[type1size]),cos(dtheta[type1size]));
	  	if(chi2cl[type1size]>best&&iscandselected &&HLT_PAMu3_v1){
	  	  best = chi2cl[type1size];
		  bestindex = type1size;
	    }
	    type1size++;
	  }  
    }  
    
    if(bestindex>-1){
      bestchi2 = bestindex;
      isbestchi2[bestindex] = 1;
    }
    
    nt0->Fill();

    // Bzero section

    size=0;
    best=-1;
    best2=10000;
    bestindex=-1;
    best2index=-1;

        
    for (int j=0;j<BInfo_size;j++){
	  if(BInfo_type[j]>7) continue;
	  if(ifchannel[BInfo_type[j]-1]!=1) continue;

	  if (!MuonInfo_passMuID[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]) continue;
	  if (!MuonInfo_passMuID[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]) continue;

	  b4Pout->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);
	  temy = b4Pout->Rapidity();
	
	  if(!IsEventSelected_pPb5TeV(REAL,PbpMC,EvtInfo_RunNo,temy)) continue;
		
	  if(BInfo_mass[j]<massmincut || BInfo_mass[j]>massmaxcut) continue;
	  if(BInfo_pt[j]<ptmincut) continue;
	  iscandselected=false;
	
	  if(BInfo_type[j]==4 || BInfo_type[j]==5){
	    fillTree(bP,bVtx,b4P,j,type4size,KAON_MASS,PION_MASS,REAL,PbpMC);
	    iscandselected=IsBzeroCandidateSelected(chi2cl[type4size],trk1Pt[type4size],trk2Pt[type4size],mumumass[type4size],(d0[type4size]/d0Err[type4size]),cos(dtheta[type4size]),tktkmass[type4size]);
	  	
	  	if(chi2cl[type4size]>best&&iscandselected &&HLT_PAMu3_v1){
	  	  best = chi2cl[type4size];
		  bestindex = type4size;
	    }
	    if(abs(tktkmass[type4size]-KSTAR_MASS)<best2&&iscandselected &&HLT_PAMu3_v1){
	  	  best2 = abs(tktkmass[type4size]-KSTAR_MASS);
		  best2index = type4size;
	    }
	    type4size++;
	  }  
    }  
    
    if(bestindex>-1){
      bestchi2 = bestindex;
      isbestchi2[bestindex] = 1;
    }
    
    if(best2index>-1){
	  besttktkmass = best2index;
	  isbesttktkmass[best2index] = 1;
    }
    
    nt3->Fill();

    if(!REAL){
	  Gensize = 0;
	  for (int j=0;j<GenInfo_size;j++){
	    
	    bGen.SetPtEtaPhiM(GenInfo_pt[j],GenInfo_eta[j],GenInfo_phi[j],GenInfo_mass[j]);
	    flag=0;
	    
	    for(type=1;type<8;type++){
		  if(signalGen(type,j)) {
            flag=type;
		    break;
          }
	    }
	    
	    Genmu1pt[j] = -1; Genmu1eta[j] = -20; Genmu1phi[j] = -20; Genmu1p[j] = -1;
	    Genmu2pt[j] = -1; Genmu2eta[j] = -20; Genmu2phi[j] = -20; Genmu2p[j] = -1;
	    Gentk1pt[j] = -1; Gentk1eta[j] = -20; Gentk1phi[j] = -20; 
	    Gentk2pt[j] = -1; Gentk2eta[j] = -20; Gentk2phi[j] = -20;

        if(flag!=0){
          Genmu1pt[j] = GenInfo_pt[GenInfo_da1[GenInfo_da1[j]]];
          Genmu1eta[j] = GenInfo_eta[GenInfo_da1[GenInfo_da1[j]]];
          Genmu1phi[j] = GenInfo_phi[GenInfo_da1[GenInfo_da1[j]]];
          Genmu1p[j] = Genmu1pt[j]*cosh(Genmu1eta[j]);
          Genmu2pt[j] = GenInfo_pt[GenInfo_da2[GenInfo_da1[j]]];
          Genmu2eta[j] = GenInfo_eta[GenInfo_da2[GenInfo_da1[j]]];
		  Genmu2phi[j] = GenInfo_phi[GenInfo_da2[GenInfo_da1[j]]];
          Genmu2p[j] = Genmu2pt[j]*cosh(Genmu2eta[j]);
		
		  if(flag==1||flag==2){
		  
		    Gentk1pt[j] = GenInfo_pt[GenInfo_da2[j]];
		    Gentk1eta[j] = GenInfo_eta[GenInfo_da2[j]];
		    Gentk1phi[j] = GenInfo_phi[GenInfo_da2[j]];
		    }else{
		    Gentk1pt[j] = GenInfo_pt[GenInfo_da1[GenInfo_da2[j]]];
		    Gentk1eta[j] = GenInfo_eta[GenInfo_da1[GenInfo_da2[j]]];
		    Gentk1phi[j] = GenInfo_phi[GenInfo_da1[GenInfo_da2[j]]];
		    Gentk2pt[j] = GenInfo_pt[GenInfo_da2[GenInfo_da2[j]]];
		    Gentk2eta[j] = GenInfo_eta[GenInfo_da2[GenInfo_da2[j]]];
		    Gentk2phi[j] = GenInfo_phi[GenInfo_da2[GenInfo_da2[j]]];
		  }
        }
	    
	    Gensize = GenInfo_size;
	    Geny[j] = bGen.Rapidity();
	    Geneta[j] = bGen.Eta();
	    Genphi[j] = bGen.Phi();
	    Genpt[j] = bGen.Pt();
	    GenpdgId[j] = GenInfo_pdgId[j];
	    GenisSignal[j] = flag;
	  }
  	  ntGen->Fill(); 
	}
  }
  outf->Write();
  outf->Close();
}

Bool_t IsEventSelected_pPb5TeV(bool myisReal, bool myispbpmc, int myevtnumber, double myrapidity){

  bool flag=false;
  
  if(myisReal){
    if((((myevtnumber>=210498&&myevtnumber<=211256&&abs(myrapidity+0.465)<1.93)||(myevtnumber>=211313&&myevtnumber<=211631&&abs(myrapidity-0.465)<1.93)))) flag=true;
  }
  if(!myisReal){
    if((myispbpmc==0)&&(!(abs(myrapidity+0.465)>=1.93))) flag=true;
    if((myispbpmc==1)&&(!(abs(myrapidity-0.465)>=1.93))) flag=true;
  }
  return flag;
}


Bool_t IsBplusCandidateSelected(double mychi2cl,double mytrk1Pt,double mymumumass,double myd0d0err,double mycostheta){
  
  bool flag=false;
  if(mychi2cl>chi2clCutBplus&&mytrk1Pt>trkPtCutBplus&&myd0d0err>d0d0ErrCutBplus&&abs(mymumumass-3.096916)<mumumassCut&&mycostheta>cosdthetaCutBplus) flag=true;
  return flag;
}

Bool_t IsBzeroCandidateSelected(double mychi2cl,double mytrk1Pt,double mytrk2Pt,double mymumumass,double myd0d0err,double mycostheta,double invmasstktk){
    
  bool flag=false;
  if(mychi2cl>chi2clCutBzero&&mytrk1Pt>trkPtCutBzero&&mytrk2Pt>trkPtCutBzero&&myd0d0err>d0d0ErrCutBzero&&abs(mymumumass-3.096916)<mumumassCut&&mycostheta>cosdthetaCutBzero&&abs(invmasstktk-0.89594)<invmasstktkBzero) flag=true;
  return flag;
}



void TrackDirectGenLabel(int myparticlegeninfo,int mypdgparticle,int myBmesongeninfo,int mypdgBmeson,
                         int pdgparticle, int pdgBmeson,int &mymindex,int &myBindex, int &mylabel){
                                                
  if(myparticlegeninfo>-1){  
    mylabel=0;
	if(mypdgparticle==pdgparticle){
	  mylabel=1;
	  if(myBmesongeninfo>-1){ 
	    mymindex=0.;
		if(mypdgBmeson==pdgBmeson){
		  mylabel=3;
		  myBindex=myBmesongeninfo;
		}
	  }
	}
  }
}

void ParticleResonanceGenLabel(int myparticlegeninfo,int mypdgparticle,
                            int mymothermesongeninfo,int mypdgmother,
                            int myBmesongeninfo,int mypdgBmeson,
                            int pdgparticle, int pdgmother, int pdgBmeson,
                            int &mymindex,int &myBindex, int &mylabel){
                                                
  if(myparticlegeninfo>-1){  
    mylabel=0;
	if(mypdgparticle==pdgparticle){
	  mylabel=1;
	  if(mymothermesongeninfo>-1){ 
		if(mypdgmother==pdgmother){
		  mylabel=2;
		  if(myBmesongeninfo>-1){
		    if(mypdgBmeson==pdgBmeson){
		    mylabel=3;
		    myBindex=myBmesongeninfo;
		    }
		  }
		  mymindex=mymothermesongeninfo;
		}
	  }
	}
  }
}
