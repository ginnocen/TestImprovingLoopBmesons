#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "utilities.h"
#include "TH1D.h"
#define NUM_BX 9000

Float_t mumumassCut=0.15;
Float_t trkPtCut=0.7;
Float_t chi2clCut=1.65e-01;
Float_t d0d0ErrCut=4.16;
Float_t cosdthetaCut=7.50e-01;
Float_t tktkmassCut=2.33e-01;

//if(chi2cl[type4size]>best&&(HLT_PAMu3_v1)&&abs(mumumass[type4size]-3.096916)<0.15&&trk1Pt[type4size]>0.7&&trk2Pt[type4size]>0.7&&chi2cl[type4size]>1.65e-01&&(d0[type4size]/d0Err[type4size])>4.16&&cos(dtheta[type4size])>7.50e-01&&abs(tktkmass[type4size]-0.89594)<2.33e-01)


/*
Float_t trkPtCut=0.7;
Float_t chi2clCut=9.94e-02;
Float_t d0d0ErrCut=6.08;
Float_t cosdthetaCut=7.93e-01;
Float_t mumumassCut=0.15;
Float_t tktkmassCut=0.10;
Float_t mintktkmassCut=0.85;
*/

const int nBins = 3;
double ptBins[nBins+1] = {10,15,20,60};


void LoopB0(TString cutconfig="newcutopt8"){

  TH1D* Loop(TTree*,double,double,int);
  TH1D* LoopGen(TTree*,double,double);
  
  TString infname_data="Original/testOriginal.root";

  TFile *inf_data = new TFile(infname_data.Data());
  
  TTree *nt_data = (TTree*) inf_data->Get("ntKstar");
    
  TH1D *hMass1=Loop(nt_data,10.,60.,1);
    
  TCanvas*canvas=new TCanvas("canvas","canvas",1000,500);
  canvas->Divide(3,1);
  canvas->cd(1);
  hMass1->Draw();
  canvas->SaveAs("canvasLoop.pdf");
  
  TFile*fB0output=new TFile(Form("ResultsLoopBzero/B0output_%s.root",cutconfig.Data()),"recreate");
  fB0output->cd();
  hMass1->SetName("hMass1");
  hMass1->Write();
  fB0output->Close();

}

TH1D* Loop(TTree* ntuple,double ptmin,double ptmax,int option=1){

  //double funcweight(double);

  bool cut_yvsRun,cut_pt,cut_mumumass,cut_HLT_PAMu3_v1,cut_mass;
  bool cut_trk1Pt,cut_trk2Pt,cut_chi2cl,cut_d0d0err;
  bool cut_dtheta,cut_tktkmass;

  int Run,HLT_PAMu3_v1,size,Event;
  Float_t pt[NUM_BX];
  Float_t mumumass[NUM_BX];
  Float_t y[NUM_BX];
  Float_t mass[NUM_BX];
  Float_t trk1Pt[NUM_BX];
  Float_t trk2Pt[NUM_BX];
  Float_t chi2cl[NUM_BX];
  Float_t d0[NUM_BX];
  Float_t d0Err[NUM_BX];
  Float_t dtheta[NUM_BX];
  Float_t tktkmass[NUM_BX];
  Int_t isbestchi2[NUM_BX];
  Float_t gen[NUM_BX];

  int bestchi2index;
  Float_t bestchi2;

  ntuple->SetBranchAddress("Run",&Run);
  ntuple->SetBranchAddress("y",y);
  ntuple->SetBranchAddress("HLT_PAMu3_v1",&HLT_PAMu3_v1);
  ntuple->SetBranchAddress("size",&size);
  ntuple->SetBranchAddress("mumumass",mumumass);
  ntuple->SetBranchAddress("mass",mass);
  ntuple->SetBranchAddress("isbestchi2",isbestchi2);
  ntuple->SetBranchAddress("trk1Pt",trk1Pt);
  ntuple->SetBranchAddress("trk2Pt",trk2Pt);
  ntuple->SetBranchAddress("chi2cl",chi2cl);
  ntuple->SetBranchAddress("d0Err",d0Err);
  ntuple->SetBranchAddress("d0",d0);
  ntuple->SetBranchAddress("dtheta",dtheta);
  ntuple->SetBranchAddress("tktkmass",tktkmass);
  ntuple->SetBranchAddress("Event",&Event);
  ntuple->SetBranchAddress("pt",pt);
  if(option==2 || option==3 || option==4) ntuple->SetBranchAddress("gen",gen);


  TH1D *h = new TH1D("h","h",30,5.03,5.93);
  TH1D *hPtMC = new TH1D("hPtMC","",nBins,ptBins);

  Int_t entries = (Int_t)ntuple->GetEntries();
  
  for (int i=0; i<entries; i++) {
    ntuple->GetEntry(i);
        
    bestchi2index=-999;
    bestchi2=-999.;
    
    for(int j=0;j<size;j++){
    
      cut_pt=((pt[j]>ptmin)&&(pt[j]<ptmax));
      cut_d0d0err=((d0[j]/d0Err[j])>d0d0ErrCut);
      if(!cut_pt) continue;
      if(!cut_d0d0err) continue;
      
      cut_yvsRun=false;
	  if(option==1 || option==2) cut_yvsRun=((Run>1&&Run<12&&abs(y[j]-0.465)<1.93)||(Run<=1&&abs(y[j]+0.465)<1.93)||(Run>=210498&&Run<=211256&&abs(y[j]+0.465)<1.93)||(Run>=211313&&Run<=211631&&abs(y[j]-0.465)<1.93));
   	  if(option==3) cut_yvsRun=((Run>1&&Run<12&&abs(y[j]-0.465)<1.93)||(Run<=1&&abs(y[j]+0.465)<1.93)||(Run>=210498&&Run<=211256&&abs(y[j]+0.465)<1.93)||(Run>=211313&&Run<=211631&&abs(y[j]-0.465)<1.93))&&(gen[j]==23333||gen[j]==41000);
   	  if(option==4) cut_yvsRun=((Run>1&&Run<12&&abs(y[j]-0.465)<1.93)||(Run<=1&&abs(y[j]+0.465)<1.93)||(Run>=210498&&Run<=211256&&abs(y[j]+0.465)<1.93)||(Run>=211313&&Run<=211631&&abs(y[j]-0.465)<1.93))&&(!(gen[j]==23333||gen[j]==41000));

	  cut_HLT_PAMu3_v1=(HLT_PAMu3_v1);
	  cut_mass=((mass[j]>5.)&&(mass[j]<6.));
	  cut_trk1Pt=(trk1Pt[j]>trkPtCut);
	  cut_trk2Pt=(trk2Pt[j]>trkPtCut);
	  cut_chi2cl=(chi2cl[j]>chi2clCut);
	  
  	  cut_dtheta=(cos(dtheta[j])>cosdthetaCut);
  	  cut_mumumass=abs(mumumass[j]-3.096916)<mumumassCut;
  	  cut_tktkmass=TMath::Abs(tktkmass[j]-0.89594)<tktkmassCut;
  	    	  
      //if(cut_yvsRun&&cut_HLT_PAMu3_v1&&cut_mass&&cut_trk1Pt&&cut_trk2Pt&&cut_chi2cl&&cut_dtheta&&cut_mumumass&&cut_tktkmass){
       if(cut_HLT_PAMu3_v1&&cut_mass&&cut_mumumass&&cut_trk1Pt&&cut_trk2Pt&&cut_chi2cl&&cut_dtheta&&cut_mumumass&&cut_tktkmass){
     
        if(chi2cl[j]>bestchi2) {bestchi2=chi2cl[j]; bestchi2index=j;}
       // h->Fill(mass[j]);  	  
  	  }//candidate seleection
  	}//loop over candidates
  	
  	if(bestchi2index>-1){
  	//if(bestchi2index>-1){
  	  h->Fill(mass[bestchi2index]);
    }
  }// loop over events
  if(option==1 || option==2 || option==4) return h;
  else return hPtMC; 
}


TH1D* LoopGen(TTree* ntuple,double ptmin,double ptmax){

  int Run,size;
  bool cut_yvsRun,cut_pt;
  Float_t gen[NUM_BX];
  Float_t pdgId[NUM_BX];
  Float_t isSignal[NUM_BX];
  Float_t y[NUM_BX];
  Float_t pt[NUM_BX];
  

  ntuple->SetBranchAddress("size",&size);
  ntuple->SetBranchAddress("Run",&Run);
  ntuple->SetBranchAddress("gen",gen);
  ntuple->SetBranchAddress("pdgId",pdgId);
  ntuple->SetBranchAddress("isSignal",isSignal);
  ntuple->SetBranchAddress("y",y);
  ntuple->SetBranchAddress("pt",pt);

  TH1D *hPtMCGen = new TH1D("hPtMCGen","",nBins,ptBins);

  Int_t entries = (Int_t)ntuple->GetEntries();
  
  for (int i=0; i<entries; i++) {
    ntuple->GetEntry(i);
        
    for(int j=0;j<size;j++){
    
      cut_yvsRun=((Run<=1&&abs(y[j]+0.465)<1.93)||(Run>1&&abs(y[j]-0.465)<1.93))&&abs(pdgId[j])==511&&isSignal[j]!=0;
      cut_pt=((pt[j]>ptmin)&&(pt[j]<ptmax));
      if(cut_yvsRun&&cut_pt) hPtMCGen->Fill(pt[j]);
          	  
  	}//loop over candidates
  }// loop over events

  return hPtMCGen; 
}

