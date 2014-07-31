#include "utilities.h"

double luminosity=34.8*1e-3;
double setparam0=100.;
double setparam1=5.28;
double setparam2=0.03;
double fixparam1=5.279;
double setparam3=0.03;
double fixparam2=0.04;

TString inputdata="Original/testOriginal.root";

//TString cut_kpi="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&mass>5&&mass<6&& isbestchi2&&trk1Pt>0.7&&trk2Pt>0.7&&chi2cl>1.65e-01&&(d0/d0Err)>4.16&&cos(dtheta)>7.50e-01&&abs(tktkmass-0.89594)<2.33e-01"; 
//TString cut_kpi="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&mass>5&&mass<6&&trk1Pt>0.7&&trk2Pt>0.7&&chi2cl>1.65e-01&&(d0/d0Err)>4.16&&cos(dtheta)>7.50e-01&&abs(tktkmass-0.89594)<2.33e-01"; 

TString cut_kpi="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&mass>5&&mass<6&& isbestchi2&&trk1Pt>0.7&&trk2Pt>0.7&&chi2cl>1.65e-01&&(d0/d0Err)>4.16&&cos(dtheta)>7.50e-01&&abs(tktkmass-0.89594)<2.33e-01"; 

TString seldata_2y_kpi=Form("((Run>1&&Run<12&&abs(y-0.465)<1.93)||(Run<=1&&abs(y+0.465)<1.93)||(Run>=210498&&Run<=211256&&abs(y+0.465)<1.93)||(Run>=211313&&Run<=211631&&abs(y-0.465)<1.93))&&%s",cut_kpi.Data());


void fitB0(TString infname="",bool doweight = 1){
  if (doweight==0) weight="1";
  if (infname=="") infname=inputdata.Data();
  TFile *inf = new TFile(infname.Data());
  TTree *nt = (TTree*) inf->Get("ntKstar");
  
   TH1D *h = new TH1D("h","",30,5.03,5.93);
   nt->Project("h","mass",Form("%s&&pt>%f&&pt<%f",seldata_2y_kpi.Data(),10.,60.));   
   
  TCanvas*canvas=new TCanvas("canvas","canvas",1000,500);
  canvas->Divide(3,1);
  canvas->cd(1);
  h->Draw();
  canvas->SaveAs("canvasfitB0.pdf");

}
