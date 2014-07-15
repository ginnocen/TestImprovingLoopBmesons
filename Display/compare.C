void compare(){

TFile *infOriginal = new TFile("../Original/testOriginal.root");
TFile *infModified = new TFile("../Modified/testModified.root");
TTree *ntOriginal = (TTree*) infOriginal->Get("ntKstar");
TTree *ntModified = (TTree*) infModified->Get("ntKstar");


TH1D *hchi2bestOriginal = new TH1D("hchi2bestOriginal","",100,0.,2.);
TH1D *hchi2bestModified = new TH1D("hchi2bestModified","",100,0.,2.);


ntOriginal->Project("hchi2bestOriginal","bestchi2","");   
ntModified->Project("hchi2bestModified","bestchi2","");   
TCanvas*canvas=new TCanvas("canvas","canvas",800,500);
canvas->Divide(2,1);
canvas->cd(1);
hchi2bestOriginal->Draw();
canvas->cd(2);
hchi2bestModified->SetLineColor(2);
hchi2bestModified->Draw();
canvas->SaveAs("canvascompare.pdf");
canvas->SaveAs("canvascompare.root");

}