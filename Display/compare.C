void compare(){

TFile *infOriginal = new TFile("../Original/testOriginal.root");
TFile *infModified = new TFile("../Modified/testModified.root");
TTree *ntOriginal = (TTree*) infOriginal->Get("ntKp");
TTree *ntModified = (TTree*) infModified->Get("ntKp");


TH1D *hchi2bestOriginal = new TH1D("hchi2bestOriginal","",100,0.,20.);
TH1D *hchi2bestModified = new TH1D("hchi2bestModified","",100,0.,20.);

TH1D *hmassOriginal = new TH1D("hmassOriginal","",100,5.,6.);
TH1D *hmassModified = new TH1D("hmassModified","",100,5.,6.);



ntOriginal->Project("hchi2bestOriginal","isbestchi2","");   
ntModified->Project("hchi2bestModified","isbestchi2","");   


ntOriginal->Project("hmassOriginal","mass","");   
ntModified->Project("hmassModified","mass","");   


TCanvas*canvas=new TCanvas("canvas","canvas",800,500);
canvas->Divide(2,2);
canvas->cd(1);
hchi2bestOriginal->Draw();
canvas->cd(2);
hchi2bestModified->Draw();
canvas->cd(3);
hmassOriginal->Draw();
canvas->cd(4);
hmassModified->Draw();

canvas->SaveAs("canvascompare.pdf");
canvas->SaveAs("canvascompare.root");

}