void compare(){

TFile *infOriginal = new TFile("../Original/testOriginal.root");
TFile *infModified = new TFile("../Modified/testModified.root");
TTree *ntOriginal = (TTree*) infOriginal->Get("ntKp");
TTree *ntModified = (TTree*) infModified->Get("ntKp");

TTree *ntOriginalGen = (TTree*) infOriginal->Get("ntGen");
TTree *ntModifiedGen = (TTree*) infModified->Get("ntGen");


TH1D *hchi2bestOriginal = new TH1D("hchi2bestOriginal","",100,0.,20.);
TH1D *hchi2bestModified = new TH1D("hchi2bestModified","",100,0.,20.);

TH1D *hmassOriginal = new TH1D("hmassOriginal","",100,5.,6.);
TH1D *hmassModified = new TH1D("hmassModified","",100,5.,6.);

TH1D *hsizeGenOriginal = new TH1D("hsizeGenOriginal","",100,0.,500.);
TH1D *hsizeGenModified = new TH1D("hsizeGenModified","",100,0.,500.);


ntOriginal->Project("hchi2bestOriginal","isbestchi2","");   
ntModified->Project("hchi2bestModified","isbestchi2","");   
ntOriginal->Project("hmassOriginal","mass","");   
ntModified->Project("hmassModified","mass","");   

ntOriginalGen->Project("hsizeGenOriginal","size","");   
ntModifiedGen->Project("hsizeGenModified","size","");   


TCanvas*canvas=new TCanvas("canvas","canvas",800,500);
canvas->Divide(3,2);
canvas->cd(1);
hchi2bestOriginal->Draw();
canvas->cd(2);
hchi2bestModified->Draw();
canvas->cd(3);
hmassOriginal->Draw();
canvas->cd(4);
hmassModified->Draw();
canvas->cd(5);
hsizeGenOriginal->Draw();
canvas->cd(6);
hsizeGenModified->Draw();


canvas->SaveAs("canvascompare.pdf");
canvas->SaveAs("canvascompare.root");

}