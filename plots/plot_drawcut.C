#include "SetCut.h"
#include "GetTrees.h"


void plot()
{
     TString filename;
     cout<<"Input file name: "<<endl;
     cin>>filename;

     TChain *T=new TChain("T");
     T=GetFiles(filename,"T");
  
     Double_t nentries=T->GetEntries(CK+Ep+TRK+trigger2+beta+FP_ACC);
     cout<<nentries<<endl;

     gStyle->SetOptStat(1111110);

     TCanvas *c1=new TCanvas("c1","c1");

     TH2F *x_th =new TH2F("x_th","L.tr.x:L.tr.th(CK+Ep+TRK+trigger2+beta+FP_ACC)",200,-0.2,0.2,200,-1.0,1.0);
     T->Draw("L.tr.x:L.tr.th>>x_th",CK+Ep+TRK+trigger2+beta+FP_ACC,"COLZ");
     TCutG* cutg;
     cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));
     c1->Update();
     cutg->SetName("cut_FP");
     cutg->SetVarX("L.tr.th");
     cutg->SetVarY("L.tr.x");
     cutg->SetLineColor(kMagenta);
     cutg->SetLineWidth(5);
     cutg->Draw("PL");
     c1->Update();

     TCut FP_ACC_draw=FP_ACC+TCut("cut_FP");
     (TCutG*)gROOT->FindObject("cut_FP");

     TCanvas *c2=new TCanvas("c2","c2");
     c2->Divide(2,2);
     c2->cd(1);
     TH1F *x_fp=new TH1F("x_fp","x_fp(CK+Ep+TRK+trigger2+beta+FP_ACC_draw)",200,-0.7,0.7);
     T->Draw("L.tr.x>>x_fp",CK+Ep+TRK+trigger2+beta+FP_ACC_draw);

     c2->cd(2);
     TH1F *y_fp=new TH1F("y_fp","y_fp(CK+Ep+TRK+trigger2+beta+FP_ACC_draw)",200,-0.06,0.05);
     T->Draw("L.tr.y>>y_fp",CK+Ep+TRK+trigger2+beta+FP_ACC_draw);
     
     c2->cd(3);
     TH1F *th_fp=new TH1F("th_fp","th_fp(CK+Ep+TRK+trigger2+beta+FP_ACC_draw)",200,-0.14,0.12);
     T->Draw("L.tr.th>>th_fp",CK+Ep+TRK+trigger2+beta+FP_ACC_draw);

     c2->cd(4);
     TH1F *ph_fp=new TH1F("ph_fp","ph_fp(CK+Ep+TRK+trigger2+beta+FP_ACC_draw)",200,-0.05,0.055);
     T->Draw("L.tr.ph>>ph_fp",CK+Ep+TRK+trigger2+beta+FP_ACC_draw);

     TCanvas *c3=new TCanvas("c3","c3");
     TH2F *y_ph =new TH2F("y_ph","L.tr.y:L.tr.ph(CK+Ep+TRK+trigger2+beta+FP_ACC_draw)",200,-0.07,0.07,200,-0.09,0.09);
     T->Draw("L.tr.y:L.tr.ph>>y_ph",CK+Ep+TRK+trigger2+beta+FP_ACC_draw,"COLZ");

     TCanvas *c4=new TCanvas("c4","c4");
     c4->Divide(2,2);
     c4->cd(1); 
     TH1F *targ_th=new TH1F("targ_th","targ_th(CK+Ep+TRK+trigger2+beta+FP_ACC_draw)",200,-0.1,0.1);
     T->Draw("L.tr.tg_th>>targ_th",CK+Ep+TRK+trigger2+beta+FP_ACC_draw);

     c4->cd(2);
     TH1F *targ_ph=new TH1F("targ_ph","targ_ph(CK+Ep+TRK+trigger2+beta+FP_ACC_draw)",200,-0.05,0.05);
     T->Draw("L.tr.tg_ph>>targ_ph",CK+Ep+TRK+trigger2+beta+FP_ACC_draw);

     c4->cd(3);
     TH1F *targ_dp=new TH1F("targ_dp","dp/p(CK+Ep+TRK+trigger2+beta+FP_ACC_draw)",200,-0.06,0.06);
     T->Draw("L.tr.tg_dp>>targ_dp",CK+Ep+TRK+trigger2+beta+FP_ACC_draw);

     c4->cd(4);
     TH2F *th_ph=new TH2F("th_ph","targ th:ph (CK+Ep+TRK+trigger2+beta+FP_ACC_draw)",200,-0.05,0.06,200,-0.1,0.1);
     T->Draw("L.tr.tg_th:L.tr.tg_ph>>th_ph",CK+Ep+TRK+trigger2+beta+FP_ACC_draw,"COLZ");

     c1->Print(Form("%s.pdf[",filename.Data()));
     c1->Print(Form("%s.pdf",filename.Data()));
     c2->Print(Form("%s.pdf",filename.Data()));
     c3->Print(Form("%s.pdf",filename.Data()));
     c4->Print(Form("%s.pdf",filename.Data()));
     c4->Print(Form("%s.pdf]",filename.Data()));
     

}
