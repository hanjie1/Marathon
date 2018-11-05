#include "SetCut.h"
#include "GetTrees.h"
void plot_kin11(){
     TChain *T1;
     T1=GetTree(1545,11,"T");

     TChain *T2;
     T2=GetTree(1546,11,"T");

     TCanvas *c1=new TCanvas("c1");
     c1->Divide(2,2);
     c1->cd(1);
     TH1F *hxbj=new TH1F("hxbj","xbj distribution",100,0,1);
     TH1F *hxbj1=new TH1F("hxbj1","xbj distribution",100,0,1);
     T1->Draw("EKLx.x_bj>>hxbj",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("EKLx.x_bj>>hxbj1",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hxbj->SetLineColor(2);
    
//     TCanvas *c2=new TCanvas("c2");
//     c2->Divide(2,1);
     c1->cd(2);
     TH1F *hNu=new TH1F("hNu","Nu distribution",500,7,8);
     TH1F *hNu1=new TH1F("hNu1","Nu distribution",500,7,8);
     T1->Draw("EKLx.nu>>hNu",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("EKLx.nu>>hNu1",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hNu->SetLineColor(2);
/*
     c2->cd(2);
     TH1F *hangle=new TH1F("hangle","angle distribution",500,27,32);
     TH1F *hangle1=new TH1F("hangle1","angle distribution",500,27,32);
     T1->Draw("EKLx.angle*180/3.14>>hangle",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("EKLx.angle*180/3.14>>hangle1",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hangle->SetLineColor(2);

     TCanvas *c3=new TCanvas("c3");
     c3->Divide(2,1);
     c3->cd(1);
     TH1F *hQsq=new TH1F("hQsq","Qsq distribution",500,7,10);
     TH1F *hQsq1=new TH1F("hQsq1","Qsq distribution",500,7,10);
     T1->Draw("EKLx.Q2>>hQsq",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("EKLx.Q2>>hQsq1",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hQsq->SetLineColor(2);
*/
     c1->cd(3);
     TH1F *hEp=new TH1F("hEp","L.gold.p distribution",200,2.5,4);
     TH1F *hEp1=new TH1F("hEp1","L.gold.p distribution",200,2.5,4);
     T1->Draw("L.gold.p>>hEp",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("L.gold.p>>hEp1",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hEp->SetLineColor(2);

     TCanvas *c4=new TCanvas("c4");
     c4->Divide(2,2);
     c4->cd(1);
     TH1F *hQ1=new TH1F("hQ1","Q1 distribution",1000,-1.15,-1);
     TH1F *hQ11=new TH1F("hQ11","Q1 distribution",1000,-1.15,-1);
     T1->Draw("HacL_Q1_LS450_FLD_DATA>>hQ1",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("HacL_Q1_LS450_FLD_DATA>>hQ11",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hQ1->SetLineColor(2);

     c4->cd(2);
     TH1F *hQ2=new TH1F("hQ2","Q2 distribution",1000,0.6,0.8);
     TH1F *hQ21=new TH1F("hQ21","Q2 distribution",1000,0.6,0.8);
     T1->Draw("HacL_Q2_LS450_FLD_DATA>>hQ2",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("HacL_Q2_LS450_FLD_DATA>>hQ21",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hQ2->SetLineColor(2);

     c4->cd(3);
     TH1F *hQ3=new TH1F("hQ3","Q3 distribution",1000,-0.7,-0.5);
     TH1F *hQ31=new TH1F("hQ31","Q3 distribution",1000,-0.7,-0.5);
     T1->Draw("HacL_Q3_LS450_FLD_DATA>>hQ3",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("HacL_Q3_LS450_FLD_DATA>>hQ31",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hQ3->SetLineColor(2);

     c4->cd(4);
     TH1F *hE0=new TH1F("hE0","E0 distribution",500,10580,10600);
     TH1F *hE01=new TH1F("hE01","E0 distribution",500,10580,10600);
     T1->Draw("HALLA_p>>hE0",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("HALLA_p>>hE01",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hE0->SetLineColor(2);

     TCanvas *c5=new TCanvas("c5");
     c5->Divide(2,2);
     c5->cd(1);
     TH1F *htrx=new TH1F("htrx","trx distribution",500,-1,1);
     TH1F *htrx1=new TH1F("htrx1","trx distribution",500,-1,1);
     T1->Draw("L.tr.x>>htrx",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("L.tr.x>>htrx1",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     htrx->SetLineColor(2);

     c5->cd(2);
     TH1F *htry=new TH1F("htry","try distribution",500,-0.5,0.5);
     TH1F *htry1=new TH1F("htry1","try distribution",500,-0.5,0.5);
     T1->Draw("L.tr.y>>htry",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("L.tr.y>>htry1",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     htry->SetLineColor(2);

     c5->cd(3);
     TH1F *hth=new TH1F("hth","th distribution",500,-0.5,0.5);
     TH1F *hth1=new TH1F("hth1","th distribution",500,-0.5,0.5);
     T1->Draw("L.tr.th>>hth",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("L.tr.th>>hth1",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hth->SetLineColor(2);
     
     c5->cd(4);
     TH1F *hph=new TH1F("hph","ph distribution",500,-0.5,0.5);
     TH1F *hph1=new TH1F("hph1","ph distribution",500,-0.5,0.5);
     T1->Draw("L.tr.ph>>hph",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     T2->Draw("L.tr.ph>>hph1",trigger2+CK+Ep+beta+ACC+VZ+TRK,"same");
     hph->SetLineColor(2);
}
