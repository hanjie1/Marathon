#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp()
{
   Double_t x[MAXBIN]={0.0},Ratio[MAXBIN]={0.0},Rerr[MAXBIN]={0.0};
   Double_t x1[MAXBIN]={0.0},Ratio1[MAXBIN]={0.0},Rerr1[MAXBIN]={0.0};
   Double_t x2[MAXBIN]={0.0},Ratio2[MAXBIN]={0.0},Rerr2[MAXBIN]={0.0};
   Double_t x3[MAXBIN]={0.0},Ratio3[MAXBIN]={0.0},Rerr3[MAXBIN]={0.0};
   Double_t x4[MAXBIN]={0.0},Ratio4[MAXBIN]={0.0},Rloerr4[MAXBIN]={0.0},Rhierr4[MAXBIN]={0.0};
   Double_t x_KP[1100]={0.0},F2p_KP[1100]={0.0},F2d_KP[1100]={0.0};

   TString Rfile="newbin/Dp_final.dat";
   int nbin=ReadYield(Rfile,x,Ratio,Rerr); 
   Rfile="Model/F2dp_Whitlow.out";
   int nbin1=ReadModel(Rfile,x1,Ratio1,Rerr1); 
   Rfile="Model/F2dp_Bodek.out";
   int nbin2=ReadModel(Rfile,x2,Ratio2,Rerr2); 
   Rfile="Model/CJ_dp.csv";
   int nbin3=ReadCJ(Rfile,x3,Ratio3,Rerr3); 
   Rfile="Model/F2dp_NMC.out";
   int nbin4=ReadNMC(Rfile,x4,Ratio4,Rloerr4,Rhierr4); 
   Rfile="Model/F2dis_os1tm1ht1mec1_Dav18_He3Salme";
   int nbin5=ReadKP(Rfile,x_KP,F2p_KP,F2d_KP); 

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hr_norm=new TGraphErrors(); //MARATHON normalization error
   TGraphErrors *hratio1=new TGraphErrors(); //Whtilow F2d/F2p without normalization error
   TGraphErrors *hr1_norm=new TGraphErrors(); //Whtilow normalization error
   TGraphErrors *hratio2=new TGraphErrors(); //Bodek F2d/F2p
   TGraphErrors *hratio3=new TGraphErrors(); //CJ F2d/F2p
   TGraphAsymmErrors *hratio4=new TGraphAsymmErrors(); //NMC F2d/F2p
   TGraph *hratio5=new TGraph();  //KP F2d/F2p
   
   for(int ii=0;ii<nbin;ii++){
	hratio->SetPoint(ii,x[ii],Ratio[ii]);
	hratio->SetPointError(ii,0,Rerr[ii]);
	hr_norm->SetPoint(ii,x[ii],1.4);
	hr_norm->SetPointError(ii,0,0.79/100*Ratio[ii]);
   } 

   for(int ii=0;ii<nbin1;ii++){
	hratio1->SetPoint(ii,x[ii],Ratio1[ii]);
	hratio1->SetPointError(ii,0,Rerr1[ii]);
	hr1_norm->SetPoint(ii,x[ii],1.4);
	hr1_norm->SetPointError(ii,0,Ratio1[ii]*1.0/100.0);
   } 

   for(int ii=0;ii<nbin2;ii++){
	hratio2->SetPoint(ii,x[ii],Ratio2[ii]);
	hratio2->SetPointError(ii,0,Rerr2[ii]);
   } 
   for(int ii=0;ii<nbin3;ii++){
	hratio3->SetPoint(ii,x[ii],2.0*Ratio3[ii]);
	hratio3->SetPointError(ii,0,2.0*Rerr3[ii]);
   } 
   for(int ii=0;ii<nbin4;ii++){
	hratio4->SetPoint(ii,x[ii],Ratio4[ii]);
//	hratio4->SetPointEYlow(ii,Rloerr4[ii]);
//	hratio4->SetPointEYhigh(ii,Rhierr4[ii]);
	hratio4->SetPointEYlow(ii,0.0);
	hratio4->SetPointEYhigh(ii,0.0);
   } 
   int nn=0;
   for(int ii=0;ii<nbin5;ii++){
	if(x_KP[ii]<x[0] || x_KP[ii]>x[7])continue;
	hratio5->SetPoint(nn,x_KP[ii],2.0*F2d_KP[ii]/F2p_KP[ii]);
	nn++;
   } 

   TCanvas *c1=new TCanvas("c1","c1",1500,1200);
   TMultiGraph *mg=new TMultiGraph();
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(2);
   hratio->SetMarkerSize(2);
   hratio1->SetFillStyle(3002);
   hratio1->SetFillColor(kAzure-3);
   hr1_norm->SetFillStyle(3004);
   hr1_norm->SetFillColor(1);
   hr_norm->SetFillStyle(3002);
   hr_norm->SetFillColor(1);
   hratio2->SetLineStyle(9);
   hratio2->SetLineColor(kBlue);
   hratio2->SetLineWidth(2);
   hratio3->SetFillStyle(3001);
   hratio3->SetFillColor(kAzure);
   hratio5->SetLineStyle(1);
   hratio5->SetLineColor(8);
   hratio5->SetLineWidth(2);
//   hratio4->SetFillStyle(3005);
//   hratio4->SetFillColor(kBlue-7);
   hratio4->SetLineStyle(7);
   hratio4->SetLineColor(kRed);
   hratio4->SetLineWidth(2);
   mg->Add(hr1_norm,"E3");
   mg->Add(hr_norm,"E3");
   mg->Add(hratio1,"E3");
   mg->Add(hratio3,"E3");
   mg->Add(hratio2,"L");
   mg->Add(hratio4,"L");
 //  mg->Add(hratio3,"L");
   mg->Add(hratio5,"L");
   mg->Add(hratio,"P");
   mg->Draw("A");
   mg->SetTitle(";Bjorken x;#sigma({}^{2}H)/#sigma({}^{1}H)");
   mg->GetYaxis()->SetRangeUser(1.3,1.9);

   auto leg1=new TLegend(0.4,0.75,0.9,0.9);
   leg1->SetNColumns(2);
   leg1->AddEntry(hratio,"#scale[1.5]{MARATHON}","P");
   leg1->AddEntry(hratio1,"#scale[1.5]{Whitlow}","F");
   leg1->AddEntry(hratio2,"#scale[1.5]{Bodek}","L");
   leg1->AddEntry(hratio3,"#scale[1.5]{CJ15}","F");
   leg1->AddEntry(hratio4,"#scale[1.5]{NMC}","L");
   leg1->AddEntry(hratio5,"#scale[1.5]{K&P}","L");
   leg1->AddEntry(hr_norm,"#scale[1.5]{MARATHON norm. uncer.}","F");
   leg1->AddEntry(hr1_norm,"#scale[1.5]{Whitlow norm. uncer.}","F");
   leg1->SetMargin(0.4);
   leg1->Draw();

   c1->Print("Dp_final.pdf");
}
