#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He()
{
   Double_t x1[MAXBIN]={0.0},Ratio1[MAXBIN]={0.0},Rerr1[MAXBIN]={0.0};
   Double_t x2[MAXBIN]={0.0},Ratio2[MAXBIN]={0.0},Rerr2[MAXBIN]={0.0};

   TString Rfile="newbin/H3He_final.dat";
   int nbin1=ReadYield(Rfile,x1,Ratio1,Rerr1); 
   Rfile="bin003/H3He_final.dat";
   int nbin2=ReadYield(Rfile,x2,Ratio2,Rerr2); 

   TGraphErrors *hratio1=new TGraphErrors();
   TGraphErrors *hratio2=new TGraphErrors();
   
   for(int ii=0;ii<nbin1;ii++){
	hratio1->SetPoint(ii,x1[ii],Ratio1[ii]);
	hratio1->SetPointError(ii,0,Rerr1[ii]);
   } 

   for(int ii=0;ii<nbin2;ii++){
	hratio2->SetPoint(ii,x2[ii],Ratio2[ii]);
	hratio2->SetPointError(ii,0,Rerr2[ii]);
   } 



   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   TMultiGraph *mg1=new TMultiGraph();
   hratio1->SetMarkerStyle(8);
   hratio1->SetMarkerColor(2);
   hratio2->SetMarkerStyle(8);
   hratio2->SetMarkerColor(1);
   mg1->Add(hratio1);
   mg1->Add(hratio2);
   mg1->Draw("AP");
   mg1->SetTitle("H3/He3 yield ratio;xbj;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio1,"newbin","P");
   leg1->AddEntry(hratio2,"bin003","P");
   leg1->Draw();

}
