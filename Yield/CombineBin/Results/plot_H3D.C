#include <fstream>
#include "ReadFile.h"
using namespace std;

Double_t KP_EMC(Double_t x){
	Double_t ratio=1.0;
	ratio=1.07251-1.61648*x+12.3626*pow(x,2)-65.5932*pow(x,3)+213.311*pow(x,4)
              -423.943*pow(x,5)+499.994*pow(x,6)-321.304*pow(x,7)+87.0596*pow(x,8);
	return ratio;
}

void plot_H3D()
{
   Double_t x1[MAXBIN]={0.0},Ratio1[MAXBIN]={0.0},Rerr1[MAXBIN]={0.0};
   Double_t x2[MAXBIN]={0.0},Ratio2[MAXBIN]={0.0},Rerr2[MAXBIN]={0.0};

   auto f1=new TF1("f1","KP_EMC(x)",0.16,0.84);

   TString Rfile="newbin/H3D_final.dat";
   int nbin1=ReadYield(Rfile,x1,Ratio1,Rerr1); 
   Rfile="bin003/H3D_final.dat";
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
   mg1->SetTitle("H3/D2 yield ratio;xbj;");

   f1->SetLineColor(4);
   f1->SetLineStyle(9);
   f1->Draw("same");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio1,"newbin","P");
   leg1->AddEntry(hratio2,"bin003","P");
   leg1->Draw();

}
