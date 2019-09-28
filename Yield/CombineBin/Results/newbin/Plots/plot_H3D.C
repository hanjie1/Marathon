#include <fstream>
#include "ReadFile.h"
using namespace std;

Double_t KP_EMC(Double_t x){
	Double_t ratio=1.0;
	ratio=1.07251-1.61648*x+12.3626*pow(x,2)-65.5932*pow(x,3)+213.311*pow(x,4)
              -423.943*pow(x,5)+499.994*pow(x,6)-321.304*pow(x,7)+87.0596*pow(x,8);
	ratio=ratio*3.0/2.0;
	return ratio;
}

void plot_H3D()
{
   Double_t x[MAXBIN]={0.0},Ratio[MAXBIN]={0.0},Rerr[MAXBIN]={0.0};

   auto f1=new TF1("f1","KP_EMC(x)",0.16,0.84);

   TString Rfile="newbin/H3D_final.dat";
   int nbin=ReadYield(Rfile,x,Ratio,Rerr); 

   TGraphErrors *hratio=new TGraphErrors();
   
   for(int ii=0;ii<nbin;ii++){
	hratio->SetPoint(ii,x[ii],Ratio[ii]);
	hratio->SetPointError(ii,0,Rerr[ii]);
   } 

   TCanvas *c1=new TCanvas("c1","c1",1500,1200);
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(4);
   hratio->SetMarkerSize(2);
   hratio->Draw("AP");
   hratio->SetTitle(";Bjorken x;#sigma({}^{3}H)/#sigma({}^{2}H)");
   hratio->GetYaxis()->SetRangeUser(1.22,1.5);


   f1->SetLineColor(4);
   f1->SetLineStyle(9);
//   f1->Draw("same");

   c1->Print("H3D_final.pdf");
}
