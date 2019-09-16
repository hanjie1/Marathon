#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp()
{
   Double_t x[MAXBIN]={0.0},Ratio[MAXBIN]={0.0},Rerr[MAXBIN]={0.0};

   TString Rfile="newbin/Dp_final.dat";
   int nbin=ReadYield(Rfile,x,Ratio,Rerr); 

   TGraphErrors *hratio=new TGraphErrors();
   
   for(int ii=0;ii<nbin;ii++){
	hratio->SetPoint(ii,x[ii],Ratio[ii]);
	hratio->SetPointError(ii,0,Rerr[ii]);
   } 

   TCanvas *c1=new TCanvas("c1","c1",1500,1200);
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(2);
   hratio->Draw("AP");
   hratio->SetTitle(";Bjorken x;#sigma({}^{2}H)/#sigma({}^{1}He)");
 //  hratio->GetYaxis()->SetRangeUser(1.22,1.5);


   c1->Print("Dp_final.pdf");
}
