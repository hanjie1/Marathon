#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp_data()
{
   Double_t x[MAXBIN]={0.0},Ratio[MAXBIN]={0.0},Rerr[MAXBIN]={0.0};
   Double_t nexp[2000]={0.0},x1[2000]={0.0},Q2[2000]={0.0},Ratio1[2000]={0.0},Rerr1[2000]={0.0},Rnorm1[2000]={0.0};

   TString Rfile="newbin/Dp_final.dat";
   int nbin=ReadYield(Rfile,x,Ratio,Rerr); 
   Rfile="Other_Data/F2_Dp_GLOBAL_DATA.out";
   int nbin1=ReadGlobal(Rfile,nexp,x1,Q2,Ratio1,Rerr1,Rnorm1); 

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hr_norm=new TGraphErrors(); //MARATHON normalization error
   TGraphErrors *hR[9];
   for(int ii=0;ii<9;ii++){
       hR[ii]=new TGraphErrors();
   }
   
   for(int ii=0;ii<nbin;ii++){
	hratio->SetPoint(ii,x[ii],Ratio[ii]);
	hratio->SetPointError(ii,0,Rerr[ii]);
	hr_norm->SetPoint(ii,x[ii],1.4);
	hr_norm->SetPointError(ii,0,0.79/100*Ratio[ii]);
   } 

   int nn[9]={0};
   for(int ii=0;ii<nbin1;ii++){
	if(Q2[ii]<3.0 || Q2[ii]>15)continue;
	if(x1[ii]<0.15 || x1[ii]>0.4)continue;
        int nnexp=nexp[ii]-1;
	hR[nnexp]->SetPoint(nn[nnexp],x1[ii],Ratio1[ii]*2.0);
	hR[nnexp]->SetPointError(nn[nnexp],0,Rerr1[ii]*Ratio1[ii]*2.0);
	nn[nnexp]++;
   } 


   TCanvas *c1=new TCanvas("c1","c1",1500,1200);
   TMultiGraph *mg=new TMultiGraph();
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(2);
   hratio->SetMarkerSize(2);

   
   for(int ii=0;ii<9;ii++){
       hR[ii]->SetMarkerStyle(8);
       hR[ii]->SetMarkerColor(kViolet+ii);
       mg->Add(hR[ii]);
   }
   mg->Add(hratio,"P");
   mg->Draw("AP");
   mg->SetTitle(";Bjorken x;#sigma({}^{2}H)/#sigma({}^{1}He)");
   //mg->GetYaxis()->SetRangeUser(1.3,1.9);

   auto leg1=new TLegend(0.55,0.75,0.9,0.9);
   leg1->SetNColumns(3);
   leg1->AddEntry(hratio,"#scale[1]{MARATHON}","P");
   leg1->AddEntry(hR[0],"#scale[1]{E49a}","P");
   leg1->AddEntry(hR[1],"#scale[1]{E49b}","P");
   leg1->AddEntry(hR[2],"#scale[1]{E61}","P");
   leg1->AddEntry(hR[3],"#scale[1]{E87}","P");
   leg1->AddEntry(hR[4],"#scale[1]{E89a}","P");
   leg1->AddEntry(hR[5],"#scale[1]{E89b}","P");
   leg1->AddEntry(hR[6],"#scale[1]{NMC}","P");
   leg1->AddEntry(hR[7],"#scale[1]{JLab}","P");
   leg1->AddEntry(hR[8],"#scale[1]{BCDMS}","P");
   leg1->Draw();

   c1->Print("Dp_data.pdf");
}
