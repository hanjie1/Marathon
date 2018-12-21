#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp()
{
     Double_t H1_x[5][MAXBIN],H1_Q2[5][MAXBIN],H1_Yield[5][MAXBIN],H1_Yerr[5][MAXBIN];
     Double_t D2_x[5][MAXBIN],D2_Q2[5][MAXBIN],D2_Yield[5][MAXBIN],D2_Yerr[5][MAXBIN];
     Double_t H1_xavg[5][MAXBIN],D2_xavg[5][MAXBIN];
     Double_t Dp_ratio[5][MAXBIN],Dp_err[5][MAXBIN];

     for(int ii=0;ii<5;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_Yield[ii][jj]=0.0; H1_Yerr[ii][jj]=0.0;
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
             H1_xavg[ii][jj]=0.0;D2_xavg[ii][jj]=0.0;
             Dp_ratio[ii][jj]=0.0;Dp_err[ii][jj]=0.0;
     }}

   TString Yfile;
   Int_t KIN=0;
   while(KIN<5){
       Yfile=Form("H1_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H1_x,H1_xavg,H1_Q2,H1_Yield,H1_Yerr); 
       Yfile=Form("D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x,D2_xavg,D2_Q2,D2_Yield,D2_Yerr); 
       KIN+=1;
   }
  TGraphErrors *gDpRaw[5];
  for(int ii=0;ii<5;ii++){
      gDpRaw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||H1_x[ii][jj]==0)continue;
          Dp_ratio[ii][jj]=D2_Yield[ii][jj]/H1_Yield[ii][jj];
	  Dp_err[ii][jj]=Dp_ratio[ii][jj]*sqrt(pow(D2_Yerr[ii][jj]/D2_Yield[ii][jj],2)+pow(H1_Yerr[ii][jj]/H1_Yield[ii][jj],2));
          if(Dp_err[ii][jj]>0.1)continue;
          gDpRaw[ii]->SetPoint(nn,D2_x[ii][jj],Dp_ratio[ii][jj]);
          gDpRaw[ii]->SetPointError(nn,0.0,Dp_err[ii][jj]);
          nn++;
      }
  }
    

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<5;ii++){
      gDpRaw[ii]->SetMarkerStyle(8);
      gDpRaw[ii]->SetMarkerColor(ii+1);
      mg1->Add(gDpRaw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("D/p Raw Yield ratio;x;D/p");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  for(int ii=0;ii<5;ii++){
      leg1->AddEntry(gDpRaw[ii],Form("D/p kin%d",ii),"P");
  }
  leg1->Draw();


}
