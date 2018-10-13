#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He()
{
     Double_t H3_x[11][MAXBIN],H3_Q2[11][MAXBIN],H3_Yield[11][MAXBIN],H3_Yerr[11][MAXBIN];
     Double_t He3_x[11][MAXBIN],He3_Q2[11][MAXBIN],He3_Yield[11][MAXBIN],He3_Yerr[11][MAXBIN];
     Double_t H3_xavg[11][MAXBIN],He3_xavg[11][MAXBIN];
     Double_t H3He_ratio[11][MAXBIN],H3He_err[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Yield[ii][jj]=0.0; H3_Yerr[ii][jj]=0.0;
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Yield[ii][jj]=0.0; He3_Yerr[ii][jj]=0.0;
             H3_xavg[ii][jj]=0.0;He3_xavg[ii][jj]=0.0;
             H3He_ratio[ii][jj]=0.0;H3He_err[ii][jj]=0.0;
     }}

   TString Yfile;
   Int_t KIN=0;
   while(KIN<16){
       Yfile=Form("H3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H3_x,H3_xavg,H3_Q2,H3_Yield,H3_Yerr); 
       Yfile=Form("He3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,He3_x,He3_xavg,He3_Q2,He3_Yield,He3_Yerr); 
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gH3HeRaw[11];
  for(int ii=0;ii<11;ii++){
      gH3HeRaw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(He3_x[ii][jj]==0||H3_x[ii][jj]==0)continue;
          H3He_ratio[ii][jj]=H3_Yield[ii][jj]/He3_Yield[ii][jj];
	  H3He_err[ii][jj]=H3He_ratio[ii][jj]*sqrt(pow(He3_Yerr[ii][jj]/He3_Yield[ii][jj],2)+pow(H3_Yerr[ii][jj]/H3_Yield[ii][jj],2));
          //if(H3He_err[ii][jj]>0.1)continue;
          gH3HeRaw[ii]->SetPoint(nn,He3_x[ii][jj],H3He_ratio[ii][jj]);
          gH3HeRaw[ii]->SetPointError(nn,0.0,H3He_err[ii][jj]);
          nn++;
      }
  }
    
  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gH3HeRaw[ii]->SetMarkerStyle(8);
      if(ii<9) gH3HeRaw[ii]->SetMarkerColor(ii+1);
      else gH3HeRaw[ii]->SetMarkerColor(ii+2);
      mg1->Add(gH3HeRaw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("H3/He Raw Yield ratio;x;H3/He");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<6)leg1->AddEntry(gH3HeRaw[ii],Form("H3/He kin%d",ii),"P");
      else {leg1->AddEntry(gH3HeRaw[ii],Form("H3/He kin%d",ii+nn),"P");
             nn++;
           }
  }
  leg1->Draw();


}
