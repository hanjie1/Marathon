#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_D2_yield()
{
     Double_t D2_x[11][MAXBIN],D2_xavg[11][MAXBIN],D2_Q2[11][MAXBIN],D2_Yield[11][MAXBIN],D2_Yerr[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
	     D2_xavg[ii][jj]=0.0;
     }}

   TString Yfile;
   Int_t KIN=0;
   while(KIN<16){
       Yfile=Form("D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x,D2_xavg,D2_Q2,D2_Yield,D2_Yerr); 
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gD2Raw[11];
  TGraph *gD2xbj[11];
  for(int ii=0;ii<11;ii++){
      gD2xbj[ii]=new TGraph();
      gD2Raw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_xavg[ii][jj]==0)continue;
          gD2Raw[ii]->SetPoint(nn,D2_x[ii][jj],D2_Yield[ii][jj]);
          gD2Raw[ii]->SetPointError(nn,0.0,D2_Yerr[ii][jj]);
          Double_t tmpx=D2_x[ii][jj]+dBin/2.0;
          Double_t tmpxavg=(int)((D2_xavg[ii][jj]+0.0005)/0.001)*0.001;
//  cout<<tmpx<<"  "<<tmpxavg<<endl;
          gD2xbj[ii]->SetPoint(nn,D2_xavg[ii][jj],D2_Yield[ii][jj]);
          nn++;
      }
  }

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gD2Raw[ii]->SetMarkerStyle(8);
      if(ii==9)gD2Raw[ii]->SetMarkerColor(30);
      else gD2Raw[ii]->SetMarkerColor(ii+1);
      mg1->Add(gD2Raw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("D2 Data Yield;x;Yield");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<=5)leg1->AddEntry(gD2Raw[ii],Form("D2 kin%d",ii),"P");
      else {leg1->AddEntry(gD2Raw[ii],Form("D2 kin%d",ii+nn),"P");
            nn++;
      }
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gD2xbj[ii]->SetMarkerStyle(8);
      if(ii==9)gD2xbj[ii]->SetMarkerColor(30);
      else gD2xbj[ii]->SetMarkerColor(ii+1);
      mg2->Add(gD2xbj[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("D2 Data Yield;x;Yield");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<=5)leg2->AddEntry(gD2xbj[ii],Form("D2 kin%d",ii),"P");
      else {leg2->AddEntry(gD2xbj[ii],Form("D2 kin%d",ii+nn),"P");
            nn++;
      }
  }
  leg2->Draw();



}
