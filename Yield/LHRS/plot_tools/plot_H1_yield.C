#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H1_yield()
{
     Double_t H1_x[11][MAXBIN],H1_xavg[11][MAXBIN],H1_Q2[11][MAXBIN],H1_Yield[11][MAXBIN],H1_Yerr[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_Yield[ii][jj]=0.0; H1_Yerr[ii][jj]=0.0;
	     H1_xavg[ii][jj]=0.0;
     }}

   TString Yfile;
   Int_t KIN=1;
   while(KIN<2){
       Yfile=Form("H1_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H1_x,H1_xavg,H1_Q2,H1_Yield,H1_Yerr); 
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gH1Raw[11];
  TGraph *gH1xbj[11];
  for(int ii=1;ii<2;ii++){
      gH1xbj[ii]=new TGraph();
      gH1Raw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(H1_xavg[ii][jj]==0)continue;
          gH1Raw[ii]->SetPoint(nn,H1_x[ii][jj],H1_Yield[ii][jj]);
          gH1Raw[ii]->SetPointError(nn,0.0,H1_Yerr[ii][jj]);
          Double_t tmpx=H1_x[ii][jj]+dBin/2.0;
          Double_t tmpxavg=(int)((H1_xavg[ii][jj]+0.0005)/0.001)*0.001;
//  cout<<tmpx<<"  "<<tmpxavg<<endl;
          gH1xbj[ii]->SetPoint(nn,tmpx,tmpxavg-tmpx);
          nn++;
      }
  }

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=1;ii<2;ii++){
      gH1Raw[ii]->SetMarkerStyle(8);
      if(ii==9)gH1Raw[ii]->SetMarkerColor(30);
      else gH1Raw[ii]->SetMarkerColor(ii+1);
      mg1->Add(gH1Raw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("H1 Data Yield;x;Yield");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=1;
  for(int ii=1;ii<2;ii++){
      if(ii<=5)leg1->AddEntry(gH1Raw[ii],Form("H1 kin%d",ii),"P");
      else {leg1->AddEntry(gH1Raw[ii],Form("H1 kin%d",ii+nn),"P");
            nn++;
      }
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gH1xbj[ii]->SetMarkerStyle(8);
      if(ii==9)gH1xbj[ii]->SetMarkerColor(30);
      else gH1xbj[ii]->SetMarkerColor(ii+1);
      mg2->Add(gH1xbj[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("H1 Data Yield;x;Yield");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  nn=1;
  for(int ii=1;ii<2;ii++){
      if(ii<=5)leg2->AddEntry(gH1xbj[ii],Form("H1 kin%d",ii),"P");
      else {leg2->AddEntry(gH1xbj[ii],Form("H1 kin%d",ii+nn),"P");
            nn++;
      }
  }
  leg2->Draw();



}
