#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3_yield()
{
     Double_t H3_x[11][MAXBIN],H3_xavg[11][MAXBIN],H3_Q2[11][MAXBIN],H3_Yield[11][MAXBIN],H3_Yerr[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Yield[ii][jj]=0.0; H3_Yerr[ii][jj]=0.0;
	     H3_xavg[ii][jj]=0.0;
     }}

   TString Yfile;
   Int_t KIN=0;
   while(KIN<16){
       Yfile=Form("H3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H3_x,H3_xavg,H3_Q2,H3_Yield,H3_Yerr); 
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gH3Raw[11];
  TGraph *gH3xbj[11];
  for(int ii=0;ii<11;ii++){
      gH3xbj[ii]=new TGraph();
      gH3Raw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(H3_xavg[ii][jj]==0)continue;
          gH3Raw[ii]->SetPoint(nn,H3_x[ii][jj],H3_Yield[ii][jj]);
          gH3Raw[ii]->SetPointError(nn,0.0,H3_Yerr[ii][jj]);
          Double_t tmpx=H3_x[ii][jj]+dBin/2.0;
          Double_t tmpxavg=(int)((H3_xavg[ii][jj]+0.0005)/0.001)*0.001;
//  cout<<tmpx<<"  "<<tmpxavg<<endl;
          gH3xbj[ii]->SetPoint(nn,tmpx,tmpxavg-tmpx);
          nn++;
      }
  }

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gH3Raw[ii]->SetMarkerStyle(8);
      if(ii==9)gH3Raw[ii]->SetMarkerColor(30);
      else gH3Raw[ii]->SetMarkerColor(ii+1);
      mg1->Add(gH3Raw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("H3 Data Yield;x;Yield");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<=5)leg1->AddEntry(gH3Raw[ii],Form("H3 kin%d",ii),"P");
      else {leg1->AddEntry(gH3Raw[ii],Form("H3 kin%d",ii+nn),"P");
            nn++;
      }
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gH3xbj[ii]->SetMarkerStyle(8);
      if(ii==9)gH3xbj[ii]->SetMarkerColor(30);
      else gH3xbj[ii]->SetMarkerColor(ii+1);
      mg2->Add(gH3xbj[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("H3 Data Yield;x;Yield");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<=5)leg2->AddEntry(gH3xbj[ii],Form("H3 kin%d",ii),"P");
      else {leg2->AddEntry(gH3xbj[ii],Form("H3 kin%d",ii+nn),"P");
            nn++;
      }
  }
  leg2->Draw();



}
