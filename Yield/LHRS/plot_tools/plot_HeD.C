#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_HeD()
{
     Double_t He3_x[11][MAXBIN],He3_Q2[11][MAXBIN],He3_Yield[11][MAXBIN],He3_Yerr[11][MAXBIN];
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_Yield[11][MAXBIN],D2_Yerr[11][MAXBIN];
     Double_t He3_xavg[11][MAXBIN],D2_xavg[11][MAXBIN];
     Double_t He3D_ratio[11][MAXBIN],He3D_err[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Yield[ii][jj]=0.0; He3_Yerr[ii][jj]=0.0;
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
             He3_xavg[ii][jj]=0.0;D2_xavg[ii][jj]=0.0;
             He3D_ratio[ii][jj]=0.0;He3D_err[ii][jj]=0.0;
     }}

   TString Yfile;
   Int_t KIN=0;
   while(KIN<16){
       Yfile=Form("He3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,He3_x,He3_xavg,He3_Q2,He3_Yield,He3_Yerr); 
       Yfile=Form("D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x,D2_xavg,D2_Q2,D2_Yield,D2_Yerr); 
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gHe3DRaw[11];
  for(int ii=0;ii<11;ii++){
      gHe3DRaw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||He3_x[ii][jj]==0)continue;
          //cout<<D2_x[ii][jj]<<"  "<<He3_x[ii][jj]<<"  "<<He3_Yield[ii][jj]<<"  "<<D2_Yield[ii][jj]<<endl;
          He3D_ratio[ii][jj]=He3_Yield[ii][jj]/D2_Yield[ii][jj];
	  He3D_err[ii][jj]=He3D_ratio[ii][jj]*sqrt(pow(D2_Yerr[ii][jj]/D2_Yield[ii][jj],2)+pow(He3_Yerr[ii][jj]/He3_Yield[ii][jj],2));
          //if(He3D_err[ii][jj]>0.1)continue;
          gHe3DRaw[ii]->SetPoint(nn,D2_x[ii][jj],He3D_ratio[ii][jj]);
          gHe3DRaw[ii]->SetPointError(nn,0.0,He3D_err[ii][jj]);
          nn++;
      }
  }
    
  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gHe3DRaw[ii]->SetMarkerStyle(8);
      if(ii<9) gHe3DRaw[ii]->SetMarkerColor(ii+1);
      else gHe3DRaw[ii]->SetMarkerColor(ii+2);
      mg1->Add(gHe3DRaw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("He3/D Raw Yield ratio;x;He3/D");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<6)leg1->AddEntry(gHe3DRaw[ii],Form("He3/D kin%d",ii),"P");
      else {leg1->AddEntry(gHe3DRaw[ii],Form("He3/D kin%d",ii+nn),"P");
             nn++;
           }
  }
  leg1->Draw();


}
