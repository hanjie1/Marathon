#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp()
{
     Double_t H1_x[5][MAXBIN],H1_Q2[5][MAXBIN],H1_Yield[5][MAXBIN],H1_Yerr[5][MAXBIN];
     Double_t H1_x1[5][MAXBIN],H1_Q21[5][MAXBIN],H1_Yield1[5][MAXBIN],H1_Yerr1[5][MAXBIN];
     Double_t D2_x[5][MAXBIN],D2_Q2[5][MAXBIN],D2_Yield[5][MAXBIN],D2_Yerr[5][MAXBIN];
     Double_t D2_x1[5][MAXBIN],D2_Q21[5][MAXBIN],D2_Yield1[5][MAXBIN],D2_Yerr1[5][MAXBIN];
     Double_t H1_xavg[5][MAXBIN],D2_xavg[5][MAXBIN];
     Double_t H1_xavg1[5][MAXBIN],D2_xavg1[5][MAXBIN];
     Double_t Dp_ratio[5][MAXBIN],Dp_err[5][MAXBIN];
     Double_t Dp_ratio1[5][MAXBIN],Dp_err1[5][MAXBIN];
     Double_t newx[5][10];

     for(int ii=0;ii<5;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_Yield[ii][jj]=0.0; H1_Yerr[ii][jj]=0.0;
             H1_x1[ii][jj]=0.0; H1_Q21[ii][jj]=0.0; H1_Yield1[ii][jj]=0.0; H1_Yerr1[ii][jj]=0.0;
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
             D2_x1[ii][jj]=0.0; D2_Q21[ii][jj]=0.0; D2_Yield1[ii][jj]=0.0; D2_Yerr1[ii][jj]=0.0;
             H1_xavg[ii][jj]=0.0;D2_xavg[ii][jj]=0.0;
             H1_xavg1[ii][jj]=0.0;D2_xavg1[ii][jj]=0.0;
             Dp_ratio[ii][jj]=0.0;Dp_err[ii][jj]=0.0;
             Dp_ratio1[ii][jj]=0.0;Dp_err1[ii][jj]=0.0;
             newx[ii][jj]=0.0;
     }}

   TString Yfile;
   Int_t KIN=0;
   while(KIN<5){
       Yfile=Form("vz009/H1_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H1_x,H1_xavg,H1_Q2,H1_Yield,H1_Yerr); 
       Yfile=Form("vz009/D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x,D2_xavg,D2_Q2,D2_Yield,D2_Yerr); 
       Yfile=Form("vz009_new_25per/H1_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H1_x1,H1_xavg1,H1_Q21,H1_Yield1,H1_Yerr1); 
       Yfile=Form("vz009_new_25per/D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x1,D2_xavg1,D2_Q21,D2_Yield1,D2_Yerr1); 
       KIN+=1;
   }

  for(int ii=0;ii<5;ii++){
      for(int jj=0;jj<10;jj++){
          if(H1_x1[ii][jj+1]!=0)
	     newx[ii][jj]=(H1_x1[ii][jj]+H1_x1[ii][jj+1])/2.0;
          else{
               if(H1_x1[ii][jj]==0)break;
               else {
		    newx[ii][jj]=H1_x1[ii][jj]+(H1_x1[ii][1]-H1_x1[ii][0])/2.0;
               }
          }
          if(ii==0)cout<<newx[ii][jj]<<", ";
      }
  }
  cout<<endl;

  ofstream myfile;
  myfile.open("DP_ratio.txt");
  myfile<<"-------- old bin ------------"<<endl;
  TGraphErrors *gDpRaw[5];
  for(int ii=0;ii<5;ii++){
      gDpRaw[ii]=new TGraphErrors();
      int nn=0;
      myfile<<"kin "<<ii<<endl;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||H1_x[ii][jj]==0)continue;
          Dp_ratio[ii][jj]=D2_Yield[ii][jj]/H1_Yield[ii][jj];
	  Dp_err[ii][jj]=Dp_ratio[ii][jj]*sqrt(pow(D2_Yerr[ii][jj]/D2_Yield[ii][jj],2)+pow(H1_Yerr[ii][jj]/H1_Yield[ii][jj],2));
          if(Dp_err[ii][jj]>0.1)continue;
          gDpRaw[ii]->SetPoint(nn,D2_xavg[ii][jj],Dp_ratio[ii][jj]);
          gDpRaw[ii]->SetPointError(nn,0.0,Dp_err[ii][jj]);
          myfile<<D2_x[ii][jj]+0.01<<"  "<<D2_xavg[ii][jj]<<"  "<<Dp_ratio[ii][jj]<<"  "<<Dp_err[ii][jj]<<endl;
          nn++;
      }
  }
  myfile.close();
    
  ofstream myfile1;
  myfile1.open("DP_ratio_new.txt");
  myfile1<<"-------- new bin ------------"<<endl;
  TGraphErrors *gDpRaw1[5];
  for(int ii=0;ii<5;ii++){
      gDpRaw1[ii]=new TGraphErrors();
      myfile1<<"kin "<<ii<<endl;
      int nn=0;
      for(int jj=0;jj<10;jj++){
          if(newx[ii][jj]==0)continue;
          Dp_ratio1[ii][jj]=D2_Yield1[ii][jj]/H1_Yield1[ii][jj];
	  Dp_err1[ii][jj]=Dp_ratio1[ii][jj]*sqrt(pow(D2_Yerr1[ii][jj]/D2_Yield1[ii][jj],2)+pow(H1_Yerr1[ii][jj]/H1_Yield1[ii][jj],2));
          gDpRaw1[ii]->SetPoint(nn,newx[ii][jj],Dp_ratio1[ii][jj]);
          gDpRaw1[ii]->SetPointError(nn,0.0,Dp_err1[ii][jj]);
          myfile1<<newx[ii][jj]<<"  "<<D2_xavg1[ii][jj]<<"  "<<Dp_ratio1[ii][jj]<<"  "<<Dp_err1[ii][jj]<<endl;
          nn++;
      }
  }
  myfile1.close();
  int color[5]={1,2,4,6,46};
  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<5;ii++){
      gDpRaw[ii]->SetMarkerStyle(8);
      gDpRaw[ii]->SetMarkerColor(color[ii]);
      gDpRaw1[ii]->SetMarkerStyle(22);
      gDpRaw1[ii]->SetMarkerColor(color[ii]);
      mg1->Add(gDpRaw[ii]);
      mg1->Add(gDpRaw1[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("D/p Raw Yield ratio;x;D/p");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  for(int ii=0;ii<5;ii++){
      leg1->AddEntry(gDpRaw[ii],Form("D/p kin%d",ii),"P");
      leg1->AddEntry(gDpRaw1[ii],Form("new D/p kin%d",ii),"P");
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  gDpRaw[0]->SetMarkerStyle(8);
  gDpRaw[0]->SetMarkerColor(1);
  gDpRaw1[0]->SetMarkerStyle(22);
  gDpRaw1[0]->SetMarkerColor(2);
  mg2->Add(gDpRaw[0]);
  mg2->Add(gDpRaw1[0]);
  mg2->Draw("AP");
  mg2->SetTitle("D/p Raw Yield ratio Kin0;x;D/p");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  leg2->AddEntry(gDpRaw[0],"D/p kin0","P");
  leg2->AddEntry(gDpRaw1[0],"new D/p kin0","P");
  leg2->Draw();

  TCanvas *c3=new TCanvas("c3");
  TMultiGraph *mg3=new TMultiGraph();
  gDpRaw[1]->SetMarkerStyle(8);
  gDpRaw[1]->SetMarkerColor(1);
  gDpRaw1[1]->SetMarkerStyle(22);
  gDpRaw1[1]->SetMarkerColor(2);
  mg3->Add(gDpRaw[1]);
  mg3->Add(gDpRaw1[1]);
  mg3->Draw("AP");
  mg3->SetTitle("D/p Raw Yield ratio Kin1;x;D/p");

  auto leg3=new TLegend(0.7,0.6,0.85,0.85);
  leg3->AddEntry(gDpRaw[1],"D/p kin1","P");
  leg3->AddEntry(gDpRaw1[1],"new D/p kin1","P");
  leg3->Draw();

  TCanvas *c4=new TCanvas("c4");
  TMultiGraph *mg4=new TMultiGraph();
  gDpRaw[2]->SetMarkerStyle(8);
  gDpRaw[2]->SetMarkerColor(1);
  gDpRaw1[2]->SetMarkerStyle(22);
  gDpRaw1[2]->SetMarkerColor(2);
  mg4->Add(gDpRaw[2]);
  mg4->Add(gDpRaw1[2]);
  mg4->Draw("AP");
  mg4->SetTitle("D/p Raw Yield ratio Kin2;x;D/p");

  auto leg4=new TLegend(0.7,0.6,0.85,0.85);
  leg4->AddEntry(gDpRaw[2],"D/p kin2","P");
  leg4->AddEntry(gDpRaw1[2],"new D/p kin2","P");
  leg4->Draw();


}
