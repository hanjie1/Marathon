#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He()
{
     Double_t He3_x[11][MAXBIN],He3_Q2[11][MAXBIN],He3_Yield[11][MAXBIN],He3_Yerr[11][MAXBIN];
     Double_t He3_x1[11][MAXBIN],He3_Q21[11][MAXBIN],He3_Yield1[11][MAXBIN],He3_Yerr1[11][MAXBIN];
     Double_t H3_x[11][MAXBIN],H3_Q2[11][MAXBIN],H3_Yield[11][MAXBIN],H3_Yerr[11][MAXBIN];
     Double_t H3_x1[11][MAXBIN],H3_Q21[11][MAXBIN],H3_Yield1[11][MAXBIN],H3_Yerr1[11][MAXBIN];
     Double_t He3_xavg[11][MAXBIN],H3_xavg[11][MAXBIN];
     Double_t He3_xavg1[11][MAXBIN],H3_xavg1[11][MAXBIN];
     Double_t H3He_ratio[11][MAXBIN],H3He_err[11][MAXBIN];
     Double_t H3He_ratio1[11][MAXBIN],H3He_err1[11][MAXBIN];
     Double_t newx[11][10];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Yield[ii][jj]=0.0; He3_Yerr[ii][jj]=0.0;
             He3_x1[ii][jj]=0.0; He3_Q21[ii][jj]=0.0; He3_Yield1[ii][jj]=0.0; He3_Yerr1[ii][jj]=0.0;
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Yield[ii][jj]=0.0; H3_Yerr[ii][jj]=0.0;
             H3_x1[ii][jj]=0.0; H3_Q21[ii][jj]=0.0; H3_Yield1[ii][jj]=0.0; H3_Yerr1[ii][jj]=0.0;
             He3_xavg[ii][jj]=0.0;H3_xavg[ii][jj]=0.0;
             He3_xavg1[ii][jj]=0.0;H3_xavg1[ii][jj]=0.0;
             H3He_ratio[ii][jj]=0.0;H3He_err[ii][jj]=0.0;
             H3He_ratio1[ii][jj]=0.0;H3He_err1[ii][jj]=0.0;
             newx[ii][jj]=0.0;
     }}

   TString Yfile;
   Int_t KIN=0;
   while(KIN<16){
       Yfile=Form("vz009/He3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,He3_x,He3_xavg,He3_Q2,He3_Yield,He3_Yerr); 
       Yfile=Form("vz009/H3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H3_x,H3_xavg,H3_Q2,H3_Yield,H3_Yerr); 
       Yfile=Form("vz009_new_25per/He3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,He3_x1,He3_xavg1,He3_Q21,He3_Yield1,He3_Yerr1); 
       Yfile=Form("vz009_new_25per/H3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H3_x1,H3_xavg1,H3_Q21,H3_Yield1,H3_Yerr1); 
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  for(int ii=0;ii<11;ii++){
      for(int jj=0;jj<10;jj++){
          if(He3_x1[ii][jj+1]!=0)
	     newx[ii][jj]=(He3_x1[ii][jj]+He3_x1[ii][jj+1])/2.0;
          else{
               if(He3_x1[ii][jj]==0)break;
               else {
		    newx[ii][jj]=He3_x1[ii][jj]+(He3_x1[ii][1]-He3_x1[ii][0])/2.0;
               }
          }
          if(ii==0)cout<<newx[ii][jj]<<", ";
      }
  }
  cout<<endl;

  ofstream myfile;
  myfile.open("H3He_ratio.txt");
  myfile<<"-------- old bin ------------"<<endl;

  TGraphErrors *gH3HeRaw[11];
  for(int ii=0;ii<11;ii++){
      gH3HeRaw[ii]=new TGraphErrors();
      int nn=0;
      myfile<<"kin "<<ii<<endl;
      for(int jj=0;jj<MAXBIN;jj++){
          if(H3_x[ii][jj]==0||He3_x[ii][jj]==0)continue;
          H3He_ratio[ii][jj]=H3_Yield[ii][jj]/He3_Yield[ii][jj];
	  H3He_err[ii][jj]=H3He_ratio[ii][jj]*sqrt(pow(H3_Yerr[ii][jj]/H3_Yield[ii][jj],2)+pow(He3_Yerr[ii][jj]/He3_Yield[ii][jj],2));
          if(H3He_err[ii][jj]>0.1)continue;
          gH3HeRaw[ii]->SetPoint(nn,H3_xavg[ii][jj],H3He_ratio[ii][jj]);
          gH3HeRaw[ii]->SetPointError(nn,0.0,H3He_err[ii][jj]);
          myfile<<H3_x[ii][jj]+0.01<<"  "<<H3_xavg[ii][jj]<<"  "<<H3He_ratio[ii][jj]<<"  "<<H3He_err[ii][jj]<<endl;
          nn++;
      }
  }
  myfile.close();    

  ofstream myfile1;
  myfile1.open("H3He_ratio_new.txt");
  myfile1<<"-------- new bin ------------"<<endl;
  TGraphErrors *gH3HeRaw1[11];
  for(int ii=0;ii<11;ii++){
      gH3HeRaw1[ii]=new TGraphErrors();
      myfile1<<"kin "<<ii<<endl;
      int nn=0;
      for(int jj=0;jj<10;jj++){
          if(newx[ii][jj]==0)continue;
          H3He_ratio1[ii][jj]=H3_Yield1[ii][jj]/He3_Yield1[ii][jj];
	  H3He_err1[ii][jj]=H3He_ratio1[ii][jj]*sqrt(pow(H3_Yerr1[ii][jj]/H3_Yield1[ii][jj],2)+pow(He3_Yerr1[ii][jj]/He3_Yield1[ii][jj],2));
          gH3HeRaw1[ii]->SetPoint(nn,newx[ii][jj],H3He_ratio1[ii][jj]);
          gH3HeRaw1[ii]->SetPointError(nn,0.0,H3He_err1[ii][jj]);
          myfile1<<newx[ii][jj]<<"  "<<H3_xavg1[ii][jj]<<"  "<<H3He_ratio1[ii][jj]<<"  "<<H3He_err1[ii][jj]<<endl;
          nn++;
      }
  }
  myfile1.close();    


  int color[11]={1,30,3,4,6,7,8,9,2,36,46};
  int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gH3HeRaw[ii]->SetMarkerStyle(8);
      gH3HeRaw[ii]->SetMarkerColor(color[ii]);
      gH3HeRaw1[ii]->SetMarkerStyle(22);
      gH3HeRaw1[ii]->SetMarkerColor(color[ii]);
      mg1->Add(gH3HeRaw[ii]);
      mg1->Add(gH3HeRaw1[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("H3/He3 Raw Yield ratio;x;H3/He3");

  auto leg1=new TLegend(0.7,0.6,0.811,0.811);
  for(int ii=0;ii<11;ii++){
      leg1->AddEntry(gH3HeRaw[ii],Form("H3/He3 kin%d",kin[ii]),"P");
      leg1->AddEntry(gH3HeRaw1[ii],Form("new H3/He3 kin%d",kin[ii]),"P");
  }
  leg1->Draw();
/*
  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  gH3HeRaw[0]->SetMarkerStyle(8);
  gH3HeRaw[0]->SetMarkerColor(1);
  gH3HeRaw1[0]->SetMarkerStyle(22);
  gH3HeRaw1[0]->SetMarkerColor(2);
  mg2->Add(gH3HeRaw[0]);
  mg2->Add(gH3HeRaw1[0]);
  mg2->Draw("AP");
  mg2->SetTitle("D/p Raw Yield ratio Kin0;x;D/p");

  auto leg2=new TLegend(0.7,0.6,0.811,0.811);
  leg2->AddEntry(gH3HeRaw[0],"D/p kin0","P");
  leg2->AddEntry(gH3HeRaw1[0],"new D/p kin0","P");
  leg2->Draw();

  TCanvas *c3=new TCanvas("c3");
  TMultiGraph *mg3=new TMultiGraph();
  gH3HeRaw[1]->SetMarkerStyle(8);
  gH3HeRaw[1]->SetMarkerColor(1);
  gH3HeRaw1[1]->SetMarkerStyle(22);
  gH3HeRaw1[1]->SetMarkerColor(2);
  mg3->Add(gH3HeRaw[1]);
  mg3->Add(gH3HeRaw1[1]);
  mg3->Draw("AP");
  mg3->SetTitle("D/p Raw Yield ratio Kin1;x;D/p");

  auto leg3=new TLegend(0.7,0.6,0.811,0.811);
  leg3->AddEntry(gH3HeRaw[1],"D/p kin1","P");
  leg3->AddEntry(gH3HeRaw1[1],"new D/p kin1","P");
  leg3->Draw();

  TCanvas *c4=new TCanvas("c4");
  TMultiGraph *mg4=new TMultiGraph();
  gH3HeRaw[2]->SetMarkerStyle(8);
  gH3HeRaw[2]->SetMarkerColor(1);
  gH3HeRaw1[2]->SetMarkerStyle(22);
  gH3HeRaw1[2]->SetMarkerColor(2);
  mg4->Add(gH3HeRaw[2]);
  mg4->Add(gH3HeRaw1[2]);
  mg4->Draw("AP");
  mg4->SetTitle("D/p Raw Yield ratio Kin2;x;D/p");

  auto leg4=new TLegend(0.7,0.6,0.811,0.811);
  leg4->AddEntry(gH3HeRaw[2],"D/p kin2","P");
  leg4->AddEntry(gH3HeRaw1[2],"new D/p kin2","P");
  leg4->Draw();
*/

}
