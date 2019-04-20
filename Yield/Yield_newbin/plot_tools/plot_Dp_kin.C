#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp_kin()
{
     Double_t H1_x[5][MAXBIN],H1_Q2[5][MAXBIN],H1_xavg[5][MAXBIN];
     Double_t D2_x[5][MAXBIN],D2_Q2[5][MAXBIN],D2_xavg[5][MAXBIN];
     Double_t H1_Y[5][MAXBIN],H1_YE[5][MAXBIN]; 
     Double_t D2_Y[5][MAXBIN],D2_YE[5][MAXBIN]; 
     Double_t Dp[5][MAXBIN],Dp_E[5][MAXBIN]; 

     Double_t H1_x1[5][MAXBIN],H1_Q21[5][MAXBIN],H1_xavg1[5][MAXBIN];
     Double_t D2_x1[5][MAXBIN],D2_Q21[5][MAXBIN],D2_xavg1[5][MAXBIN];
     Double_t H1_Y1[5][MAXBIN],H1_YE1[5][MAXBIN]; 
     Double_t D2_Y1[5][MAXBIN],D2_YE1[5][MAXBIN]; 
     Double_t Dp1[5][MAXBIN],Dp_E1[5][MAXBIN]; 

     Double_t H1_x2[5][MAXBIN],H1_Q22[5][MAXBIN],H1_xavg2[5][MAXBIN];
     Double_t D2_x2[5][MAXBIN],D2_Q22[5][MAXBIN],D2_xavg2[5][MAXBIN];
     Double_t H1_Y2[5][MAXBIN],H1_YE2[5][MAXBIN]; 
     Double_t D2_Y2[5][MAXBIN],D2_YE2[5][MAXBIN]; 
     Double_t Dp2[5][MAXBIN],Dp_E2[5][MAXBIN]; 


     for(int ii=0;ii<5;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_xavg[ii][jj]=0.0; 
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_xavg[ii][jj]=0.0;
             H1_Y[ii][jj]=0.0;  H1_YE[ii][jj]=0.0; 
     	     D2_Y[ii][jj]=0.0;  D2_YE[ii][jj]=0.0; 
             Dp[ii][jj]=0.0;   Dp_E[ii][jj]=0.0; 

             H1_x1[ii][jj]=0.0; H1_Q21[ii][jj]=0.0; H1_xavg1[ii][jj]=0.0; 
             D2_x1[ii][jj]=0.0; D2_Q21[ii][jj]=0.0; D2_xavg1[ii][jj]=0.0;
             H1_Y1[ii][jj]=0.0;  H1_YE1[ii][jj]=0.0; 
     	     D2_Y1[ii][jj]=0.0;  D2_YE1[ii][jj]=0.0; 
             Dp1[ii][jj]=0.0;   Dp_E1[ii][jj]=0.0; 

             H1_x2[ii][jj]=0.0; H1_Q22[ii][jj]=0.0; H1_xavg2[ii][jj]=0.0; 
             D2_x2[ii][jj]=0.0; D2_Q22[ii][jj]=0.0; D2_xavg2[ii][jj]=0.0;
             H1_Y2[ii][jj]=0.0;  H1_YE2[ii][jj]=0.0; 
     	     D2_Y2[ii][jj]=0.0;  D2_YE2[ii][jj]=0.0; 
             Dp2[ii][jj]=0.0;   Dp_E2[ii][jj]=0.0; 
     }}

   TString Yfile;
   int kin[5]={0,1,2,3,4};
   for(int ii=0;ii<5;ii++){
       Yfile=Form("Yield_newbin/RawYield/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,D2_x,D2_xavg,D2_Q2,D2_Y,D2_YE); 
       Yfile=Form("Yield_final/RawYield/bin003/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,D2_x1,D2_xavg1,D2_Q21,D2_Y1,D2_YE1); 

       Yfile=Form("Yield_newbin/RawYield/H1_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,H1_x,H1_xavg,H1_Q2,H1_Y,H1_YE);
       Yfile=Form("Yield_final/RawYield/bin003/H1_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,H1_x1,H1_xavg1,H1_Q21,H1_Y1,H1_YE1);
   }

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hratio1=new TGraphErrors();
   TGraphErrors *hKin[5];
   TGraphErrors *hKin1[5];
    
   int nn=0;
   for(int ii=0;ii<5;ii++){
       int nnn=0;
       hKin[ii]=new TGraphErrors();
       for(int jj=0;jj<MAXBIN;jj++){
	   if(D2_x[ii][jj]==0)continue;
           if(abs(D2_xavg[ii][jj]-H1_xavg[ii][jj])>0.001){cout<<"1 Point mismatch !!!"<<endl; continue;}

           if(H1_Y[ii][jj]>0){
             Dp[ii][jj]=D2_Y[ii][jj]/H1_Y[ii][jj];
	     Dp_E[ii][jj]=Dp[ii][jj]*sqrt(pow(D2_YE[ii][jj]/D2_Y[ii][jj],2)+pow(H1_YE[ii][jj]/H1_Y[ii][jj],2));
	   }
           hratio->SetPoint(nn,D2_xavg[ii][jj],Dp[ii][jj]);
           hratio->SetPointError(nn,0.0,Dp_E[ii][jj]);
 
           hKin[ii]->SetPoint(nnn,D2_xavg[ii][jj],Dp[ii][jj]);
           hKin[ii]->SetPointError(nnn,0.0,Dp_E[ii][jj]);
           nn++;
           nnn++;
       }
   } 

   int nn1=0;
   for(int ii=0;ii<5;ii++){
       int nnn=0;
       hKin1[ii]=new TGraphErrors();
       for(int jj=0;jj<MAXBIN;jj++){
	   if(D2_x1[ii][jj]==0)continue;
           if(abs(D2_xavg1[ii][jj]-H1_xavg1[ii][jj])>0.001){cout<<"2 Point mismatch !!!"; continue;}

           if(H1_Y1[ii][jj]>0){
             Dp1[ii][jj]=D2_Y1[ii][jj]/H1_Y1[ii][jj];
	     Dp_E1[ii][jj]=Dp1[ii][jj]*sqrt(pow(D2_YE1[ii][jj]/D2_Y1[ii][jj],2)+pow(H1_YE1[ii][jj]/H1_Y1[ii][jj],2));
	   }
           hratio1->SetPoint(nn1,D2_xavg1[ii][jj],Dp1[ii][jj]);
           hratio1->SetPointError(nn1,0.0,Dp_E1[ii][jj]);

           hKin1[ii]->SetPoint(nnn,D2_xavg1[ii][jj],Dp1[ii][jj]);
           hKin1[ii]->SetPointError(nnn,0.0,Dp_E1[ii][jj]);
           nn1++;
           nnn++;
       }
   } 

   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   TMultiGraph *mg1=new TMultiGraph();
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(1);
   hratio1->SetMarkerStyle(8);
   hratio1->SetMarkerColor(2);
   mg1->Add(hratio);
   mg1->Add(hratio1);
   mg1->Draw("AP");
   mg1->SetTitle("D/p yield ratio;xbj;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio,"newbin","P");
   leg1->AddEntry(hratio1,"bin003","P");
   leg1->Draw();

   TCanvas *c2=new TCanvas("c2","c2",1500,1500);
   TMultiGraph *mg2=new TMultiGraph();
   int color[5]={1,2,3,4,6};
   for(int ii=0;ii<5;ii++){
	hKin[ii]->SetMarkerStyle(8);
	hKin[ii]->SetMarkerColor(color[ii]);
	hKin1[ii]->SetMarkerStyle(22);
	hKin1[ii]->SetMarkerColor(color[ii]);
  	mg2->Add(hKin[ii]);
  	mg2->Add(hKin1[ii]);
   }
   mg2->Draw("AP");
   mg2->SetTitle("D/p yield ratio;xbj;");

   auto leg2=new TLegend(0.7,0.6,0.85,0.85);
   leg2->SetNColumns(2);
   for(int ii=0;ii<5;ii++){
      leg2->AddEntry(hKin[ii],Form("newbin kin%d",ii),"P");
      leg2->AddEntry(hKin1[ii],Form("bin003 kin%d",ii),"P");
   }
   leg2->Draw();
   

}
