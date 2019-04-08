#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp()
{
     Double_t H1_x[2][MAXBIN],H1_Q2[2][MAXBIN],H1_xavg[2][MAXBIN];
     Double_t D2_x[2][MAXBIN],D2_Q2[2][MAXBIN],D2_xavg[2][MAXBIN];
     Double_t H1_Y[2][MAXBIN],H1_YE[2][MAXBIN]; 
     Double_t D2_Y[2][MAXBIN],D2_YE[2][MAXBIN]; 
     Double_t Dp[2][MAXBIN],Dp_E[2][MAXBIN]; 

     Double_t H1_x1[2][MAXBIN],H1_Q21[2][MAXBIN],H1_xavg1[2][MAXBIN];
     Double_t D2_x1[2][MAXBIN],D2_Q21[2][MAXBIN],D2_xavg1[2][MAXBIN];
     Double_t H1_Y1[2][MAXBIN],H1_YE1[2][MAXBIN]; 
     Double_t D2_Y1[2][MAXBIN],D2_YE1[2][MAXBIN]; 
     Double_t Dp1[2][MAXBIN],Dp_E1[2][MAXBIN]; 

     for(int ii=0;ii<2;ii++){
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
     }}

   TString Yfile;
   int kin[2]={1,4};
   for(int ii=0;ii<2;ii++){
       Yfile=Form("Nocut/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,D2_x,D2_xavg,D2_Q2,D2_Y,D2_YE); 
       Yfile=Form("cut1/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,D2_x1,D2_xavg1,D2_Q21,D2_Y1,D2_YE1); 
       Yfile=Form("Nocut/H1_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,H1_x,H1_xavg,H1_Q2,H1_Y,H1_YE);
       Yfile=Form("cut1/H1_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,H1_x1,H1_xavg1,H1_Q21,H1_Y1,H1_YE1);
   }

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hratio1=new TGraphErrors();
    
   int nn=0;
   for(int ii=0;ii<2;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(D2_x[ii][jj]==0)continue;
           cout<<D2_xavg[ii][jj]<<"  "<<H1_xavg[ii][jj]<<endl;
           if(abs(D2_xavg[ii][jj]-H1_xavg[ii][jj])>0.001){cout<<"1 Point mismatch !!!"<<endl; continue;}

           if(H1_Y[ii][jj]>0){
             Dp[ii][jj]=D2_Y[ii][jj]/H1_Y[ii][jj];
	     Dp_E[ii][jj]=Dp[ii][jj]*sqrt(pow(D2_YE[ii][jj]/D2_Y[ii][jj],2)+pow(H1_YE[ii][jj]/H1_Y[ii][jj],2));
	   }
           hratio->SetPoint(nn,D2_xavg[ii][jj],Dp[ii][jj]);
           hratio->SetPointError(nn,0.0,Dp_E[ii][jj]);
	
	   if(D2_x1[ii][jj]==0){cout<<"Something wrong here!"<<endl;nn++;continue;}
           if(abs(D2_xavg1[ii][jj]-H1_xavg1[ii][jj])>0.001){cout<<"2 Point mismatch !!!"; nn++;continue;}

           if(H1_Y1[ii][jj]>0){
             Dp1[ii][jj]=D2_Y1[ii][jj]/H1_Y1[ii][jj];
	     Dp_E1[ii][jj]=Dp1[ii][jj]*sqrt(pow(D2_YE1[ii][jj]/D2_Y1[ii][jj],2)+pow(H1_YE1[ii][jj]/H1_Y1[ii][jj],2));
	   }
           hratio1->SetPoint(nn,D2_xavg1[ii][jj],Dp1[ii][jj]);
           hratio1->SetPointError(nn,0.0,Dp_E1[ii][jj]);

           nn++;
       }
   } 

   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   TMultiGraph *mg1=new TMultiGraph();
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(2);
   hratio1->SetMarkerStyle(8);
   hratio1->SetMarkerColor(4);
   mg1->Add(hratio);
   mg1->Add(hratio1);
   mg1->Draw("AP");
   mg1->SetTitle("D/p yield ratio;xbj;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio,"Nocut","P");
   leg1->AddEntry(hratio1,"cut1","P");
   leg1->Draw();

}
