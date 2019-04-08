#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_HeD()
{
     Double_t D2_x[4][MAXBIN],D2_Q2[4][MAXBIN],D2_xavg[4][MAXBIN];
     Double_t He_x[4][MAXBIN],He_Q2[4][MAXBIN],He_xavg[4][MAXBIN];
     Double_t D2_Y[4][MAXBIN],D2_YE[4][MAXBIN]; 
     Double_t He_Y[4][MAXBIN],He_YE[4][MAXBIN]; 
     Double_t HeD[4][MAXBIN],HeD_E[4][MAXBIN]; 

     Double_t D2_x1[4][MAXBIN],D2_Q21[4][MAXBIN],D2_xavg1[4][MAXBIN];
     Double_t He_x1[4][MAXBIN],He_Q21[4][MAXBIN],He_xavg1[4][MAXBIN];
     Double_t D2_Y1[4][MAXBIN],D2_YE1[4][MAXBIN]; 
     Double_t He_Y1[4][MAXBIN],He_YE1[4][MAXBIN]; 
     Double_t HeD1[4][MAXBIN],HeD_E1[4][MAXBIN]; 

     for(int ii=0;ii<4;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_xavg[ii][jj]=0.0; 
             He_x[ii][jj]=0.0; He_Q2[ii][jj]=0.0; He_xavg[ii][jj]=0.0;
             D2_Y[ii][jj]=0.0;  D2_YE[ii][jj]=0.0; 
     	     He_Y[ii][jj]=0.0;  He_YE[ii][jj]=0.0; 
             HeD[ii][jj]=0.0;   HeD_E[ii][jj]=0.0; 

             D2_x1[ii][jj]=0.0; D2_Q21[ii][jj]=0.0; D2_xavg1[ii][jj]=0.0; 
             He_x1[ii][jj]=0.0; He_Q21[ii][jj]=0.0; He_xavg1[ii][jj]=0.0;
             D2_Y1[ii][jj]=0.0;  D2_YE1[ii][jj]=0.0; 
     	     He_Y1[ii][jj]=0.0;  He_YE1[ii][jj]=0.0; 
             HeD1[ii][jj]=0.0;   HeD_E1[ii][jj]=0.0; 
     }}

   TString Yfile;
   int kin[4]={1,4,9,15};
   for(int ii=0;ii<4;ii++){
       Yfile=Form("Nocut/He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,He_x,He_xavg,He_Q2,He_Y,He_YE); 
       Yfile=Form("cut1/He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,He_x1,He_xavg1,He_Q21,He_Y1,He_YE1); 
       Yfile=Form("Nocut/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,D2_x,D2_xavg,D2_Q2,D2_Y,D2_YE);
       Yfile=Form("cut1/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,D2_x1,D2_xavg1,D2_Q21,D2_Y1,D2_YE1);
   }

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hratio1=new TGraphErrors();
    
   int nn=0,nn1=0;
   for(int ii=0;ii<4;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(He_x[ii][jj]==0)continue;
           if(abs(He_xavg[ii][jj]-D2_xavg[ii][jj])>0.001){cout<<"1 Point mismatch !!!  "<<He_x[ii][jj]<<endl; continue;}

           if(D2_Y[ii][jj]>0){
             HeD[ii][jj]=He_Y[ii][jj]/D2_Y[ii][jj];
	     HeD_E[ii][jj]=HeD[ii][jj]*sqrt(pow(He_YE[ii][jj]/He_Y[ii][jj],2)+pow(D2_YE[ii][jj]/D2_Y[ii][jj],2));
	   }
           hratio->SetPoint(nn,He_xavg[ii][jj],HeD[ii][jj]);
           hratio->SetPointError(nn,0.0,HeD_E[ii][jj]);
           nn++;
	
	   if(He_x1[ii][jj]==0)continue;
           if(abs(He_xavg1[ii][jj]-D2_xavg1[ii][jj])>0.001){cout<<"2 Point mismatch !!!"<<He_x1[ii][jj]<<endl; continue;}

           if(D2_Y1[ii][jj]>0){
             HeD1[ii][jj]=He_Y1[ii][jj]/D2_Y1[ii][jj];
	     HeD_E1[ii][jj]=HeD1[ii][jj]*sqrt(pow(He_YE1[ii][jj]/He_Y1[ii][jj],2)+pow(D2_YE1[ii][jj]/D2_Y1[ii][jj],2));
	   }
           hratio1->SetPoint(nn1,He_xavg1[ii][jj],HeD1[ii][jj]);
           hratio1->SetPointError(nn1,0.0,HeD_E1[ii][jj]);

           nn1++;
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
   mg1->SetTitle("He3/D yield ratio;xbj;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio,"Nocut","P");
   leg1->AddEntry(hratio1,"cut1","P");
   leg1->Draw();

}
