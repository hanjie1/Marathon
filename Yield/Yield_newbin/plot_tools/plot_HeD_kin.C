#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_HeD_kin()
{
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_xavg[11][MAXBIN];
     Double_t He_x[11][MAXBIN],He_Q2[11][MAXBIN],He_xavg[11][MAXBIN];
     Double_t D2_Y[11][MAXBIN],D2_YE[11][MAXBIN]; 
     Double_t He_Y[11][MAXBIN],He_YE[11][MAXBIN]; 
     Double_t HeD[11][MAXBIN],HeD_E[11][MAXBIN]; 

     Double_t D2_x1[11][MAXBIN],D2_Q21[11][MAXBIN],D2_xavg1[11][MAXBIN];
     Double_t He_x1[11][MAXBIN],He_Q21[11][MAXBIN],He_xavg1[11][MAXBIN];
     Double_t D2_Y1[11][MAXBIN],D2_YE1[11][MAXBIN]; 
     Double_t He_Y1[11][MAXBIN],He_YE1[11][MAXBIN]; 
     Double_t HeD1[11][MAXBIN],HeD_E1[11][MAXBIN]; 

     Double_t D2_x2[11][MAXBIN],D2_Q22[11][MAXBIN],D2_xavg2[11][MAXBIN];
     Double_t He_x2[11][MAXBIN],He_Q22[11][MAXBIN],He_xavg2[11][MAXBIN];
     Double_t D2_Y2[11][MAXBIN],D2_YE2[11][MAXBIN]; 
     Double_t He_Y2[11][MAXBIN],He_YE2[11][MAXBIN]; 
     Double_t HeD2[11][MAXBIN],HeD_E2[11][MAXBIN]; 

     for(int ii=0;ii<11;ii++){
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

             D2_x2[ii][jj]=0.0; D2_Q22[ii][jj]=0.0; D2_xavg2[ii][jj]=0.0; 
             He_x2[ii][jj]=0.0; He_Q22[ii][jj]=0.0; He_xavg2[ii][jj]=0.0;
             D2_Y2[ii][jj]=0.0;  D2_YE2[ii][jj]=0.0; 
     	     He_Y2[ii][jj]=0.0;  He_YE2[ii][jj]=0.0; 
             HeD2[ii][jj]=0.0;   HeD_E2[ii][jj]=0.0; 
     }}

   TString Yfile;
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   for(int ii=0;ii<11;ii++){
       Yfile=Form("Yield_newbin/RawYield/He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],He_x,He_xavg,He_Q2,He_Y,He_YE); 
       Yfile=Form("Yield_final/RawYield/bin003/He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],He_x1,He_xavg1,He_Q21,He_Y1,He_YE1); 

       Yfile=Form("Yield_newbin/RawYield/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_xavg,D2_Q2,D2_Y,D2_YE);
       Yfile=Form("Yield_final/RawYield/bin003/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x1,D2_xavg1,D2_Q21,D2_Y1,D2_YE1);
   }

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hratio1=new TGraphErrors();
   TGraphErrors *hKin[11];
   TGraphErrors *hKin1[11];

   int nn=0;
   for(int ii=0;ii<11;ii++){
       hKin[ii]=new TGraphErrors();
       int nnn=0;
       for(int jj=0;jj<MAXBIN;jj++){
	   if(He_x[ii][jj]==0)continue;
           if(abs(He_xavg[ii][jj]-D2_xavg[ii][jj])>0.001){cout<<"1 Point mismatch !!!  "<<He_x[ii][jj]<<endl; continue;}

           if(D2_Y[ii][jj]>0){
             HeD[ii][jj]=He_Y[ii][jj]/D2_Y[ii][jj];
	     HeD_E[ii][jj]=HeD[ii][jj]*sqrt(pow(He_YE[ii][jj]/He_Y[ii][jj],2)+pow(D2_YE[ii][jj]/D2_Y[ii][jj],2));
	   }
           hratio->SetPoint(nn,He_xavg[ii][jj],HeD[ii][jj]);
           hratio->SetPointError(nn,0.0,HeD_E[ii][jj]);

           hKin[ii]->SetPoint(nnn,He_xavg[ii][jj],HeD[ii][jj]);
           hKin[ii]->SetPointError(nnn,0.0,HeD_E[ii][jj]);
           nn++;
           nnn++;
       }
   } 

   int nn1=0;
   for(int ii=0;ii<11;ii++){
       hKin1[ii]=new TGraphErrors();
       int nnn=0;
       for(int jj=0;jj<MAXBIN;jj++){
	   if(He_x1[ii][jj]==0)continue;
           if(abs(He_xavg1[ii][jj]-D2_xavg1[ii][jj])>0.001){cout<<"2 Point mismatch !!!"<<He_x1[ii][jj]<<endl; continue;}

           if(D2_Y1[ii][jj]>0){
             HeD1[ii][jj]=He_Y1[ii][jj]/D2_Y1[ii][jj];
	     HeD_E1[ii][jj]=HeD1[ii][jj]*sqrt(pow(He_YE1[ii][jj]/He_Y1[ii][jj],2)+pow(D2_YE1[ii][jj]/D2_Y1[ii][jj],2));
	   }
           hratio1->SetPoint(nn1,He_xavg1[ii][jj],HeD1[ii][jj]);
           hratio1->SetPointError(nn1,0.0,HeD_E1[ii][jj]);

           hKin1[ii]->SetPoint(nnn,He_xavg1[ii][jj],HeD1[ii][jj]);
           hKin1[ii]->SetPointError(nnn,0.0,HeD_E1[ii][jj]);
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
   mg1->SetTitle("He3/D yield ratio;xbj;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio,"newbin","P");
   leg1->AddEntry(hratio1,"bin003","P");
   leg1->Draw();

   TCanvas *c2=new TCanvas("c2","c2",1500,1500);
   TMultiGraph *mg2=new TMultiGraph();
   int color[11]={1,2,3,4,6,7,8,9,46,30,28};
   for(int ii=0;ii<11;ii++){
        hKin[ii]->SetMarkerStyle(8);
        hKin[ii]->SetMarkerColor(color[ii]);
        hKin1[ii]->SetMarkerStyle(22);
        hKin1[ii]->SetMarkerColor(color[ii]);
        mg2->Add(hKin[ii]);
        mg2->Add(hKin1[ii]);
   }
   mg2->Draw("AP");
   mg2->SetTitle("He3/D yield ratio;xbj;");

   auto leg2=new TLegend(0.7,0.6,0.85,0.85);
   leg2->SetNColumns(2);
   for(int ii=0;ii<11;ii++){
      leg2->AddEntry(hKin[ii],Form("newbin kin%d",kin[ii]),"P");
      leg2->AddEntry(hKin1[ii],Form("bin003 kin%d",kin[ii]),"P");
   }
   leg2->Draw();

}
