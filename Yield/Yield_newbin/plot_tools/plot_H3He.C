#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He()
{
     Double_t He_x[11][MAXBIN],He_Q2[11][MAXBIN],He_xavg[11][MAXBIN];
     Double_t H3_x[11][MAXBIN],H3_Q2[11][MAXBIN],H3_xavg[11][MAXBIN];
     Double_t He_Y[11][MAXBIN],He_YE[11][MAXBIN]; 
     Double_t H3_Y[11][MAXBIN],H3_YE[11][MAXBIN]; 
     Double_t H3He[11][MAXBIN],H3He_E[11][MAXBIN]; 

     Double_t He_x1[11][MAXBIN],He_Q21[11][MAXBIN],He_xavg1[11][MAXBIN];
     Double_t H3_x1[11][MAXBIN],H3_Q21[11][MAXBIN],H3_xavg1[11][MAXBIN];
     Double_t He_Y1[11][MAXBIN],He_YE1[11][MAXBIN]; 
     Double_t H3_Y1[11][MAXBIN],H3_YE1[11][MAXBIN]; 
     Double_t H3He1[11][MAXBIN],H3He_E1[11][MAXBIN]; 

     Double_t He_x2[11][MAXBIN],He_Q22[11][MAXBIN],He_xavg2[11][MAXBIN];
     Double_t H3_x2[11][MAXBIN],H3_Q22[11][MAXBIN],H3_xavg2[11][MAXBIN];
     Double_t He_Y2[11][MAXBIN],He_YE2[11][MAXBIN]; 
     Double_t H3_Y2[11][MAXBIN],H3_YE2[11][MAXBIN]; 
     Double_t H3He2[11][MAXBIN],H3He_E2[11][MAXBIN]; 

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             He_x[ii][jj]=0.0; He_Q2[ii][jj]=0.0; He_xavg[ii][jj]=0.0; 
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_xavg[ii][jj]=0.0;
             He_Y[ii][jj]=0.0;  He_YE[ii][jj]=0.0; 
     	     H3_Y[ii][jj]=0.0;  H3_YE[ii][jj]=0.0; 
             H3He[ii][jj]=0.0;   H3He_E[ii][jj]=0.0; 

             He_x1[ii][jj]=0.0; He_Q21[ii][jj]=0.0; He_xavg1[ii][jj]=0.0; 
             H3_x1[ii][jj]=0.0; H3_Q21[ii][jj]=0.0; H3_xavg1[ii][jj]=0.0;
             He_Y1[ii][jj]=0.0;  He_YE1[ii][jj]=0.0; 
     	     H3_Y1[ii][jj]=0.0;  H3_YE1[ii][jj]=0.0; 
             H3He1[ii][jj]=0.0;   H3He_E1[ii][jj]=0.0; 

             He_x2[ii][jj]=0.0; He_Q22[ii][jj]=0.0; He_xavg2[ii][jj]=0.0; 
             H3_x2[ii][jj]=0.0; H3_Q22[ii][jj]=0.0; H3_xavg2[ii][jj]=0.0;
             He_Y2[ii][jj]=0.0;  He_YE2[ii][jj]=0.0; 
     	     H3_Y2[ii][jj]=0.0;  H3_YE2[ii][jj]=0.0; 
             H3He2[ii][jj]=0.0;   H3He_E2[ii][jj]=0.0; 
     }}

   TString Yfile;
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   for(int ii=0;ii<11;ii++){
       Yfile=Form("Yield_newbin/RawYield/H3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x,H3_xavg,H3_Q2,H3_Y,H3_YE); 
       Yfile=Form("Yield_final/RawYield/bin003/H3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x1,H3_xavg1,H3_Q21,H3_Y1,H3_YE1); 
//       Yfile=Form("H3_kin%d.txt",kin[ii]);
//       ReadYield(Yfile,kin[ii],H3_x2,H3_xavg2,H3_Q22,H3_Y2,H3_YE2); 

       Yfile=Form("Yield_newbin/RawYield/He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],He_x,He_xavg,He_Q2,He_Y,He_YE);
       Yfile=Form("Yield_final/RawYield/bin003/He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],He_x1,He_xavg1,He_Q21,He_Y1,He_YE1);
//       Yfile=Form("He_kin%d.txt",kin[ii]);
//       ReadYield(Yfile,kin[ii],He_x2,He_xavg2,He_Q22,He_Y2,He_YE2);
   }

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hratio1=new TGraphErrors();
   TGraphErrors *hratio2=new TGraphErrors();

   ofstream outfile1;
   outfile1.open("Ratio_H3He.dat");
   outfile1<<"newbin: "<<endl;    
   int nn=0;
   for(int ii=0;ii<11;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x[ii][jj]==0)continue;
           if(abs(H3_xavg[ii][jj]-He_xavg[ii][jj])>0.001){cout<<"1 Point mismatch !!!  "<<H3_x[ii][jj]<<endl; continue;}

           if(He_Y[ii][jj]>0){
             H3He[ii][jj]=H3_Y[ii][jj]/He_Y[ii][jj];
	     H3He_E[ii][jj]=H3He[ii][jj]*sqrt(pow(H3_YE[ii][jj]/H3_Y[ii][jj],2)+pow(He_YE[ii][jj]/He_Y[ii][jj],2));
	   }
           hratio->SetPoint(nn,H3_xavg[ii][jj],H3He[ii][jj]);
           hratio->SetPointError(nn,0.0,H3He_E[ii][jj]);
           outfile1<<H3_xavg[ii][jj]<<"  "<<H3He[ii][jj]<<"  "<<H3He_E[ii][jj]/H3He[ii][jj]<<endl;
           nn++;
       }
   } 

   outfile1<<"------------------------------------"<<endl;
   outfile1<<"bin003: "<<endl;
   int nn1=0;
   for(int ii=0;ii<11;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x1[ii][jj]==0)continue;
           if(abs(H3_xavg1[ii][jj]-He_xavg1[ii][jj])>0.001){cout<<"2 Point mismatch !!!"<<H3_x1[ii][jj]<<endl; continue;}

           if(He_Y1[ii][jj]>0){
             H3He1[ii][jj]=H3_Y1[ii][jj]/He_Y1[ii][jj];
	     H3He_E1[ii][jj]=H3He1[ii][jj]*sqrt(pow(H3_YE1[ii][jj]/H3_Y1[ii][jj],2)+pow(He_YE1[ii][jj]/He_Y1[ii][jj],2));
	   }
           hratio1->SetPoint(nn1,H3_xavg1[ii][jj],H3He1[ii][jj]);
           hratio1->SetPointError(nn1,0.0,H3He_E1[ii][jj]);
           outfile1<<H3_xavg1[ii][jj]<<"  "<<H3He1[ii][jj]<<"  "<<H3He_E1[ii][jj]/H3He1[ii][jj]<<endl;
           nn1++;
       }
   }
/*
   outfile1<<"-----------------------------------"<<endl;
   outfile1<<"cut4: "<<endl;
   int nn2=0;
   for(int ii=0;ii<4;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x2[ii][jj]==0)continue;
           if(abs(H3_xavg2[ii][jj]-He_xavg2[ii][jj])>0.001){cout<<"3 point mismatch !!!"<<H3_x2[ii][jj]<<endl; continue;}

           if(He_Y2[ii][jj]>0){
             H3He2[ii][jj]=H3_Y2[ii][jj]/He_Y2[ii][jj];
	     H3He_E2[ii][jj]=H3He2[ii][jj]*sqrt(pow(H3_YE2[ii][jj]/H3_Y2[ii][jj],2)+pow(He_YE2[ii][jj]/He_Y2[ii][jj],2));
	   }
           hratio2->SetPoint(nn2,H3_xavg2[ii][jj],H3He2[ii][jj]);
           hratio2->SetPointError(nn2,0.0,H3He_E2[ii][jj]);
           outfile1<<H3_xavg2[ii][jj]<<"  "<<H3He2[ii][jj]<<"  "<<H3He_E2[ii][jj]/H3He2[ii][jj]<<endl;
           nn2++;
       }
   }
*/   outfile1.close(); 

   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   TMultiGraph *mg1=new TMultiGraph();
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(1);
   hratio1->SetMarkerStyle(8);
   hratio1->SetMarkerColor(2);
   hratio2->SetMarkerStyle(8);
   hratio2->SetMarkerColor(4);
   mg1->Add(hratio);
   mg1->Add(hratio1);
//   mg1->Add(hratio2);
   mg1->Draw("AP");
   mg1->SetTitle("H3/D yield ratio;xbj;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio,"newbin","P");
   leg1->AddEntry(hratio1,"bin003","P");
//   leg1->AddEntry(hratio2,"cut4","P");
   leg1->Draw();

}
