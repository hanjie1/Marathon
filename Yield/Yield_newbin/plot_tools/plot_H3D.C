#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3D()
{
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_xavg[11][MAXBIN];
     Double_t H3_x[11][MAXBIN],H3_Q2[11][MAXBIN],H3_xavg[11][MAXBIN];
     Double_t D2_Y[11][MAXBIN],D2_YE[11][MAXBIN]; 
     Double_t H3_Y[11][MAXBIN],H3_YE[11][MAXBIN]; 
     Double_t H3D[11][MAXBIN],H3D_E[11][MAXBIN]; 

     Double_t D2_x1[11][MAXBIN],D2_Q21[11][MAXBIN],D2_xavg1[11][MAXBIN];
     Double_t H3_x1[11][MAXBIN],H3_Q21[11][MAXBIN],H3_xavg1[11][MAXBIN];
     Double_t D2_Y1[11][MAXBIN],D2_YE1[11][MAXBIN]; 
     Double_t H3_Y1[11][MAXBIN],H3_YE1[11][MAXBIN]; 
     Double_t H3D1[11][MAXBIN],H3D_E1[11][MAXBIN]; 

     Double_t D2_x2[11][MAXBIN],D2_Q22[11][MAXBIN],D2_xavg2[11][MAXBIN];
     Double_t H3_x2[11][MAXBIN],H3_Q22[11][MAXBIN],H3_xavg2[11][MAXBIN];
     Double_t D2_Y2[11][MAXBIN],D2_YE2[11][MAXBIN]; 
     Double_t H3_Y2[11][MAXBIN],H3_YE2[11][MAXBIN]; 
     Double_t H3D2[11][MAXBIN],H3D_E2[11][MAXBIN]; 

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_xavg[ii][jj]=0.0; 
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_xavg[ii][jj]=0.0;
             D2_Y[ii][jj]=0.0;  D2_YE[ii][jj]=0.0; 
     	     H3_Y[ii][jj]=0.0;  H3_YE[ii][jj]=0.0; 
             H3D[ii][jj]=0.0;   H3D_E[ii][jj]=0.0; 

             D2_x1[ii][jj]=0.0; D2_Q21[ii][jj]=0.0; D2_xavg1[ii][jj]=0.0; 
             H3_x1[ii][jj]=0.0; H3_Q21[ii][jj]=0.0; H3_xavg1[ii][jj]=0.0;
             D2_Y1[ii][jj]=0.0;  D2_YE1[ii][jj]=0.0; 
     	     H3_Y1[ii][jj]=0.0;  H3_YE1[ii][jj]=0.0; 
             H3D1[ii][jj]=0.0;   H3D_E1[ii][jj]=0.0; 

             D2_x2[ii][jj]=0.0; D2_Q22[ii][jj]=0.0; D2_xavg2[ii][jj]=0.0; 
             H3_x2[ii][jj]=0.0; H3_Q22[ii][jj]=0.0; H3_xavg2[ii][jj]=0.0;
             D2_Y2[ii][jj]=0.0;  D2_YE2[ii][jj]=0.0; 
     	     H3_Y2[ii][jj]=0.0;  H3_YE2[ii][jj]=0.0; 
             H3D2[ii][jj]=0.0;   H3D_E2[ii][jj]=0.0; 
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

       Yfile=Form("Yield_newbin/RawYield/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_xavg,D2_Q2,D2_Y,D2_YE);
       Yfile=Form("Yield_final/RawYield/bin003/D2_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x1,D2_xavg1,D2_Q21,D2_Y1,D2_YE1);
//       Yfile=Form("D2_kin%d.txt",kin[ii]);
//       ReadYield(Yfile,kin[ii],D2_x2,D2_xavg2,D2_Q22,D2_Y2,D2_YE2);
   }

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hratio1=new TGraphErrors();
   TGraphErrors *hratio2=new TGraphErrors();

   ofstream outfile1;
   outfile1.open("Ratio_H3D.dat");
   outfile1<<"newbin: "<<endl;    
   int nn=0;
   for(int ii=0;ii<11;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x[ii][jj]==0)continue;
           if(abs(H3_xavg[ii][jj]-D2_xavg[ii][jj])>0.001){cout<<"1 Point mismatch !!!  "<<H3_x[ii][jj]<<endl; continue;}

           if(D2_Y[ii][jj]>0){
             H3D[ii][jj]=H3_Y[ii][jj]/D2_Y[ii][jj];
	     H3D_E[ii][jj]=H3D[ii][jj]*sqrt(pow(H3_YE[ii][jj]/H3_Y[ii][jj],2)+pow(D2_YE[ii][jj]/D2_Y[ii][jj],2));
	   }
           hratio->SetPoint(nn,H3_xavg[ii][jj],H3D[ii][jj]);
           hratio->SetPointError(nn,0.0,H3D_E[ii][jj]);
           outfile1<<H3_xavg[ii][jj]<<"  "<<H3D[ii][jj]<<"  "<<H3D_E[ii][jj]/H3D[ii][jj]<<endl;
           nn++;
       }
   } 

   outfile1<<"------------------------------------"<<endl;
   outfile1<<"bin003: "<<endl;
   int nn1=0;
   for(int ii=0;ii<11;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x1[ii][jj]==0)continue;
           if(abs(H3_xavg1[ii][jj]-D2_xavg1[ii][jj])>0.001){cout<<"2 Point mismatch !!!"<<H3_x1[ii][jj]<<endl; continue;}

           if(D2_Y1[ii][jj]>0){
             H3D1[ii][jj]=H3_Y1[ii][jj]/D2_Y1[ii][jj];
	     H3D_E1[ii][jj]=H3D1[ii][jj]*sqrt(pow(H3_YE1[ii][jj]/H3_Y1[ii][jj],2)+pow(D2_YE1[ii][jj]/D2_Y1[ii][jj],2));
	   }
           hratio1->SetPoint(nn1,H3_xavg1[ii][jj],H3D1[ii][jj]);
           hratio1->SetPointError(nn1,0.0,H3D_E1[ii][jj]);
           outfile1<<H3_xavg1[ii][jj]<<"  "<<H3D1[ii][jj]<<"  "<<H3D_E1[ii][jj]/H3D1[ii][jj]<<endl;
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
           if(abs(H3_xavg2[ii][jj]-D2_xavg2[ii][jj])>0.001){cout<<"3 point mismatch !!!"<<H3_x2[ii][jj]<<endl; continue;}

           if(D2_Y2[ii][jj]>0){
             H3D2[ii][jj]=H3_Y2[ii][jj]/D2_Y2[ii][jj];
	     H3D_E2[ii][jj]=H3D2[ii][jj]*sqrt(pow(H3_YE2[ii][jj]/H3_Y2[ii][jj],2)+pow(D2_YE2[ii][jj]/D2_Y2[ii][jj],2));
	   }
           hratio2->SetPoint(nn2,H3_xavg2[ii][jj],H3D2[ii][jj]);
           hratio2->SetPointError(nn2,0.0,H3D_E2[ii][jj]);
           outfile1<<H3_xavg2[ii][jj]<<"  "<<H3D2[ii][jj]<<"  "<<H3D_E2[ii][jj]/H3D2[ii][jj]<<endl;
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
