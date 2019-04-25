#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He_kin()
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
       Yfile=Form("H3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x,H3_xavg,H3_Q2,H3_Y,H3_YE); 
       Yfile=Form("H3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x1,H3_xavg1,H3_Q21,H3_Y1,H3_YE1); 

       Yfile=Form("He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],He_x,He_xavg,He_Q2,He_Y,He_YE);
       Yfile=Form("He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,kin[ii],He_x1,He_xavg1,He_Q21,He_Y1,He_YE1);
   }

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hratio1=new TGraphErrors();
   TGraphErrors *hKin[11];
   TGraphErrors *hKin1[11];

   ofstream outfile;
   outfile.open("Ratio_H3He.dat");
   int nn=0;
   for(int ii=0;ii<11;ii++){
       hKin[ii]=new TGraphErrors();
       int nnn=0;
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x[ii][jj]==0)continue;
           if(abs(H3_xavg[ii][jj]-He_xavg[ii][jj])>0.001){cout<<"1 Point mismatch !!!  "<<H3_x[ii][jj]<<endl; continue;}

           if(He_Y[ii][jj]>0){
             H3He[ii][jj]=H3_Y[ii][jj]/He_Y[ii][jj];
	     H3He_E[ii][jj]=H3He[ii][jj]*sqrt(pow(H3_YE[ii][jj]/H3_Y[ii][jj],2)+pow(He_YE[ii][jj]/He_Y[ii][jj],2));
	   }
           hratio->SetPoint(nn,H3_xavg[ii][jj],H3He[ii][jj]);
           hratio->SetPointError(nn,0.0,H3He_E[ii][jj]);

           hKin[ii]->SetPoint(nnn,H3_xavg[ii][jj],H3He[ii][jj]);
           hKin[ii]->SetPointError(nnn,0.0,H3He_E[ii][jj]);

	   outfile<<H3_xavg[ii][jj]<<"  "<<H3He[ii][jj]<<"  "<<H3He_E[ii][jj]<<endl;
           nn++;
           nnn++;
       }
   } 
   outfile.close();

   int nn1=0;
   for(int ii=0;ii<11;ii++){
       hKin1[ii]=new TGraphErrors();
       int nnn=0;
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x1[ii][jj]==0)continue;
           if(abs(H3_xavg1[ii][jj]-He_xavg1[ii][jj])>0.001){cout<<"2 Point mismatch !!!"<<H3_x1[ii][jj]<<endl; continue;}

           if(He_Y1[ii][jj]>0){
             H3He1[ii][jj]=H3_Y1[ii][jj]/He_Y1[ii][jj];
	     H3He_E1[ii][jj]=H3He1[ii][jj]*sqrt(pow(H3_YE1[ii][jj]/H3_Y1[ii][jj],2)+pow(He_YE1[ii][jj]/He_Y1[ii][jj],2));
	   }
           hratio1->SetPoint(nn1,H3_xavg1[ii][jj],H3He1[ii][jj]);
           hratio1->SetPointError(nn1,0.0,H3He_E1[ii][jj]);

           hKin1[ii]->SetPoint(nnn,H3_xavg1[ii][jj],H3He1[ii][jj]);
           hKin1[ii]->SetPointError(nnn,0.0,H3He_E1[ii][jj]);
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
   mg1->SetTitle("H3/D yield ratio;xbj;");

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
   mg2->SetTitle("H3/He yield ratio;xbj;");

   auto leg2=new TLegend(0.7,0.6,0.85,0.85);
   leg2->SetNColumns(2);
   for(int ii=0;ii<11;ii++){
      leg2->AddEntry(hKin[ii],Form("newbin kin%d",kin[ii]),"P");
      leg2->AddEntry(hKin1[ii],Form("bin003 kin%d",kin[ii]),"P");
   }
   leg2->Draw();

}
