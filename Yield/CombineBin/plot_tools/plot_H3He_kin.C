#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He_kin()
{
     Double_t He3_x[12][MAXBIN],He3_Q2[12][MAXBIN],He3_xavg[12][MAXBIN];
     Double_t H3_x[12][MAXBIN],H3_Q2[12][MAXBIN],H3_xavg[12][MAXBIN];
     Double_t He3_Y[12][MAXBIN],He3_YE[12][MAXBIN]; 
     Double_t H3_Y[12][MAXBIN],H3_YE[12][MAXBIN]; 
     Double_t H3He[12][MAXBIN],H3He_E[12][MAXBIN]; 

     Double_t He3_x1[12][MAXBIN],He3_Q21[12][MAXBIN],He3_xavg1[12][MAXBIN];
     Double_t H3_x1[12][MAXBIN],H3_Q21[12][MAXBIN],H3_xavg1[12][MAXBIN];
     Double_t He3_Y1[12][MAXBIN],He3_YE1[12][MAXBIN]; 
     Double_t H3_Y1[12][MAXBIN],H3_YE1[12][MAXBIN]; 
     Double_t H3He1[12][MAXBIN],H3He_E1[12][MAXBIN]; 

     Double_t He3_Radx[12][MAXBIN],He3_RadQ2[12][MAXBIN];
     Double_t H3_Radx[12][MAXBIN],H3_RadQ2[12][MAXBIN];
     Double_t He3_RadCor[12][MAXBIN]; 
     Double_t H3_RadCor[12][MAXBIN]; 
     Double_t H3He_RadCor[12][MAXBIN]; 

     Double_t H3_CouX1[12][MAXBIN],H3_CouX2[12][MAXBIN];
     Double_t He3_CouX1[12][MAXBIN],He3_CouX2[12][MAXBIN];
     Double_t H3_born[12][MAXBIN],H3_bornEff[12][MAXBIN];
     Double_t He3_born[12][MAXBIN],He3_bornEff[12][MAXBIN];
     Double_t H3He_CouCor[12][MAXBIN];


     for(int ii=0;ii<12;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_xavg[ii][jj]=0.0; 
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_xavg[ii][jj]=0.0;
             He3_Y[ii][jj]=0.0;  He3_YE[ii][jj]=0.0; 
     	     H3_Y[ii][jj]=0.0;  H3_YE[ii][jj]=0.0; 
             H3He[ii][jj]=0.0;   H3He_E[ii][jj]=0.0; 

             He3_x1[ii][jj]=0.0; He3_Q21[ii][jj]=0.0; He3_xavg1[ii][jj]=0.0; 
             H3_x1[ii][jj]=0.0; H3_Q21[ii][jj]=0.0; H3_xavg1[ii][jj]=0.0;
             He3_Y1[ii][jj]=0.0;  He3_YE1[ii][jj]=0.0; 
     	     H3_Y1[ii][jj]=0.0;  H3_YE1[ii][jj]=0.0; 
             H3He1[ii][jj]=0.0;   H3He_E1[ii][jj]=0.0; 

             He3_Radx[ii][jj]=0.0; He3_RadQ2[ii][jj]=0.0; 
             H3_Radx[ii][jj]=0.0; H3_RadQ2[ii][jj]=0.0;
             He3_RadCor[ii][jj]=0.0; 
     	     H3_RadCor[ii][jj]=0.0; 
             H3He_RadCor[ii][jj]=0.0;  

             H3_CouX1[ii][jj]=0.0; H3_CouX2[ii][jj]=0.0;
             He3_CouX1[ii][jj]=0.0; He3_CouX2[ii][jj]=0.0;
             H3_born[ii][jj]=0.0; H3_bornEff[ii][jj]=0.0;
             He3_born[ii][jj]=0.0; He3_bornEff[ii][jj]=0.0;
             H3He_CouCor[ii][jj]=0.0;
     }}

   TString Yfile;
   int kin[12]={0,1,2,3,4,5,7,9,11,13,15,16};
   for(int ii=0;ii<12;ii++){
       Yfile=Form("newbin/RawYield/H3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,H3_x,H3_xavg,H3_Q2,H3_Y,H3_YE); 
       Yfile=Form("bin003/RawYield/H3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,H3_x1,H3_xavg1,H3_Q21,H3_Y1,H3_YE1); 

       Yfile=Form("newbin/RawYield/He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,He3_x,He3_xavg,He3_Q2,He3_Y,He3_YE);
       Yfile=Form("bin003/RawYield/He3_kin%d.txt",kin[ii]);
       ReadYield(Yfile,ii,He3_x1,He3_xavg1,He3_Q21,He3_Y1,He3_YE1);

       Yfile=Form("newbin/RadCor/He3_kin%d_xs.out",kin[ii]);
       ReadRadCor(Yfile,ii,He3_Radx,He3_RadQ2,He3_RadCor);

       Yfile=Form("newbin/RadCor/H3_kin%d_xs.out",kin[ii]);
       ReadRadCor(Yfile,ii,H3_Radx,H3_RadQ2,H3_RadCor);

       Yfile=Form("newbin/RadCor/H3_kin%d_xs.out",kin[ii]);
       ReadCoulomb(Yfile,ii,H3_CouX1,H3_born);
       Yfile=Form("newbin/Coulomb_Q2eff/H3_kin%d_xs.out",kin[ii]);
       ReadCoulomb(Yfile,ii,H3_CouX2,H3_bornEff);

       Yfile=Form("newbin/RadCor/He3_kin%d_xs.out",kin[ii]);
       ReadCoulomb(Yfile,ii,He3_CouX1,He3_born);
       Yfile=Form("newbin/Coulomb_Q2eff/He3_kin%d_xs.out",kin[ii]);
       ReadCoulomb(Yfile,ii,He3_CouX2,He3_bornEff);

   }

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hratio1=new TGraphErrors();
   TGraphErrors *hKin[12];
   TGraphErrors *hKin1[12];
   
   ofstream outfile;
   outfile.open("newbin/Ratio_H3He.dat"); 
   ofstream outfile1;
   outfile1.open("newbin/Ratio_H3He_statE.dat");
   ofstream outfile2;
   outfile2.open("newbin/Ratio_H3He_Coulomb.dat");
   int nn=0;
   for(int ii=0;ii<12;ii++){
       int nnn=0;
       hKin[ii]=new TGraphErrors();
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x[ii][jj]==0)continue;
           if(abs(H3_xavg[ii][jj]-He3_xavg[ii][jj])>0.001){cout<<"1 Point mismatch !!!"<<endl; continue;}

           if(He3_Y[ii][jj]>0){
             H3He[ii][jj]=H3_Y[ii][jj]/He3_Y[ii][jj];
	     H3He_E[ii][jj]=H3He[ii][jj]*sqrt(pow(H3_YE[ii][jj]/H3_Y[ii][jj],2)+pow(He3_YE[ii][jj]/He3_Y[ii][jj],2));
	   }
           hratio->SetPoint(nn,H3_xavg[ii][jj],H3He[ii][jj]);
           hratio->SetPointError(nn,0.0,H3He_E[ii][jj]);
 
           hKin[ii]->SetPoint(nnn,H3_xavg[ii][jj],H3He[ii][jj]);
           hKin[ii]->SetPointError(nnn,0.0,H3He_E[ii][jj]);

	   if(abs(H3_xavg[ii][jj]-H3_Radx[ii][jj])<0.0001){
	      H3He_RadCor[ii][jj]=He3_RadCor[ii][jj]/H3_RadCor[ii][jj];
	   }
	   else cout<<"Something wrong with RadCor !!"<<endl;

           Double_t tmpH3_Ccor=0.0,tmpHe3_Ccor=0.0;
           if((abs(He3_CouX1[ii][jj]-He3_CouX2[ii][jj])<0.0001) && (abs(H3_CouX1[ii][jj]-H3_CouX2[ii][jj])<0.0001)){
               tmpH3_Ccor=H3_born[ii][jj]/H3_bornEff[ii][jj];
               tmpHe3_Ccor=He3_born[ii][jj]/He3_bornEff[ii][jj];
              if(abs(He3_CouX1[ii][jj]-H3_CouX2[ii][jj])<0.001){
                  H3He_CouCor[ii][jj]=tmpH3_Ccor/tmpHe3_Ccor;
              }
              else cout<<"Something wrong with Coulomb 1 !!"<<endl;
              cout<<He3_CouX1[ii][jj]<<"  "<<H3_CouX2[ii][jj]<<endl;
           }
           else cout<<"Something wrong with Coulomb 2 !!"<<endl;

	   outfile<<H3_xavg[ii][jj]<<"  "<<H3_Q2[ii][jj]<<"  "<<H3He[ii][jj]<<"  "<<H3He_E[ii][jj]<<"  "<<H3He_RadCor[ii][jj]<<"  "<<H3He_CouCor[ii][jj]<<"  "<<kin[ii]<<endl;
           outfile1<<H3_xavg[ii][jj]<<"  "<<H3_Q2[ii][jj]<<"  "<<H3He[ii][jj]<<"  "<<H3He_E[ii][jj]/H3He[ii][jj]<<"  "<<kin[ii]<<endl;
           outfile2<<H3_xavg[ii][jj]<<"  "<<tmpH3_Ccor<<"  "<<tmpHe3_Ccor<<"  "<<kin[ii]<<endl;

           nn++;
           nnn++;
       }
   } 
   outfile.close();
   outfile1.close();
   outfile2.close();

   int nn1=0;
   for(int ii=0;ii<12;ii++){
       int nnn=0;
       hKin1[ii]=new TGraphErrors();
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x1[ii][jj]==0)continue;
           if(abs(H3_xavg1[ii][jj]-He3_xavg1[ii][jj])>0.001){cout<<"2 Point mismatch !!!"; continue;}

           if(He3_Y1[ii][jj]>0){
             H3He1[ii][jj]=H3_Y1[ii][jj]/He3_Y1[ii][jj];
	     H3He_E1[ii][jj]=H3He1[ii][jj]*sqrt(pow(H3_YE1[ii][jj]/H3_Y1[ii][jj],2)+pow(He3_YE1[ii][jj]/He3_Y1[ii][jj],2));
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
   mg1->SetTitle("H3/He3 yield ratio;xbj;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio,"newbin","P");
   leg1->AddEntry(hratio1,"bin003","P");
   leg1->Draw();

   TCanvas *c2=new TCanvas("c2","c2",1500,1500);
   TMultiGraph *mg2=new TMultiGraph();
   int color[12]={1,2,3,4,6,7,8,9,46,30,12,38};
   for(int ii=0;ii<12;ii++){
	hKin[ii]->SetMarkerStyle(8);
	hKin[ii]->SetMarkerColor(color[ii]);
	hKin1[ii]->SetMarkerStyle(47);
	hKin1[ii]->SetMarkerColor(color[ii]);
	hKin1[ii]->SetMarkerSize(1.5);
  	mg2->Add(hKin[ii]);
  	mg2->Add(hKin1[ii]);
   }
   mg2->Draw("AP");
   mg2->SetTitle("H3/He3 yield ratio;xbj;");

   auto leg2=new TLegend(0.7,0.6,0.85,0.85);
   leg2->SetNColumns(2);
   for(int ii=0;ii<12;ii++){
      leg2->AddEntry(hKin[ii],Form("newbin kin%d",kin[ii]),"P");
      leg2->AddEntry(hKin1[ii],Form("bin003 kin%d",kin[ii]),"P");
   }
   leg2->Draw();
   
   TCanvas *c3=new TCanvas("c3","c3",1500,1500);
   TMultiGraph *mg3=new TMultiGraph();
   for(int ii=0;ii<12;ii++){
  	mg3->Add(hKin[ii]);
   }
   mg3->Draw("AP");
   mg3->SetTitle("H3/He3 yield ratio;Bjorken x;");

   auto leg3=new TLegend(0.7,0.6,0.85,0.85);
   for(int ii=0;ii<12;ii++){
      leg3->AddEntry(hKin[ii],Form("kin%d",kin[ii]),"P");
   }
   leg3->Draw();
/*
   TCanvas *c4=new TCanvas("c4","c4",1500,1500);
   TMultiGraph *mg4=new TMultiGraph();
   for(int ii=0;ii<12;ii++){
  	mg4->Add(hKin1[ii]);
   }
   mg4->Draw("AP");
   mg4->SetTitle("D/p yield ratio;xbj;");

   auto leg4=new TLegend(0.7,0.6,0.812,0.812);
   for(int ii=0;ii<12;ii++){
      leg4->AddEntry(hKin1[ii],Form("bin003 kin%d",ii),"P");
   }
   leg4->Draw();
*/
}
