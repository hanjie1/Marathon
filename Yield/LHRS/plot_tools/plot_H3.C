#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3()
{
     Double_t H3_x[10][MAXBIN],H3_Q2[10][MAXBIN],H3_Yield[10][MAXBIN],H3_Yerr[10][MAXBIN];
     Double_t H3R_x[10][MAXBIN],H3R_Q2[10][MAXBIN],H3_Born1[10][MAXBIN],H3_Rad1[10][MAXBIN];
     Double_t H3_Born2[10][MAXBIN],H3_Rad2[10][MAXBIN];

     for(int ii=0;ii<10;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Yield[ii][jj]=0.0; H3_Yerr[ii][jj]=0.0;
             H3R_x[ii][jj]=0.0; H3R_Q2[ii][jj]=0.0; H3_Born1[ii][jj]=0.0; H3_Rad1[ii][jj]=0.0;
             H3_Born2[ii][jj]=0.0; H3_Rad2[ii][jj]=0.0;
     }}

   TString Yfile;
   TString RCfile;
   Int_t KIN=1;
   while(KIN<16){
       Yfile=Form("H3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H3_x,H3_Q2,H3_Yield,H3_Yerr); 
       RCfile=Form("gsmearing/H3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,H3R_x,H3R_Q2,H3_Born1,H3_Rad1);  
       RCfile=Form("Bodek/H3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,H3R_x,H3R_Q2,H3_Born2,H3_Rad2);  
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gH3Raw[10];
  for(int ii=0;ii<10;ii++){
      gH3Raw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(H3_x[ii][jj]==0)continue;
          gH3Raw[ii]->SetPoint(nn,H3_x[ii][jj],H3_Yield[ii][jj]);
          gH3Raw[ii]->SetPointError(nn,0.0,H3_Yerr[ii][jj]);
          nn++;
      }
  }
  TGraph *gH3Rad1[10];
  TGraph *gH3Born1[10];
  TGraph *gH3RC1[10];
  TGraph *gH3Rad2[10];
  TGraph *gH3Born2[10];
  TGraph *gH3RC2[10];
  TGraph *gRatio_1[10];
  TGraph *gRatio_2[10];
  for(int ii=0;ii<10;ii++){
      gH3Rad1[ii]=new TGraph();
      gH3Born1[ii]=new TGraph();
      gH3RC1[ii]=new TGraph();
      gH3Rad2[ii]=new TGraph();
      gH3Born2[ii]=new TGraph();
      gH3RC2[ii]=new TGraph();
      gRatio_1[ii]=new TGraph();
      gRatio_2[ii]=new TGraph();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(H3R_x[ii][jj]==0)continue;
          gH3Rad1[ii]->SetPoint(nn,H3R_x[ii][jj],H3_Rad1[ii][jj]);
          gH3Born1[ii]->SetPoint(nn,H3R_x[ii][jj],H3_Born1[ii][jj]);
          gH3RC1[ii]->SetPoint(nn,H3R_x[ii][jj],H3_Born1[ii][jj]/H3_Rad1[ii][jj]);
          gH3Rad2[ii]->SetPoint(nn,H3R_x[ii][jj],H3_Rad2[ii][jj]);
          gH3Born2[ii]->SetPoint(nn,H3R_x[ii][jj],H3_Born2[ii][jj]);
          gH3RC2[ii]->SetPoint(nn,H3R_x[ii][jj],H3_Born2[ii][jj]/H3_Rad2[ii][jj]);
          gRatio_1[ii]->SetPoint(nn,H3R_x[ii][jj],(H3_Born1[ii][jj]/H3_Rad1[ii][jj])/(H3_Born2[ii][jj]/H3_Rad2[ii][jj]));
          gRatio_2[ii]->SetPoint(nn,H3R_x[ii][jj],H3_Born1[ii][jj]/H3_Born2[ii][jj]);
          nn++;
      }
  }
    

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gH3Raw[ii]->SetMarkerStyle(8);
      if(ii==9)gH3Raw[ii]->SetMarkerColor(30);
      else gH3Raw[ii]->SetMarkerColor(ii+1);
      mg1->Add(gH3Raw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("H3 Data Yield;x;Yield");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=2;
  for(int ii=0;ii<10;ii++){
      if(ii<5)leg1->AddEntry(gH3Raw[ii],Form("H3 kin%d",ii+1),"P");
      else {leg1->AddEntry(gH3Raw[ii],Form("H3 kin%d",ii+nn),"P");
	    nn++;	    
      }
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gH3RC1[ii]->SetMarkerStyle(22);
      gH3RC1[ii]->SetMarkerSize(1.5);
      gH3RC1[ii]->SetMarkerColor(4);
      mg2->Add(gH3RC1[ii]);
      gH3RC2[ii]->SetMarkerStyle(22);
      gH3RC2[ii]->SetMarkerSize(1.5);
      gH3RC2[ii]->SetMarkerColor(2);
      mg2->Add(gH3RC2[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("H3 RC factor;x;Born/Rad");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  leg2->AddEntry(gH3RC1[0],"gsmearing","P");
  leg2->AddEntry(gH3RC2[0],"INEFT","P");
  leg2->Draw();

  TCanvas *c3=new TCanvas("c3");
  TMultiGraph *mg3=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gH3Born1[ii]->SetMarkerStyle(22);
      gH3Born1[ii]->SetMarkerSize(1.5);
      gH3Born1[ii]->SetMarkerColor(4);
      mg3->Add(gH3Born1[ii]);
      gH3Born2[ii]->SetMarkerStyle(22);
      gH3Born2[ii]->SetMarkerSize(1.5);
      gH3Born2[ii]->SetMarkerColor(2);
      mg3->Add(gH3Born2[ii]);
  }
  mg3->Draw("AP");
  mg3->SetTitle("H3 Born cross section;x;Born XS");

  auto leg3=new TLegend(0.7,0.6,0.85,0.85);
  leg3->AddEntry(gH3Born1[0],"gsmearing","P");
  leg3->AddEntry(gH3Born2[0],"INEFT","P");
  leg3->Draw();

  TCanvas *c4=new TCanvas("c4");
  TMultiGraph *mg4=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gH3Rad1[ii]->SetMarkerStyle(22);
      gH3Rad1[ii]->SetMarkerSize(1.5);
      gH3Rad1[ii]->SetMarkerColor(4);
      mg4->Add(gH3Rad1[ii]);
      gH3Rad2[ii]->SetMarkerStyle(22);
      gH3Rad2[ii]->SetMarkerSize(1.5);
      gH3Rad2[ii]->SetMarkerColor(2);
      mg4->Add(gH3Rad2[ii]);
  }
  mg4->Draw("AP");
  mg4->SetTitle("H3 Radiative cross section;x;Rad XS");

  auto leg4=new TLegend(0.7,0.6,0.85,0.85);
  leg4->AddEntry(gH3Rad1[0],"gsmearing","P");
  leg4->AddEntry(gH3Rad2[0],"INEFT","P");
  leg4->Draw();

  TCanvas *c5=new TCanvas("c5");
  TMultiGraph *mg5=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gRatio_1[ii]->SetMarkerStyle(22);
      gRatio_1[ii]->SetMarkerSize(1.5);
      gRatio_1[ii]->SetMarkerColor(4);
      mg5->Add(gRatio_1[ii]);
  }
  mg5->Draw("AP");
  mg5->SetTitle("H3 RC factor ratio;x;");

  auto leg5=new TLegend(0.65,0.7,0.85,0.85);
  leg5->AddEntry(gRatio_1[0],"gsmearing/INEFT","P");
  leg5->Draw();

  TCanvas *c6=new TCanvas("c6");
  TMultiGraph *mg6=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gRatio_2[ii]->SetMarkerStyle(22);
      gRatio_2[ii]->SetMarkerSize(1.5);
      gRatio_2[ii]->SetMarkerColor(4);
      mg6->Add(gRatio_2[ii]);
  }
  mg6->Draw("AP");
  mg6->SetTitle("H3 born cross section ratio;x;");

  auto leg6=new TLegend(0.65,0.7,0.85,0.85);
  leg6->AddEntry(gRatio_2[0],"gsmearing/INEFT","P");
  leg6->Draw();



/*
  for(int ii=0;ii<10;ii++){
     cout<<"-------  "<<ii<<"  -------"<<endl;
     for(int jj=0;jj<MAXBIN;jj++){
         if(H3_x[ii][jj]==0)continue;
         //cout<<H3R_x[ii][jj]<<", "<<H3R_Q2[ii][jj]<<", "<<H3_Born[ii][jj]<<", "<<H3_Rad[ii][jj]<<endl;
         cout<<H3_x[ii][jj]<<", "<<H3_Q2[ii][jj]<<", "<<H3_Yield[ii][jj]<<", "<<H3_Yerr[ii][jj]<<endl;
     }
  }
*/


}
