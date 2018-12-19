#include "GetRunList.h"
#include "SetCut.h"

void CheckXQ2_all()
{
   gStyle->SetOptStat(1111111);
   TString target[4]={"H1","D2","He3","H3"};
   int kin[3]={0,4,15};

   TFile *f1=new TFile("X_Q2_all.root","RECREATE");
   
   for(int ii=0;ii<4;ii++){
    int maxkin;
    if(ii==0)maxkin=2;
    else maxkin=3;
    for(int jj=0;jj<maxkin;jj++){   
     int nrun=0;
     vector<Int_t> runList;
     nrun=GetRunList(runList,kin[jj],target[ii]);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

     TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass1/";
     TChain *T=new TChain("T");
     for(int ii=0;ii<runList.size();ii++){
        TString File=rootpath+Form("kin%d/tritium_%d.root",kin[jj],runList[ii]);
        if(!gSystem->AccessPathName(File)){
           T->Add(File);
           int index=1;
           File=rootpath+Form("kin%d/tritium_%d_%d.root",kin[jj],runList[ii],index);
           while(!gSystem->AccessPathName(File)){
                  T->Add(File);
                  index++;
                  File=rootpath+Form("kin%d/tritium_%d_%d.root",kin[jj],runList[ii],index);
           }
        }
        else {cout<<runList[ii]<<" rootfile can't be found"<<endl;}
     }

     //TCanvas *c1=new TCanvas("c1");
     TString hname=Form("%s_kin%d",target[ii].Data(),kin[jj]);
     TH2F *hX_Q2=new TH2F(hname.Data(),"Q2 vs. x for one kin histogram",830,0.12,0.95,2400,1.5,13.5);//1000,0,1,1000,0,14);
     T->Draw(Form("EKLx.Q2:EKLx.x_bj>>%s",hname.Data()),ACC+CK+Ep+trigger2+VZ+beta+TRK,"COLZ");
   }
  }

  f1->Write(); 



}
