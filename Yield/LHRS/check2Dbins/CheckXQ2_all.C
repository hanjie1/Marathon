#include "GetRunList.h"
#include "SetCut.h"

void CheckXQ2_all()
{
   TString target[4]={"H1","D2","He3","H3"};
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

   TFile *f1=new TFile("X_Q2_all.root","RECREATE");
   
   for(int ii=1;ii<2;ii++){
    int maxkin;
    if(ii==0)maxkin=5;
    else maxkin=11;
    for(int jj=1;jj<2;jj++){   
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
     TH2F *hX_Q2=new TH2F(hname.Data(),"xbj vs. Q2 for one kin histogram",300,0.16,0.27,900,2.3,4);//1000,0,1,1000,0,14);
     T->Draw(Form("EKLx.Q2:EKLx.x_bj>>%s",hname.Data()),ACC+CK+Ep+trigger2+VZ+beta+TRK,"surf2");
   }
  }

  f1->Write(); 



}
