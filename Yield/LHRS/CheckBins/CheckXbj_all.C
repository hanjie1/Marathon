#include "GetRunList.h"
#include "SetCut.h"

void CheckXbj_all()
{
   TString target[4]={"H1","D2","He3","H3"};
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

   TFile *f1=new TFile("Xbj_new.root","RECREATE");
   
   for(int ii=0;ii<4;ii++){
    int maxkin;
    if(ii==0)maxkin=5;
    else maxkin=11;
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

     T->Draw(">>electron",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     TEventList *electron;
     gDirectory->GetObject("electron",electron);
     T->SetEventList(electron);

     T->SetBranchStatus("*",0);
     T->SetBranchStatus("EKLx.x_bj",1);

     Double_t axbj=0.0;
     T->SetBranchAddress("EKLx.x_bj",&axbj);

     Double_t x_min=1.0,x_max=0.0;

     Int_t nentries=electron->GetN();
     for(int jj=0;jj<nentries;jj++){
         T->GetEntry(electron->GetEntry(jj));
         if(axbj<x_min)x_min=axbj;
         if(axbj>x_max)x_max=axbj;
     }
     
     TString hname=Form("%s_kin%d",target[ii].Data(),kin[jj]);
     int nbin_x=(int)((x_max-x_min+0.0005)/0.001);
     TH1F *hxbj=new TH1F(hname.Data(),"xbj for one kin histogram",nbin_x,x_min,x_max);
     T->Draw(Form("EKLx.x_bj>>%s",hname.Data()));
     delete electron;
     delete T;
   }
  }

  f1->Write(); 



}
