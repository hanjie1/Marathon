#include "GetRunList.h"
#include "SetCut.h"

void CheckXbj()
{
     TString target;
     int kin;
     cout<<"Target:  ";
     cin>>target;
     cout<<"Kin:     ";
     cin>>kin;

     int nrun=0;
     vector<Int_t> runList;
     nrun=GetRunList(runList,kin,target);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

     TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass1/";
     TChain *T=new TChain("T");
     for(int ii=0;ii<runList.size();ii++){
        TString File=rootpath+Form("kin%d/tritium_%d.root",kin,runList[ii]);
        if(!gSystem->AccessPathName(File)){
           T->Add(File);
           int index=1;
           File=rootpath+Form("kin%d/tritium_%d_%d.root",kin,runList[ii],index);
           while(!gSystem->AccessPathName(File)){
                  T->Add(File);
                  index++;
                  File=rootpath+Form("kin%d/tritium_%d_%d.root",kin,runList[ii],index);
           }
        }
        else {cout<<runList[ii]<<" rootfile can't be found"<<endl;}
     }

     TCanvas *c1=new TCanvas("c1");
     TH1F *hxbj=new TH1F("hxbj","xbj for one kin histogram",1000,0,1);
     T->Draw("EKLx.x_bj>>hxbj",ACC+CK+Ep+trigger2+VZ+beta+TRK);



}
