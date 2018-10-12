#include "GetRunList.h"
#include "SetCut.h"

void CheckXbj_all()
{
   TString target[4]={"H1","D2","He3","H3"};
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

   TFile *f1=new TFile("Phys_all.root","RECREATE");
   
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
     T->SetBranchStatus("EKLx.Q2",1);
     T->SetBranchStatus("EKLx.angle",1);
     T->SetBranchStatus("EKLx.nu",1);
     T->SetBranchStatus("EKLx.W2",1);
     T->SetBranchStatus("L.gold.dp",1);

     Double_t axbj=0.0,aQ2=0.0,aTheta=0.0,aNu=0.0,aW2=0.0,adp=0.0;
     T->SetBranchAddress("EKLx.x_bj",&axbj);
     T->SetBranchAddress("EKLx.Q2",&aQ2);
     T->SetBranchAddress("EKLx.angle",&aTheta);
     T->SetBranchAddress("EKLx.nu",&aNu);
     T->SetBranchAddress("EKLx.W2",&aW2);
     T->SetBranchAddress("L.gold.dp",&adp);

     Double_t x_min=1.0,x_max=0.0;
     Double_t Q2_min=100.0,Q2_max=0.0;
     Double_t W2_min=100.0,W2_max=0.0;
     Double_t th_min=100.0,th_max=0.0;
     Double_t Nu_min=100.0,Nu_max=0.0;
     Double_t dp_min=100.0,dp_max=-100.0;

     Int_t nentries=electron->GetN();
     for(int jj=0;jj<nentries;jj++){
         T->GetEntry(electron->GetEntry(jj));
         if(axbj<x_min)x_min=axbj;
         if(axbj>x_max)x_max=axbj;
         if(aQ2<Q2_min)Q2_min=aQ2;
         if(aQ2>Q2_max)Q2_max=aQ2;
         if(aW2<W2_min)W2_min=aW2;
         if(aW2>W2_max)W2_max=aW2;
         if(aNu<Nu_min)Nu_min=aNu;
         if(aNu>Nu_max)Nu_max=aNu;
         if(adp<dp_min)dp_min=adp;
         if(adp>dp_max)dp_max=adp;
         Double_t tmp_theta=aTheta*180/TMath::Pi();
         if(tmp_theta<th_min)th_min=tmp_theta;
         if(tmp_theta>th_max)th_max=tmp_theta;
     }
     
     TString hname=Form("%s_kin%d_xbj",target[ii].Data(),kin[jj]);
     int nbin_x=(int)((x_max-x_min)/0.001+1);
     TH1F *hxbj=new TH1F(hname.Data(),"xbj for one kin histogram",nbin_x,x_min,x_min+nbin_x*0.001);
     T->Draw(Form("EKLx.x_bj>>%s",hname.Data()));

     hname=Form("%s_kin%d_Q2",target[ii].Data(),kin[jj]);
     int nbin_Q2=(int)((Q2_max-Q2_min)/0.005+1);
     TH1F *hQ2=new TH1F(hname.Data(),"Q2 for one kin histogram",nbin_Q2,Q2_min,Q2_min+nbin_Q2*0.005);
     T->Draw(Form("EKLx.Q2>>%s",hname.Data()));

     hname=Form("%s_kin%d_W2",target[ii].Data(),kin[jj]);
     int nbin_W2=(int)((W2_max-W2_min)/0.005+1);
     TH1F *hW2=new TH1F(hname.Data(),"W2 for one kin histogram",nbin_W2,W2_min,W2_min+nbin_W2*0.005);
     T->Draw(Form("EKLx.W2>>%s",hname.Data()));

     hname=Form("%s_kin%d_Nu",target[ii].Data(),kin[jj]);
     int nbin_Nu=(int)((Nu_max-Nu_min)/0.001+1);
     TH1F *hNu=new TH1F(hname.Data(),"Nu for one kin histogram",nbin_Nu,Nu_min,Nu_min+nbin_Nu*0.001);
     T->Draw(Form("EKLx.nu>>%s",hname.Data()));

     hname=Form("%s_kin%d_th",target[ii].Data(),kin[jj]);
     int nbin_th=(int)((th_max-th_min)/0.001+1);
     TH1F *hth=new TH1F(hname.Data(),"Theta for one kin histogram",nbin_th,th_min,th_min+nbin_th*0.001);
     T->Draw(Form("EKLx.angle*180/3.14159265358979323>>%s",hname.Data()));

     hname=Form("%s_kin%d_dp",target[ii].Data(),kin[jj]);
     int nbin_dp=(int)((dp_max-dp_min)/0.001+1);
     TH1F *hdp=new TH1F(hname.Data(),"dp for one kin histogram",nbin_dp,dp_min,dp_min+nbin_dp*0.001);
     T->Draw(Form("L.gold.dp>>%s",hname.Data()));

     delete electron;
     delete T;
   }
  }

  f1->Write(); 



}
