#include "GetRunList.h"
#include "SetCut.h"

void CheckNuTheta_all()
{
   gStyle->SetOptStat(1111111);
   TString target[4]={"H1","D2","He3","H3"};
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   Double_t Nu_c=7.49;
   Double_t Theta_c[11]={16.8075,17.5717,19.1125,20.575,21.9401,23.2065,25.5858,27.7642,29.8087,31.7274,33.5552};
   Double_t Theta_c1[4]={25.5909,27.7744,29.8159,31.7325};//2nd run of kin 7,9,11,13; 2nd and 3rd kin15 are almost the same as 1st run, so use same center theta value;

   TFile *f1=new TFile("Nu_Theta_all.root","RECREATE");
   
   for(int ii=0;ii<4;ii++){
    int maxkin;
    if(ii==0)maxkin=5;
    else maxkin=11;
    for(int jj=0;jj<maxkin;jj++){   
     int nrun=0;

     TCut VZ = Form("L.tr.vz>%f && L.tr.vz<%f",vz_min[jj],vz_max[jj]);

     vector<Int_t> runList1;
     vector<Int_t> runList2;
     vector<Int_t> runList3;
     nrun=GetRunList(runList1,runList2,runList3,kin[jj],target[ii]);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

     TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass2/";

     if(kin[jj]<7){
        TChain *T=new TChain("T");
        for(int kk=0;kk<runList1.size();kk++){
            TString File=rootpath+Form("kin%d/tritium_%d.root",kin[jj],runList1[kk]);
            if(!gSystem->AccessPathName(File)){
               T->Add(File);
               int index=1;
               File=rootpath+Form("kin%d/tritium_%d_%d.root",kin[jj],runList1[kk],index);
               while(!gSystem->AccessPathName(File)){
                     T->Add(File);
                     index++;
                     File=rootpath+Form("kin%d/tritium_%d_%d.root",kin[jj],runList1[kk],index);
               }
            }  
            else {cout<<runList1[kk]<<" rootfile can't be found"<<endl;}
        }

        //TCanvas *c1=new TCanvas("c1");
        TString hname=Form("%s_kin%d",target[ii].Data(),kin[jj]);
        TH2F *hNu_th=new TH2F(hname.Data(),"Nu vs. theta for one kin histogram",800,-4.0,4.0,400,-0.2,0.2);
        T->Draw(Form("EKLxe.nu-%f:EKLxe.angle*180.0/3.14159-%f>>%s",Nu_c,Theta_c[jj],hname.Data()),ACC+CK+Ep+trigger2+VZ+beta+TRK+W2,"COLZ");
     }

     if(kin[jj]<15 && kin[jj]>6){
        TChain *T1=new TChain("T");
        for(int kk=0;kk<runList1.size();kk++){
            TString File=rootpath+Form("kin%d_1st/tritium_%d.root",kin[jj],runList1[kk]);
            if(!gSystem->AccessPathName(File)){
               T1->Add(File);

               int index=1;
               File=rootpath+Form("kin%d_1st/tritium_%d_%d.root",kin[jj],runList1[kk],index);
               while(!gSystem->AccessPathName(File)){
                     T1->Add(File);
                     index++;
                     File=rootpath+Form("kin%d_1st/tritium_%d_%d.root",kin[jj],runList1[kk],index);
               }
            }  
            else {cout<<runList1[kk]<<" rootfile can't be found"<<endl;}
        }

        //TCanvas *c1=new TCanvas("c1");
        TString hname=Form("%s_kin%d_1st",target[ii].Data(),kin[jj]);
        TH2F *hNu_th=new TH2F(hname.Data(),"Nu vs. theta for one kin histogram",800,-4.0,4.0,400,-0.2,0.2);
        T1->Draw(Form("EKLxe.nu-%f:EKLxe.angle*180.0/3.14159-%f>>%s",Nu_c,Theta_c[jj],hname.Data()),ACC+CK+Ep+trigger2+VZ+beta+TRK+W2,"COLZ");


        TChain *T2=new TChain("T");
        for(int kk=0;kk<runList2.size();kk++){
            TString File=rootpath+Form("kin%d_2nd/tritium_%d.root",kin[jj],runList2[kk]);
            if(!gSystem->AccessPathName(File)){
               T2->Add(File);
               int index=1;
               File=rootpath+Form("kin%d_2nd/tritium_%d_%d.root",kin[jj],runList2[kk],index);
               while(!gSystem->AccessPathName(File)){
                     T2->Add(File);
                     index++;
                     File=rootpath+Form("kin%d_2nd/tritium_%d_%d.root",kin[jj],runList2[kk],index);
               }
            }  
            else {cout<<runList2[kk]<<" rootfile can't be found"<<endl;}
        }

        //TCanvas *c1=new TCanvas("c1");
        TString hname1=Form("%s_kin%d_2nd",target[ii].Data(),kin[jj]);
        TH2F *hNu_th_1=new TH2F(hname1.Data(),"Nu vs. theta for one kin histogram",800,-4.0,4.0,400,-0.2,0.2);
        T2->Draw(Form("EKLxe.nu-%f:EKLxe.angle*180.0/3.14159-%f>>%s",Nu_c,Theta_c1[jj-6],hname1.Data()),ACC+CK+Ep+trigger2+VZ+beta+TRK+W2,"COLZ");
cout<<"theta_c1 "<<Theta_c1[jj-6]<<endl;

     }

     if(kin[jj]==15){
        TChain *T=new TChain("T");
        for(int kk=0;kk<runList1.size();kk++){
            TString File=rootpath+Form("kin%d_1st/tritium_%d.root",kin[jj],runList1[kk]);
            if(!gSystem->AccessPathName(File)){
               T->Add(File);
               int index=1;
               File=rootpath+Form("kin%d_1st/tritium_%d_%d.root",kin[jj],runList1[kk],index);
               while(!gSystem->AccessPathName(File)){
                     T->Add(File);
                     index++;
                     File=rootpath+Form("kin%d_1st/tritium_%d_%d.root",kin[jj],runList1[kk],index);
               }
            }  
            else {cout<<runList1[kk]<<" rootfile can't be found"<<endl;}
        }

        for(int kk=0;kk<runList2.size();kk++){
            TString File=rootpath+Form("kin%d_2nd/tritium_%d.root",kin[jj],runList2[kk]);
            if(!gSystem->AccessPathName(File)){
               T->Add(File);
               int index=1;
               File=rootpath+Form("kin%d_2nd/tritium_%d_%d.root",kin[jj],runList2[kk],index);
               while(!gSystem->AccessPathName(File)){
                     T->Add(File);
                     index++;
                     File=rootpath+Form("kin%d_2nd/tritium_%d_%d.root",kin[jj],runList2[kk],index);
               }
            }  
            else {cout<<runList2[kk]<<" rootfile can't be found"<<endl;}
        }

        for(int kk=0;kk<runList3.size();kk++){
            TString File=rootpath+Form("kin%d_3rd/tritium_%d.root",kin[jj],runList3[kk]);
            if(!gSystem->AccessPathName(File)){
               T->Add(File);
               int index=1;
               File=rootpath+Form("kin%d_3rd/tritium_%d_%d.root",kin[jj],runList3[kk],index);
               while(!gSystem->AccessPathName(File)){
                     T->Add(File);
                     index++;
                     File=rootpath+Form("kin%d_3rd/tritium_%d_%d.root",kin[jj],runList3[kk],index);
               }
            }  
            else {cout<<runList3[kk]<<" rootfile can't be found"<<endl;}
        }

        //TCanvas *c1=new TCanvas("c1");
        TString hname=Form("%s_kin%d",target[ii].Data(),kin[jj]);
        TH2F *hNu_th=new TH2F(hname.Data(),"Nu vs. theta for one kin histogram",800,-4.0,4.0,400,-0.2,0.2);
        T->Draw(Form("EKLxe.nu-%f:EKLxe.angle*180.0/3.14159-%f>>%s",Nu_c,Theta_c[jj],hname.Data()),ACC+CK+Ep+trigger2+VZ+beta+TRK+W2,"COLZ");
     }
   }
  }

  f1->Write(); 



}
