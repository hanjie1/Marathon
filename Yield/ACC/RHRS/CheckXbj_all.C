#include "GetRunList_All.h"
#include "SetCut.h"

void CheckXbj_all()
{
   gStyle->SetOptStat(1111111);
   TString target[3]={"D2","He3","H3"};
   int kin=16;
   TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass2/";

   TFile *f1=new TFile("Xbj_kin16.root","RECREATE");
   
   for(int ii=0;ii<3;ii++){
     vector<vector<Int_t> > runList;
     int run_number=0,nrun=0;
     nrun=GetRunList_All(runList,kin,target[ii]);

     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

     TChain *T=new TChain("T");
     for(int kk=0;kk<nrun;kk++){
         run_number=runList[kk][0];
         int tmp_list=runList[kk][1];

	 if(tmp_list==0){
            TString filename=rootpath+Form("kin%d/tritium_%d.root",kin,run_number);
	     if(!gSystem->AccessPathName(filename)){
		T->Add(filename);

                int index=1;
                filename=rootpath+Form("kin%d/tritium_%d_%d.root",kin,run_number,index);
                while(!gSystem->AccessPathName(filename)){
                    T->Add(filename);
                    index++;
                    filename=rootpath+Form("kin%d/tritium_%d_%d.root",kin,run_number,index);
                }
	     }
	 }

	 if(tmp_list==1){
            TString filename=rootpath+Form("kin%d_1st/tritium_%d.root",kin,run_number);
	     if(!gSystem->AccessPathName(filename)){
		T->Add(filename);

                int index=1;
                filename=rootpath+Form("kin%d_1st/tritium_%d_%d.root",kin,run_number,index);
                while(!gSystem->AccessPathName(filename)){
                    T->Add(filename);
                    index++;
                    filename=rootpath+Form("kin%d_1st/tritium_%d_%d.root",kin,run_number,index);
                }
	     }
	 }

	 if(tmp_list==2){
            TString filename=rootpath+Form("kin%d_2nd/tritium_%d.root",kin,run_number);
	     if(!gSystem->AccessPathName(filename)){
		T->Add(filename);

                int index=1;
                filename=rootpath+Form("kin%d_2nd/tritium_%d_%d.root",kin,run_number,index);
                while(!gSystem->AccessPathName(filename)){
                    T->Add(filename);
                    index++;
                    filename=rootpath+Form("kin%d_2nd/tritium_%d_%d.root",kin,run_number,index);
                }
	     }
	 }

	 if(tmp_list==3){
            TString filename=rootpath+Form("kin%d_3rd/tritium_%d.root",kin,run_number);
	     if(!gSystem->AccessPathName(filename)){
		T->Add(filename);

                int index=1;
                filename=rootpath+Form("kin%d_3rd/tritium_%d_%d.root",kin,run_number,index);
                while(!gSystem->AccessPathName(filename)){
                    T->Add(filename);
                    index++;
                    filename=rootpath+Form("kin%d_3rd/tritium_%d_%d.root",kin,run_number,index);
                }
	     }
	 }
     }
     TString hname=Form("%s_kin%d",target[ii].Data(),kin);
     TH1F *hXbj=new TH1F(hname.Data(),"Xbj for one kin histogram",750,0.15,0.9);
     T->Draw(Form("EKRxe.x_bj>>%s",hname.Data()),ACC+CK+Ep+trigger2+VZ+beta+TRK+W2);
   }

  f1->Write(); 



}
