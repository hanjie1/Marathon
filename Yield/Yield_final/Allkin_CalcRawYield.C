#include "GetTrees.h"
#include "GetRunList.h"
#include "SetCut.h"
#include "SetCons.h"
#include "CalcLT.h"
#include "CalcCharge.h"
#include "CalcLum.h"
#include "SearchXS.h"
#include <TMath.h>

void Allkin_CalcRawYield(){
   TString target[4]={"H1","D2","He3","H3"};
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

   for(int nn=2;nn<4;nn++){   
    for(int mm=0;mm<11;mm++){
     if(nn==0&&mm>4)break;
     Double_t LUM=CalcLum(kin[mm],target[nn]); //total luminosity get for this kinematics;
     cout<<"Get total Luminosity for target "<<target[nn]<<"  "<<" kin "<<kin[mm]<<" : "<<LUM<<endl;

     ofstream myfile;
     myfile.open(Form("RawYield/bin003/%s_kin%d.txt",target[nn].Data(),kin[mm]));
     myfile<<"n   xbj   Q2   Yield   Yield_err"<<endl;

     vector<Int_t> runList;
     int run_number=0,nrun=0;
     nrun=GetRunList(runList,kin[mm],target[nn]);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

     Double_t xbj[26];
    // Double_t dBin=(xmax[mm]-xmin[mm])/(nBin[mm]*1.0);
     Int_t nBin=26;
     Double_t dBin=0.03;
     for(int ii=0;ii<nBin;ii++){
        // xbj[ii]=xmin[mm]+ii*dBin;
         xbj[ii]=0.15+ii*dBin;
     }

     int KKin=0;  //kinematica series number
     if(kin[mm]<6)KKin=kin[mm];
     if(kin[mm]==7)KKin=6;
     if(kin[mm]==9)KKin=7;
     if(kin[mm]==11)KKin=8;
     if(kin[mm]==13)KKin=9;
     if(kin[mm]==15)KKin=10;
     TCut VZ=Form("L.tr.vz<%f && L.tr.vz>%f",vz_max[KKin],vz_min[KKin]);

     TString TreeName="T";
     TChain* T;
     Double_t totalNe[26]={0.0};
     Double_t RawNe[26]={0.0};
     Double_t totalQ2[26]={0.0};
     Double_t totalXbj[26]={0.0};
     Double_t totalNe_err[26]={0.0};
     for(int ii=0;ii<nrun;ii++){
         run_number=runList[ii];
         TRI_VAR LT=CalcLT(run_number,kin[mm],1);
         Double_t livetime=LT.value; 
         Double_t livetime_err=LT.err; 
         cout<<"Get LT:  "<<livetime<<"  "<<livetime_err<<endl;
         T=GetTree(run_number,kin[mm],TreeName);
         T->Draw(">>electron",trigger2+CK+Ep+beta+ACC+VZ+TRK+W2);
         TEventList *electron;
         gDirectory->GetObject("electron",electron);
         T->SetEventList(electron);

         T->SetBranchStatus("*",0);
         T->SetBranchStatus("L.gold.p",1);
         T->SetBranchStatus("EKLxe.angle",1);
         T->SetBranchStatus("EKLxe.x_bj",1);
         T->SetBranchStatus("EKLxe.Q2",1);

         Double_t aEprime=0.0,aTheta=0.0,axbj=0.0,aQ2=0.0;
         T->SetBranchAddress("L.gold.p",&aEprime);
         T->SetBranchAddress("EKLxe.angle",&aTheta);
         T->SetBranchAddress("EKLxe.x_bj",&axbj);
         T->SetBranchAddress("EKLxe.Q2",&aQ2);

         Double_t Radcor=1.0;
	 Int_t NNe[26]={0};
         Int_t nentries=electron->GetN();
         for(int jj=0;jj<nentries;jj++){
	     T->GetEntry(electron->GetEntry(jj));
             for(int kk=0;kk<nBin;kk++){
		 Double_t dxbj=axbj-xbj[kk];
                 if(dxbj<dBin && dxbj>=0){
		    totalNe[kk]+=1.0/livetime;
		    NNe[kk]++;
                    RawNe[kk]++;
		    totalQ2[kk]+=aQ2;
		    totalXbj[kk]+=axbj;
		    break;
                 }
	     }
         } 

         for(int jj=0;jj<nBin;jj++){
           Double_t tmp_err=0.0;
           if(NNe[jj]!=0){
              tmp_err=sqrt(NNe[jj])/livetime; //LT doesn't have error
           }
	   totalNe_err[jj]+=tmp_err*tmp_err;
         }
         delete T;
         cout<<"Get electron coutns"<<endl;
     }


    Double_t rawYield[26]={0.0};
    Double_t rawYield_err[26]={0.0};
    Double_t avgQ2[26]={0.0};
    Double_t avgXbj[26]={0.0};
    for(int ii=0;ii<nBin;ii++){
        totalNe_err[ii]=sqrt(totalNe_err[ii]);
        //cout<<"Before LUM: "<<ii<<"  "<<totalNe[ii]<<"  "<<LUM<<endl;
	rawYield[ii]=totalNe[ii]/LUM;
	rawYield_err[ii]=totalNe_err[ii]/LUM;
        if(RawNe[ii]!=0){
           avgQ2[ii]=totalQ2[ii]/RawNe[ii];
	   avgXbj[ii]=totalXbj[ii]/RawNe[ii];
        }
	if(avgXbj[ii]==0)continue;
        else{
           myfile<<xbj[ii]<<", "<<avgXbj[ii]<<", "<<avgQ2[ii]<<", "<<rawYield[ii]<<", "<<rawYield_err[ii]<<endl;
        }
    }

   myfile.close();
  }
 }

}
