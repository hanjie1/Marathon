#include "GetTrees.h"
#include "GetRunList.h"
#include "SetCut.h"
#include "SetCons.h"
#include "CalcLT.h"
#include "CalcCharge.h"
#include "CalcLum.h"
#include <TMath.h>

void Allkin_CalcRawYield(){
   TString target[3]={"D2","He3","H3"};
   int kin=16;

   for(int nn=0;nn<3;nn++){   
     Double_t LUM=CalcLum(kin,target[nn]); //total luminosity get for this kinematics;
     cout<<"Get total Luminosity for target "<<target[nn]<<"  "<<" kin "<<kin<<" : "<<LUM<<endl;

     ofstream myfile;
     myfile.open(Form("RawYield/%s_kin%d.txt",target[nn].Data(),kin));
     myfile<<"n   xbj   Q2   Yield   Yield_err"<<endl;

     vector<Int_t> runList;
     int run_number=0,nrun=0;
     nrun=GetRunList(runList,kin,target[nn]);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

     Double_t xbj[3];
     xbj[0]=xmin;
     for(int ii=1;ii<nBin;ii++){
         xbj[ii]=xmin+dBin[ii];
     }

     TString TreeName="T";
     TChain* T;
     Double_t totalNe[3]={0.0};
     Double_t RawNe[3]={0.0};
     Double_t totalQ2[3]={0.0};
     Double_t totalXbj[3]={0.0};
     Double_t totalNe_err[3]={0.0};
     for(int ii=0;ii<nrun;ii++){
         run_number=runList[ii];
         TRI_VAR LT=CalcLT(run_number,kin,1);
         Double_t livetime=LT.value; 
         Double_t livetime_err=LT.err; 
         cout<<"Get LT:  "<<livetime<<"  "<<livetime_err<<endl;

         T=GetTree(run_number,kin,TreeName);
         if(T==NULL)continue;
         T->Draw(">>electron",trigger2+CK+Ep+beta+ACC+VZ+TRK+W2);
         TEventList *electron;
         gDirectory->GetObject("electron",electron);
         T->SetEventList(electron);

         T->SetBranchStatus("*",0);
         T->SetBranchStatus("R.gold.p",1);
         T->SetBranchStatus("EKRxe.angle",1);
         T->SetBranchStatus("EKRxe.x_bj",1);
         T->SetBranchStatus("EKRxe.Q2",1);

         Double_t aEprime=0.0,aTheta=0.0,axbj=0.0,aQ2=0.0;
         T->SetBranchAddress("R.gold.p",&aEprime);
         T->SetBranchAddress("EKRxe.angle",&aTheta);
         T->SetBranchAddress("EKRxe.x_bj",&axbj);
         T->SetBranchAddress("EKRxe.Q2",&aQ2);

	 Int_t NNe[3]={0};
         Int_t nentries=electron->GetN();
         for(int jj=0;jj<nentries;jj++){
	     T->GetEntry(electron->GetEntry(jj));
             for(int kk=0;kk<nBin;kk++){
		 Double_t dxbj=axbj-xbj[kk];
                 if(dxbj<dBin[kk] && dxbj>=0){
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


    Double_t rawYield[3]={0.0};
    Double_t rawYield_err[3]={0.0};
    Double_t avgQ2[3]={0.0};
    Double_t avgXbj[3]={0.0};
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
