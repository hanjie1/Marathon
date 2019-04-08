#include "GetTrees.h"
#include "GetRunList.h"
#include "SetCut.h"
#include "SetCons.h"
#include "CalcLT.h"
#include "CalcCharge.h"
#include "CalcLum.h"
#include "SearchXS.h"
#include <TMath.h>

#define nBin 26
#define dBin 0.03

using namespace std;

void Allkin_CalcRawYield(){
   TString target[4]={"H1","D2","He3","H3"};
   int kin[4]={1,4,9,15};

   Double_t Nu_c=7.49;
   Double_t Theta_c[11]={16.8075,17.5717,19.1125,20.575,21.9401,23.2065,25.5858,27.7642,29.8087,31.7274,33.5552};
   Double_t Theta_c1[4]={25.5909,27.7744,29.8159,31.7325};//2nd run of kin 7,9,11,13; 2nd and 3rd kin15 are almost the same as 1st run, so use sam

   for(int nn=0;nn<4;nn++){   
    for(int mm=0;mm<4;mm++){
     if(nn==0&&kin[mm]>4)break;
     Double_t LUM=CalcLum(kin[mm],target[nn]); //total luminosity get for this kinematics;
     cout<<"Get total Luminosity for target "<<target[nn]<<"  "<<" kin "<<kin[mm]<<" : "<<LUM<<endl;

     ofstream myfile;
     myfile.open(Form("RawYield/Nocut/%s_kin%d.txt",target[nn].Data(),kin[mm]));
     myfile<<"n   xbj   Q2   Yield   Yield_err"<<endl;

     vector<vector<Int_t> > runList;
     int run_number=0,nrun=0;
     nrun=GetRunList(runList,kin[mm],target[nn]);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

     Double_t xbj[nBin];
    // Double_t dBin=(xmax[mm]-xmin[mm])/(nBin[mm]*1.0);
     for(int ii=0;ii<nBin;ii++){
        // xbj[ii]=xmin[mm]+ii*dBin;
         xbj[ii]=0.15+ii*dBin;
     }
     cout<<endl;

     TString TreeName="T";
     TChain* T;
     Double_t totalNe[nBin]={0.0};
     Double_t RawNe[nBin]={0.0};
     Double_t totalQ2[nBin]={0.0};
     Double_t totalXbj[nBin]={0.0};
     Double_t totalNe_err[nBin]={0.0};

     int KKin=0;  //kinematica series number
     if(kin[mm]<6)KKin=kin[mm];
     if(kin[mm]==7)KKin=6;
     if(kin[mm]==9)KKin=7;
     if(kin[mm]==11)KKin=8;
     if(kin[mm]==13)KKin=9;
     if(kin[mm]==15)KKin=10;
     TCut VZ=Form("L.tr.vz<%f && L.tr.vz>%f",vz_max[KKin],vz_min[KKin]);

     for(int ii=0;ii<nrun;ii++){
         run_number=runList[ii][0];
         Int_t tmp_runp=runList[ii][1];
         TCut ACC_phy;
         if(kin[mm]<6 || kin[mm]==15)
           ACC_phy=Form("abs(EKLxe.angle*180.0/3.14159-%f)<=1.6",Theta_c[KKin]);
         else{ 
           if(tmp_runp==1)
              ACC_phy=Form("abs(EKLxe.angle*180.0/3.14159-%f)<=1.6",Theta_c[KKin]);
           else 
              ACC_phy=Form("abs(EKLxe.angle*180.0/3.14159-%f)<=1.6",Theta_c1[KKin-6]);
	 }

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
	 Int_t NNe[nBin]={0};
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
             // tmp_err=sqrt(1.0/NNe[jj]+(livetime_err/livetime)*(livetime_err/livetime))*NNe[jj]/livetime;
              tmp_err=sqrt(NNe[jj])/livetime; //LT doesn't have error
           }
	   totalNe_err[jj]+=tmp_err*tmp_err;
         }
         delete T;
         cout<<"Get electron coutns"<<endl;
     }


    Double_t rawYield[nBin]={0.0};
    Double_t rawYield_err[nBin]={0.0};
    Double_t avgQ2[nBin]={0.0};
    Double_t avgXbj[nBin]={0.0};
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
