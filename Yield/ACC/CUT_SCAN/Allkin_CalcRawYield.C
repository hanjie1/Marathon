#include "GetTrees.h"
#include "GetRunList.h"
#include "SetCut.h"
#include "SetCons.h"
#include "CalcLT.h"
#include "CalcCharge.h"
#include "CalcLum.h"
#include "SearchACC.h"
#include <TMath.h>

using namespace std;

void Allkin_CalcRawYield(){
   TString target[4]={"H1","D2","He3","H3"};
   int kin[4]={1,4,9,15};

   Double_t Nu_c=7.49;
   Double_t Theta_c[11]={16.8075,17.5717,19.1125,20.575,21.9401,23.2065,25.5858,27.7642,29.8087,31.7274,33.5552};
   Double_t Theta_c1[4]={25.5909,27.7744,29.8159,31.7325};//2nd run of kin 7,9,11,13; 2nd and 3rd kin15 are almost the same as 1st run, so use sam

/* load ACC_matrix table */
   ifstream infile;
   infile.open("ACC_matrix4.dat");
   Double_t dTheta[nTh]={-100.0},ACC_table[nTh][nEp]={-100.0};
   int ACC_matrix[nTh][nEp]={0};
   Ssiz_t from=0;
   TString content,tmp;
   int xx=0,yy=0;
   Double_t tmpTh_save=-100.0;
   while(tmp.ReadLine(infile)){
         tmp.Tokenize(content,from,"  ");
         Double_t tmpTh=atof(content.Data());
         tmp.Tokenize(content,from,"  ");
         Double_t tmpEp=atof(content.Data());
         tmp.Tokenize(content,from,"  ");
         Double_t tmpACC=atoi(content.Data());

	 if(xx==0&&yy==0){
            tmpTh_save=tmpTh;
	    dTheta[xx]=tmpTh;
	 }

	 if(tmpTh!=tmpTh_save){
	    xx++;
	    yy=0;
	    dTheta[xx]=tmpTh;
	    ACC_table[xx][yy]=tmpEp;
	    ACC_matrix[xx][yy]=tmpACC;
            tmpTh_save=tmpTh;
	 }
	 else{
	    ACC_table[xx][yy]=tmpEp;
	    ACC_matrix[xx][yy]=tmpACC;
	 }
         yy++;
         from=0;
   }
   infile.close();

   for(int nn=0;nn<4;nn++){   
    for(int mm=0;mm<4;mm++){
     if(nn==0&&kin[mm]>4)break;
     Double_t LUM=CalcLum(kin[mm],target[nn]); //total luminosity get for this kinematics;
     cout<<"Get total Luminosity for target "<<target[nn]<<"  "<<" kin "<<kin[mm]<<" : "<<LUM<<endl;

     ofstream myfile;
     myfile.open(Form("RawYield/newbin/%s_kin%d.txt",target[nn].Data(),kin[mm]));
     myfile<<"n   xbj   Q2   Yield   Yield_err"<<endl;

     vector<vector<Int_t> > runList;
     int run_number=0,nrun=0;
     nrun=GetRunList(runList,kin[mm],target[nn]);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);
/*
     Double_t xbj[nBin];
     for(int ii=0;ii<nBin;ii++){
         xbj[ii]=0.15+ii*dBin;
     }
*/
     Double_t xbj[9]={0.0};
//     Double_t dBin=(xmax[mm]-xmin[mm])/(1.0*nBin[mm]);
     xbj[0]=xmin[mm];
     for(int ii=1;ii<nBin[mm];ii++){
	 xbj[ii]=xbj[ii-1]+dBin[mm][ii-1];
     }
  
     TString TreeName="T";
     TChain* T;
     Double_t totalNe[10]={0.0};
     Double_t RawNe[10]={0.0};
     Double_t totalQ2[10]={0.0};
     Double_t totalXbj[10]={0.0};
     Double_t totalNe_err[10]={0.0};

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
/*         TCut ACC_phy;
         if(kin[mm]<6 || kin[mm]==15)
           ACC_phy=Form("abs(EKLxe.angle*180.0/3.14159-%f)<=1.52",Theta_c[KKin]);
         else{ 
           if(tmp_runp==1)
              ACC_phy=Form("abs(EKLxe.angle*180.0/3.14159-%f)<=1.52",Theta_c[KKin]);
           else 
              ACC_phy=Form("abs(EKLxe.angle*180.0/3.14159-%f)<=1.52",Theta_c1[KKin-6]);
	 }
*/
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
         T->SetBranchStatus("EKLxe.nu",1);
         T->SetBranchStatus("EKLxe.x_bj",1);
         T->SetBranchStatus("EKLxe.Q2",1);

         Double_t aEprime=0.0,aTheta=0.0,axbj=0.0,aQ2=0.0,aNu=0.0;
         T->SetBranchAddress("L.gold.p",&aEprime);
         T->SetBranchAddress("EKLxe.angle",&aTheta);
         T->SetBranchAddress("EKLxe.nu",&aNu);
         T->SetBranchAddress("EKLxe.x_bj",&axbj);
         T->SetBranchAddress("EKLxe.Q2",&aQ2);

         Double_t Radcor=1.0;
	 Int_t NNe[10]={0};
         Int_t nentries=electron->GetN();
         for(int jj=0;jj<nentries;jj++){
	     T->GetEntry(electron->GetEntry(jj));
	     Double_t dTh=-100.0;
             if(kin[mm]<6 || kin[mm]==15)
	        dTh=aTheta*180.0/3.14159-Theta_c[KKin];
             else{ 
                if(tmp_runp==1)
	           dTh=aTheta*180.0/3.14159-Theta_c[KKin];
                else 
	           dTh=aTheta*180.0/3.14159-Theta_c1[KKin-6];
	     }
	     Double_t dEp=aNu-Nu_c;
             int pass_ACC=0;
//	     pass_ACC=SearchACC(dTh,dEp,dTheta,ACC_table,ACC_matrix);
//             if(pass_ACC==0)continue;

             for(int kk=0;kk<nBin[mm];kk++){
		 Double_t dxbj=axbj-xbj[kk];
                 if(dxbj<dBin[mm][kk] && dxbj>=0){
		    totalNe[kk]+=1.0/livetime;
		    NNe[kk]++;
                    RawNe[kk]++;
		    totalQ2[kk]+=aQ2;
		    totalXbj[kk]+=axbj;
		    break;
                 }
	     }
         } 

         for(int jj=0;jj<nBin[mm];jj++){
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


    Double_t rawYield[10]={0.0};
    Double_t rawYield_err[10]={0.0};
    Double_t avgQ2[10]={0.0};
    Double_t avgXbj[10]={0.0};
    for(int ii=0;ii<nBin[mm];ii++){
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
