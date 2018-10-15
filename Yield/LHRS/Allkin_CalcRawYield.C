#include "GetTrees.h"
#include "GetRunList.h"
#include "SetCut.h"
#include "SetCons.h"
#include "CalcLum.h"
#include "SearchXS.h"
#include <TMath.h>

TRI_VAR GetLT(int run_number)
{
     TString table="LHRStest";
     TSQLServer* Server = TSQLServer::Connect("mysql://halladb/triton-work","triton-user","3He3Hdata");

     TString  query=Form("select * from %s where run_number=%d;",table.Data(),run_number);
     TSQLResult* result=Server->Query(query.Data());
     TSQLRow *row;
     int nrows = result->GetRowCount();
     TRI_VAR LT;
     if (nrows==0) {
        cout<<"run "<<run_number<<" is NOT in "<<table<<endl; 
        Server->Close();
        LT.value=0;
        LT.err=0;   
        return LT;
     }

     query=Form("select livetime,LT_err from %s where run_number=%d;",table.Data(),run_number);
     result=Server->Query(query.Data());
     row=result->Next();
     Double_t livetime=atof(row->GetField(0));
     Double_t livetime_err=atof(row->GetField(1));

     LT.value=livetime;
     LT.err=livetime_err;

     Server->Close(); 
     return LT;      
}

void Allkin_CalcRawYield(){
   TString target[4]={"H1","D2","He3","H3"};
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

   for(int nn=0;nn<4;nn++){   
    for(int mm=0;mm<11;mm++){
     if(nn==0&&mm>4)break;
     Double_t LUM=CalcLum(kin[mm],target[nn]); //total luminosity get for this kinematics;
     cout<<"Get total Luminosity for target "<<target[nn]<<"  "<<" kin "<<kin[mm]<<" : "<<LUM<<endl;

    ofstream myfile;
    myfile.open(Form("RawYield/vz009_new_25per/%s_kin%d.txt",target[nn].Data(),kin[mm]));
    myfile<<"n   xbj   Q2   Yield   Yield_err"<<endl;

     vector<Int_t> runList;
     int run_number=0,nrun=0;
     nrun=GetRunList(runList,kin[mm],target[nn]);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

     Double_t xbj[35];
     Double_t dBin=(xmax[mm]-xmin[mm])/(nBin[mm]*1.0);
     for(int ii=0;ii<nBin[mm];ii++){
        // xbj[ii]=xmin[mm]+ii*dBin;
         xbj[ii]=0.14+ii*0.02;
     }
     cout<<endl;

     TString TreeName="T";
     TChain* T;
     Double_t totalNe[35]={0.0};
     Double_t RawNe[35]={0.0};
     Double_t totalQ2[35]={0.0};
     Double_t totalXbj[35]={0.0};
     Double_t totalNe_err[35]={0.0};
     for(int ii=0;ii<nrun;ii++){
         run_number=runList[ii];
         TRI_VAR LT=GetLT(run_number);
         Double_t livetime=LT.value; 
         Double_t livetime_err=LT.err; 
         cout<<"Get LT:  "<<livetime<<"  "<<livetime_err<<endl;
         T=GetTree(run_number,kin[mm],TreeName);
         T->Draw(">>electron",trigger2+CK+Ep+beta+ACC+VZ+TRK);
         TEventList *electron;
         gDirectory->GetObject("electron",electron);
         T->SetEventList(electron);

         T->SetBranchStatus("*",0);
         T->SetBranchStatus("L.gold.p",1);
         T->SetBranchStatus("EKLx.angle",1);
         T->SetBranchStatus("EKLx.x_bj",1);
         T->SetBranchStatus("EKLx.Q2",1);
         T->SetBranchStatus("LeftBCMev.isrenewed",1);

         Double_t aEprime=0.0,aTheta=0.0,axbj=0.0,aQ2=0.0;
         Double_t isrenewed=0;
         T->SetBranchAddress("L.gold.p",&aEprime);
         T->SetBranchAddress("EKLx.angle",&aTheta);
         T->SetBranchAddress("EKLx.x_bj",&axbj);
         T->SetBranchAddress("EKLx.Q2",&aQ2);
         T->SetBranchAddress("LeftBCMev.isrenewed",&isrenewed);

         Double_t Radcor=1.0;
	 Int_t NNe[35]={0};
         Int_t nentries=electron->GetN();
         for(int jj=0;jj<nentries;jj++){
	     T->GetEntry(electron->GetEntry(jj));
             for(int kk=0;kk<nBin[mm];kk++){
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

         for(int jj=0;jj<nBin[mm];jj++){
           Double_t tmp_err=0.0;
           if(NNe[jj]!=0){
              tmp_err=sqrt(1.0/NNe[jj]+(livetime_err/livetime)*(livetime_err/livetime))*NNe[jj]/livetime;
           }
	   totalNe_err[jj]+=tmp_err*tmp_err;
         }
         delete T;
         cout<<"Get electron coutns"<<endl;
     }


    Double_t rawYield[35]={0.0};
    Double_t rawYield_err[35]={0.0};
    Double_t avgQ2[35]={0.0};
    Double_t avgXbj[35]={0.0};
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
