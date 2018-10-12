#include "GetTrees.h"
#include "SetCut.h"
#include "SetCons.h"
#include "CalcLum.h"
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

void CalcRawYield(){
     TString filename;
     cout<<"Input file name: "<<endl;
     cin>>filename;
     filename="/w/halla-scifs17exp/triton/Runlist/"+filename;

     Double_t LUM=CalcLum(filename); //total luminosity get for this kinematics;
     cout<<"Get total Luminosity:  "<<LUM<<endl;

     ifstream infile;
     infile.open(filename);
     if(!infile.is_open()){cout<<"!!! run list file not found "<<endl;return 0;}

     
     TString tmp;
     TString target;
     int kin=0;
     if(tmp.ReadToken(infile))target=tmp;
     else{
          cout<<"No target type!!!"<<endl;
          exit(0);
     }

     if(tmp.ReadToken(infile))kin=atoi(tmp);
     else{
          cout<<"No kinematic!!!"<<endl;
          exit(0);
     }

     ofstream myfile;
     myfile.open(Form("Output_Yield/%s_kin%d.txt",target.Data(),kin));
     myfile<<"xbj  Yield"<<endl;

     vector<Int_t> runList;
     int run_number,nrun=0;
     Ssiz_t from=0;
     TString content;
     if(tmp.ReadLine(infile)){
        while(tmp.Tokenize(content,from,","))
         {
              run_number = atoi(content);
              if(run_number>0){runList.push_back(run_number);nrun++;}
         }
     }
 
     infile.close();

     TCut xbj[17];
     for(int ii=0;ii<17;ii++){
         double xmin=0.16+ii*0.02;
         double xmax=0.16+(ii+1)*0.02;
         xbj[ii]=Form("EKLx.x_bj>%f && EKLx.x_bj<%f",xmin,xmax);
     }


     TString TreeName="T";
     TChain* T;
     Double_t totalNe[17]={0.0};
     Double_t totalNe_err[17]={0.0};
     Double_t NNe=0.0;
     for(int ii=0;ii<nrun;ii++){
         run_number=runList[ii];
         TRI_VAR LT=GetLT(run_number);
         Double_t livetime=LT.value; 
         Double_t livetime_err=LT.err; 
         cout<<"Get LT:  "<<livetime<<"  "<<livetime_err<<endl;
         T=GetTree(run_number,TreeName);
         for(int jj=0;jj<17;jj++){
           NNe=T->GetEntries(trigger2+CK+Ep+beta+ACC+VZ+TRK+beamtrip+xbj[jj]);
           totalNe[jj]+=NNe/livetime;
           if(NNe!=0)totalNe_err[jj]+=sqrt(1.0/NNe+(livetime_err/livetime)*(livetime_err/livetime))*NNe/livetime;
         }
         delete T;
         cout<<"Get electron coutns"<<endl;
     }

    Double_t rawYield[17]={0.0};
    Double_t rawYield_err[17]={0.0};
    for(int ii=0;ii<17;ii++){
	rawYield[ii]=totalNe[ii]/LUM;
	rawYield_err[ii]=totalNe_err[ii]/LUM;
        double tmp_x=0.17+ii*0.02;
        myfile<<tmp_x<<"   "<<rawYield[ii]<<"   "<<rawYield_err[ii]<<endl;
    }

   myfile.close();
}
