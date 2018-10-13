#include <fstream>
#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "GetTrees.h"
#include "GetRunList.h"
#include "SetCut.h"
#include "SetCons.h"
#include "CheckCurrent.h"
#include "CalcLT.h"
#include "CalcCharge.h"

using namespace std;

int Update(int run_number=0,TString target="NULL",int kin=0){
  int beamcut=1;
  int ncurrent=0;
  ncurrent=CheckCurrent(run_number,kin);
  TRI_VAR livetime=CalcLT(run_number,kin,beamcut);
  Double_t charge=CalcCharge(run_number,kin);
  Double_t target_thickness;
  //int ncurrent=0;
  //ncurrent=CheckCurrent(run_number,kin);
  TString scut="T2,LeftBCMev.current_dnew>0";

  if(target=="D2")target_thickness=0.1422;
  if(target=="H1")target_thickness=0.0708;
  if(target=="He3")target_thickness=0.0534;
  if(target=="H3")target_thickness=0.077;
  if(target=="NULL"){
     cout<<"no target type"<<endl;
     return 0;
  }

  TString table="LHRStest";

  TSQLServer* Server = TSQLServer::Connect("mysql://halladb/triton-work","triton-user","3He3Hdata");
 
  TString  query=Form("select * from %s where run_number=%d;",table.Data(),run_number);
  TSQLResult* result=Server->Query(query.Data());
  TSQLRow *row;
  int nrows = result->GetRowCount();
  int update=0,insert=0;
  if (nrows==0)insert=1;
  else update=1;

  if(insert){
      query=Form("insert %s (run_number,target,kinematic,livetime,charge,target_thickness,current,LT_err,cut) values(%d,'%s',%d,%lf,%lf,%lf,%d,%lf,'%s');",table.Data(),run_number,target.Data(),kin,livetime.value,charge,target_thickness,ncurrent,livetime.err,scut.Data());
      result=Server->Query(query.Data());
      if(!result){
        cout<<"run "<<run_number<<" can't insert sql"<<endl;
        Server->Close();
        return 0;
       }
  }

  if(update){
//      query=Form("update %s set target='%s',kinematic=%d,livetime=%lf,charge=%lf,target_thickness=%lf,current=%d where run_number=%d",table.Data(),target.Data(),kin,livetime,charge,target_thickness,ncurrent,run_number);
      query=Form("update %s set charge=%lf,current=%d where run_number=%d",table.Data(),charge,ncurrent,run_number);
      result=Server->Query(query.Data());
      if(!result){
         cout<<"run "<<run_number<<" can't update sql"<<endl;
         Server->Close();
         return 0;
       }
  }

  Server->Close();
  return 1;
}

void Update_sql(){
     TString target;
     int kin;
     cout<<"Target:  ";
     cin>>target;
     cout<<"Kin:     ";
     cin>>kin;

     vector<Int_t> runList;
     int run_number=0,nrun=0;
     nrun=GetRunList(runList,kin,target);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

     int success=0;
     for(int ii=0;ii<runList.size();ii++) 
         {
             run_number = runList[ii];
             success=Update(run_number,target,kin);
         } 


}
