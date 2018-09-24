#include <fstream>
#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "GetTrees.h"
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
     TString filename;
     cout<<"Input file name: "<<endl;
     cin>>filename;
     filename="/w/halla-scifs17exp/triton/Runlist/"+filename;
     //filename="/w/halla-scifs17exp/triton/hanjie/MARATHON/analysis/Yield/Runlist/"+filename;

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

     int run_number,success=0;
     Ssiz_t from=0;
     TString content;
     if(tmp.ReadLine(infile)){
        while(tmp.Tokenize(content,from,","))
         {
              run_number = atoi(content);
              success=Update(run_number,target,kin);
         } 
    } 

       infile.close();

}
