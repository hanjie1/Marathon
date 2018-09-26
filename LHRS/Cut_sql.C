#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"

using namespace std;

int update_cut(int run_number=0){
  TString ssCut="T2,LeftBCMev.current_dnew>4";

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
      query=Form("insert %s (run_number,cut) values(%d,%s);",table.Data(),run_number,ssCut.Data());
      result=Server->Query(query.Data());
      if(!result){
        cout<<"run "<<run_number<<" can't insert sql"<<endl;
        Server->Close();
        return 0;
       }
  }
  if(update){
      query=Form("update %s set cut=%s where run_number=%d",table.Data(),ssCut.Data(),run_number);
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
              success=update_cut(run_number);
         }
    }

       infile.close();

}
       


