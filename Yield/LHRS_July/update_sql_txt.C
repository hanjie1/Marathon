#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"

using namespace std;

int update_current(int run_number=0,int ncurrent=0,double I1=0,double I2=0,double I3=0){

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
      query=Form("insert %s (run_number,current,I1,I2,I3) values(%d,%d,%lf,%lf,%lf);",table.Data(),run_number,ncurrent,I1,I2,I3);
      result=Server->Query(query.Data());
      if(!result){
        cout<<"run "<<run_number<<" can't insert sql"<<endl;
        Server->Close();
        return 0;
       }
  }

  if(update){
      query=Form("update %s set current=%d,I1=%lf,I2=%lf,I3=%lf where run_number=%d",table.Data(),ncurrent,I1,I2,I3,run_number);
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

void update_sql_txt(){
     TString filename;
     filename="current_final.txt";

     ifstream infile;
     infile.open(filename);
     if(!infile.is_open()){cout<<"!!! run list file not found "<<endl;return 0;}

     int run_number,success=0;
     double I1=0,I2=0,I3=0,ncurrent=0;
     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(infile)){
        while(tmp.Tokenize(content,from,"  "))
         {    //cout<<content<<"  "<<nn<<endl;
              if(nn==0) run_number = atoi(content);
              if(nn==1) ncurrent = atoi(content);
              if(nn==2) I1 = atof(content);
              if(nn==3) I2 = atof(content);
              if(nn==4) I3 = atof(content);
	      nn++;
         } 
        success=update_current(run_number,ncurrent,I1,I2,I3);
        nn=0;
        I1=0;I2=0;I3=0;ncurrent=0;
        from=0;
     } 

       infile.close();

}
