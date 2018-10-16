#include <TMath.h>
#include "CalcCharge.h"

const double bH3_A=0.00021;
const double bH3_B=-0.00893;

const double bHe3_A=0.00006;
const double bHe3_B=-0.00392;

const double bD2_A=0.00017;
const double bD2_B=-0.00813;

const double bH1_A=0.00022;
const double bH1_B=-0.009822;

const Double_t Qe=TMath::Qe();
const Double_t Na=TMath::Na();

const Double_t CMtoNB=1.0e33; 

void RunLum(int run_number,int kin,Double_t& Charge,Double_t& Ntarg)
{
     Charge=0;
     Ntarg=0;
     TString table="LHRStest";
 
     TSQLServer* Server = TSQLServer::Connect("mysql://halladb/triton-work","triton-user","3He3Hdata");

     TString  query=Form("select * from %s where run_number=%d;",table.Data(),run_number);
     TSQLResult* result=Server->Query(query.Data());
     TSQLRow *row;
     int nrows = result->GetRowCount();
     if (nrows==0) {
        cout<<"run "<<run_number<<" is NOT in "<<table<<endl; 
        Server->Close();
        return;
     }

     query=Form("select charge,current,I1,target_thickness,target from %s where run_number=%d;",table.Data(),run_number);
     result=Server->Query(query.Data());
     row=result->Next();
//     double aCharge=atof(row->GetField(0));
//     int ncurrent=atoi(row->GetField(1));
//     double I1=atof(row->GetField(2));
     double target_thickness=atof(row->GetField(3));
     TString target=row->GetField(4); 
     Server->Close(); 

     Double_t aCharge=CalcCharge(run_number,kin);

      
     TString TreeName="T";
     TChain* T=GetTree(run_number,kin,TreeName);

     Charge=aCharge/(Qe*1e6);
     
     T->SetBranchStatus("*",0);
     T->SetBranchStatus("LeftBCMev.BeamUp_time_v1495",1);
     T->SetBranchStatus("LeftBCMev.current_dnew",1);

     Double_t current_dnew;
     Double_t beamUp[5];
     T->SetBranchAddress("LeftBCMev.BeamUp_time_v1495",beamUp);
     T->SetBranchAddress("LeftBCMev.current_dnew",&current_dnew);
     
     Int_t nentries=T->GetEntries();
     Double_t totalI=0;
     Int_t nI=0;
     for(int ii=0;ii<nentries;ii++)
      {
	 T->GetEntry(ii);
         if(current_dnew>4){
            totalI+=current_dnew;
            nI++;
         }
      }
     
     Double_t avgI=0;
     if(nI!=0)avgI=totalI/nI;
    
     Double_t boiling_corr=1;
     Double_t massA=1.0;
     if(target=="H1") {
        boiling_corr=bH1_A*avgI*avgI+bH1_B*avgI+1.0;
        massA=1.0;
     }
     if(target=="D2"){
        boiling_corr=bD2_A*avgI*avgI+bD2_B*avgI+1.0;
        massA=2.0;
     }
     if(target=="He3"){
        boiling_corr=bHe3_A*avgI*avgI+bHe3_B*avgI+1.0;
        massA=3.0;
     }
     if(target=="H3"){
        boiling_corr=bH3_A*avgI*avgI+bH3_B*avgI+1.0;
        massA=3.0;
     }
     Ntarg=target_thickness*boiling_corr*Na/massA;

     return;
}

Double_t CalcLum(TString filename){
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
     Double_t NNtarg=0.0,Ncharge=0.0;
     Double_t LUM=0.0;
     if(tmp.ReadLine(infile)){
        while(tmp.Tokenize(content,from,","))
         {
              run_number = atoi(content);
              RunLum(run_number,kin,Ncharge,NNtarg);
              LUM+=Ncharge*NNtarg/CMtoNB;
         }
    }
       infile.close();
    
    return LUM;

}






