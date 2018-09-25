#include <TMath.h>

const double bH3_A=0.0001293;
const double bH3_B=-0.007399;

const double bHe3_A=0.00008686;
const double bHe3_B=-0.004759;

const double bD2_A=0.0001147;
const double bD2_B=-0.006651;

const double bH1_A=0.0001527; 
const double bH1_B=-0.008529;

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
     double aCharge=atof(row->GetField(0));
     int ncurrent=atoi(row->GetField(1));
     //double I1=atof(row->GetField(2));
     double target_thickness=atof(row->GetField(3));
     TString target=row->GetField(4);

 
     Server->Close(); 
      
     Charge=aCharge/(Qe*1e6);

     TString TreeName="T";
     TChain* T=GetTree(run_number,kin,TreeName);

     Double_t gain=0.0003361;
     Double_t offset=0.0217;

     Double_t dnewr;
     T->SetBranchStatus("*",0);
     T->SetBranchStatus("evLeftdnew",1);

     T->SetBranchAddress("evLeftdnew_r",&dnewr);
 
     Int_t nentries=T->GetEntries();
     Double_t totalI=0;
     Int_t nI=0;
     for(int ii=0;ii<nentries;ii++)
      {
	 T->GetEntry(ii);
         Double_t tmpI = gain*dnewr+offset;
         if(tmpI>0 && tmpI<30){
            totalI+=tmpI;
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
//cout<<"current:  "<<avgI<<endl;
     return;
}

Double_t CalcLum(int kin, TString target){
     int nrun=0;
     vector<Int_t> runList;
     nrun=GetRunList(runList,kin,target);
     if(nrun==0)exit(0);

     int run_number,success=0;
     Ssiz_t from=0;
     TString content;
     Double_t NNtarg=0.0,Ncharge=0.0;
     Double_t LUM=0.0;
     for(int ii=0;ii<runList.size();ii++){
         run_number = runList[ii];
         RunLum(run_number,kin,Ncharge,NNtarg);
         LUM+=Ncharge*NNtarg/CMtoNB;
    }
    
    return LUM;

}






