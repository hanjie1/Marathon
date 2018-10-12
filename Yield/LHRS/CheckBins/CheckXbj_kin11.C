#include "GetRunList.h"
#include "GetTrees.h"
#include "SetCut.h"

void CheckXbj_kin11()
{
   TString target[4]={"H1","D2","He3","H3"};
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

   ofstream myfile;
   myfile.open("CheckKin11.txt");
   
   for(int ii=3;ii<4;ii++){
    int maxkin;
    if(ii==0)maxkin=5;
    else maxkin=11;
    for(int jj=8;jj<9;jj++){   
     int nrun=0;
     vector<Int_t> runList;
     nrun=GetRunList(runList,kin[jj],target[ii]);
     cout<<nrun<<" runs are added "<<endl;
     if(nrun==0)exit(0);

    for(int nn=0;nn<runList.size();nn++){
     TChain *T;
     T=GetTree(runList[nn],kin[jj],"T"); 
  
     T->Draw(">>electron",trigger2+CK+Ep+beta+ACC+VZ+TRK);
     TEventList *electron;
     gDirectory->GetObject("electron",electron);
     T->SetEventList(electron);

     T->SetBranchStatus("*",0);
     T->SetBranchStatus("EKLx.nu",1);

     Double_t aNu=0.0;
     T->SetBranchAddress("EKLx.nu",&aNu);
     
     int flag=0;
     Int_t nentries=electron->GetN();
     for(int jj=0;jj<nentries;jj++){
         T->GetEntry(electron->GetEntry(jj));
         if(aNu>7.625)flag=1;
     }
     
     if(flag==1)myfile<<target[ii].Data()<<"  "<<runList[nn]<<endl;

     delete electron;
     delete T;
   }
   }
  }

  myfile.close();


}
