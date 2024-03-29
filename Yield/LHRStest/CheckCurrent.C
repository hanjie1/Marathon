#include <fstream>
#include "GetTrees.h"
using namespace std;

int CheckCurrent(const int run_number,int kin)
{
     TString TreeName="T";
     TChain *T=GetTree(run_number,kin,TreeName);

     Double_t Current[10];
     for(int ii=0;ii<10;ii++){
	 Current[ii]=0;
     }

     TCut timeup="LeftBCMev.BeamUp_time_v1495[0]>60";
     T->Draw(">>beamup",timeup);
     TEventList *beamup;
     gDirectory->GetObject("beamup",beamup);
     T->SetEventList(beamup);

     T->SetBranchStatus("*",0);
     T->SetBranchStatus("LeftBCMev.current_dnew",1);

     Double_t c_dnew;
     T->SetBranchAddress("LeftBCMev.current_dnew",&c_dnew);

     Int_t nentries=beamup->GetN();
     cout<<nentries<<endl;
     Double_t lastcurrent=0;
     int index=0;
     for(int ii=0;ii<nentries;ii++){
        //cout<<T->GetFile()<<endl;
        Int_t tmp_en=beamup->GetEntry(ii);
        T->GetEntry(tmp_en);
        //cout<<T->GetFile()<<endl;
        if(abs(c_dnew-lastcurrent)>2){
           Current[index++]=c_dnew;
           lastcurrent=c_dnew;
        }
     }

     Double_t mark[10]={0}; 
     Double_t Current_final[10]={0};
     int kk=0;
     for(int ii=0;ii<10;ii++){
         //cout<<Current[ii]<<endl;
         if(Current[ii]==0||mark[ii]==1)continue;
         Current_final[kk]=Current[ii];
         kk++;
         for(int jj=ii+1;jj<10;jj++){
             if(mark[jj]==1)continue;
             if(abs(Current[jj]-Current[ii])<=1)mark[jj]=1;
          }
     }

     ofstream myfile;
     myfile.open("current.txt",fstream::app); 
     myfile<<run_number<<"  "<<kk<<"  ";
     for(int ii=0;ii<10;ii++){
	 if(Current_final[ii]==0)continue;
         myfile<<fixed<<setprecision(1)<<Current_final[ii]<<"  ";
     }
     myfile<<endl;
     delete T;
     return kk;
}
