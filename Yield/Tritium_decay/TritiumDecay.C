#include "SetCons.h"
#include "GetTrees.h"
#include "GetRunList.h"
#include "CalcCharge.h"
#include "GetTime.h"

void TritiumDecay(){
     int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

     ofstream file1;
     file1.open("run_info.dat");
     file1<<"run No.    Charge     decay_time     fHi      Q*fHi"<<endl;

     ofstream file2;
     file2.open("TriDedcay.dat");
     file2<<"kin      total_charge      total_He3"<<endl;

     for(int ii=0;ii<11;ii++){
         vector<Int_t> runList;
         int run_number=0,nrun=0;
         nrun=GetRunList(runList,kin[ii],"H3");
         cout<<nrun<<" runs are added "<<endl;
         if(nrun==0)exit(0);
	 Double_t kin_prod=0.0,totalQ=0.0;;

	 for(int jj=0;jj<nrun;jj++){
	     run_number=runList[jj];
	     TChain *T;
	     T=GetTree(run_number,kin[ii],"T");
	     Int_t time=GetTime(T);
	     cout<<run_number<<"  "<<time<<endl;
	     Double_t Qj=CalcCharge(run_number,kin[ii]);
 	     Double_t fHj=(He3_t0+H3_t0*(1.0-exp(-1.0*time/tau)))/(He3_t0+H3_t0);
	     Double_t prodj=Qj*fHj;
	     file1<<run_number<<"  "<<Qj<<"  "<<time<<"  "<<fHj<<"  "<<prodj<<endl;

	     totalQ=totalQ+Qj;
	     kin_prod=kin_prod+prodj;
	 }

	 file2<<kin[ii]<<"  "<<totalQ<<"  "<<kin_prod<<endl; 
     }
     
     file1.close();
     file2.close();
}
