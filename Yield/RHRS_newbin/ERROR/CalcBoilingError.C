#include "GetRunList.h"
#include "GetTrees.h"
#include "SetErrorCons.h"

void CalcBoilingError()
{
     int kin[1]={16};
     TString target[4]={"H1","D2","He3","H3"};

     for(int ii=0;ii<4;ii++){
	ofstream outfile;
	outfile.open(Form("Results/%s_boiling.dat",target[ii].Data()));
	for(int jj=0;jj<1;jj++){
            int nrun=0;
     	    vector<vector<Int_t> > runList;
     	    nrun=GetRunList(runList,kin[jj],target[ii]);
     	    if(nrun==0)exit(0);

    	    int run_number=0;
	    Double_t avgI=0.0,totalI=0.0;
	    Int_t nI=0;
     	    for(int kk=0;kk<nrun;kk++){
                run_number = runList[kk][0];
		cout<<run_number<<endl;

                TChain* T=GetTree(run_number,kin[jj],"T");
		if(T==NULL)continue;
     		Double_t dnewr;
     		T->SetBranchStatus("*",0);
     		T->SetBranchStatus("evRightdnew",1);
     		T->SetBranchAddress("evRightdnew_r",&dnewr);

     		Int_t nentries=T->GetEntries();
     		for(int mm=0;mm<nentries;mm++){
         	    T->GetEntry(mm);
         	    Double_t tmpI = dnew_gain*dnewr+dnew_offset;
                    if(tmpI>dnew_offset && tmpI<30){
                       totalI+=tmpI;
                       nI++;
                    }
      		}
     		delete T;
            }
	    avgI=totalI/(nI*1.0);

	    Double_t boiling_corr=0.0,boiling_err=0.0;
            if(target[ii]=="H1") {
       	       boiling_corr=bH1_A*avgI*avgI+bH1_B*avgI+1.0;
               boiling_err=bH1_VC+pow(avgI,2)*bH1_VB+pow(avgI,4)*bH1_VA
                          +2.0*avgI*bH1_COV_BC+2.0*pow(avgI,2)*bH1_COV_AC+2.0*pow(avgI,3)*bH1_COV_AB;
     	    }
            if(target[ii]=="D2"){
               boiling_corr=bD2_A*avgI*avgI+bD2_B*avgI+1.0;
               boiling_err=bD2_VC+pow(avgI,2)*bD2_VB+pow(avgI,4)*bD2_VA
                          +2.0*avgI*bD2_COV_BC+2.0*pow(avgI,2)*bD2_COV_AC+2.0*pow(avgI,3)*bD2_COV_AB;
            }
            if(target[ii]=="He3"){
               boiling_corr=bHe3_A*avgI*avgI+bHe3_B*avgI+1.0;
               boiling_err=bHe3_VC+pow(avgI,2)*bHe3_VB+pow(avgI,4)*bHe3_VA
                          +2.0*avgI*bHe3_COV_BC+2.0*pow(avgI,2)*bHe3_COV_AC+2.0*pow(avgI,3)*bHe3_COV_AB;
            }
            if(target[ii]=="H3"){
               boiling_corr=bH3_A*avgI*avgI+bH3_B*avgI+1.0;
               boiling_err=bH3_VC+pow(avgI,2)*bH3_VB+pow(avgI,4)*bH3_VA
                          +2.0*avgI*bH3_COV_BC+2.0*pow(avgI,2)*bH3_COV_AC+2.0*pow(avgI,3)*bH3_COV_AB;
            }

	    Double_t rel_err=sqrt(boiling_err)/boiling_corr; 

	    outfile<<kin[jj]<<"  "<<avgI<<"  "<<rel_err<<endl;
	}
	outfile.close();
    }

}
