#include "GetRunList.h"
#include "GetTrees.h"
#include "SetErrorCons.h"
#include "CalcCharge.h"
#include "CalcBoilingCorr.h"

void CalcBoilingError()
{
     int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
     TString target[4]={"H1","D2","He3","H3"};

     for(int ii=0;ii<4;ii++){
        int maxkin=11;
	if(ii==0)maxkin=5;
	ofstream outfile;
	outfile.open(Form("Results/%s_boiling.dat",target[ii].Data()));
	for(int jj=0;jj<maxkin;jj++){
            int nrun=0;
     	    vector<vector<Int_t> > runList;
     	    nrun=GetRunList(runList,kin[jj],target[ii]);
     	    if(nrun==0)exit(0);

    	    int run_number=0;
	    Double_t totalQbC=0.0,totalQbC_err=0.0;
     	    for(int kk=0;kk<nrun;kk++){
                run_number = runList[kk][0];
                //run_number = 2490;
		cout<<run_number<<endl;
                Double_t NQe=CalcCharge(run_number,kin[jj]);

		Double_t I=0.0,bCorr=1.0;
		Double_t dbCorr=0.0; //square of the uncertainty
		CalcBoilingCorr(run_number,kin[jj],target[ii],I,dbCorr,bCorr);
		totalQbC += NQe*bCorr;
		totalQbC_err += NQe*NQe*dbCorr;
            }
	    totalQbC_err=sqrt(totalQbC_err);
            Double_t rel_err=totalQbC_err/totalQbC;
	    outfile<<kin[jj]<<"  "<<totalQbC<<"  "<<totalQbC_err<<"  "<<rel_err<<endl;
	}
	outfile.close();
    }

}
