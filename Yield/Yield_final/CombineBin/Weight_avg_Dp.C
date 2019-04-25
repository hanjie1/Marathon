#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_Dp()
{
    Double_t x[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0};  

    TString filename;
    int totalN = ReadFile(filename,x,Ratio,Rerr);
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}
    
    cout<<"Number of points:  "<<totalN<<endl;
    for(int ii=0;ii<totalN;ii++){

	/* Endcup correction */
	Double_t tmp_ECC=1.0+Exp(ECCA_Dp*x[ii]+ECCB_Dp);
	Ratio1[ii]=Ratio[ii]*tmp_ECC;

        /* positron correction */
	Double_t tmp_pH1=1.0-Exp(pA_H1*x[ii]+pB_H1);
	Double_t tmp_pD2=1.0-Exp(pA_D2*x[ii]+pB_D2);
	Ratio2[ii]=Ratio1[ii]*(tmp_pD2/tmp_pH1);
    }     

     

}
