#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_Dp()
{
    Double_t x[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString filename;
    filename="newbin/Ratio_Dp.dat";
    int totalN = ReadFile(filename,x,Ratio,Rerr,RadCor,kin);
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}
    
    cout<<"Number of points:  "<<totalN<<endl;
    for(int ii=0;ii<totalN;ii++){

	/* Endcup correction */
	Double_t tmp_ECC=1.0+TMath::Exp(ECCA_Dp*x[ii]+ECCB_Dp);
	Ratio1[ii]=Ratio[ii]*tmp_ECC;
	Rerr1[ii]=Rerr[ii]*tmp_ECC;

        /* positron correction */
	Double_t tmp_pH1=1.0-TMath::Exp(pA_H1*x[ii]+pB_H1);
	Double_t tmp_pD2=1.0-TMath::Exp(pA_D2*x[ii]+pB_D2);
	Ratio2[ii]=Ratio1[ii]*(tmp_pD2/tmp_pH1);
	Rerr2[ii]=Rerr1[ii]*(tmp_pD2/tmp_pH1);

	/* radiative correction */
	Ratio3[ii]=Ratio2[ii]*RadCor[ii];	
	Rerr3[ii]=Rerr2[ii]*RadCor[ii];
    }     

    Double_t binx[26]={0.0};
    for(int ii=0;ii<26;ii++) binx[ii]=0.15+ii*0.03;

    Double_t x_final[26]={0.0},Ratio_final[26]={0.0},Rerr_final[26]={0.0};

    TGraphErrors *gDp=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/newbin/Dp_final.dat");
    outfile<<"x     Ratio     Ratio_err"<<endl;
    for(int ii=0;ii<26;ii++){
	int nn=0;
        Double_t tmpY[5]={0.0},tmpYerr[5]={0.0},tmpX[5]={0.0};;
	for(int jj=0;jj<totalN;jj++){
	    if(x[jj]==0)continue;
	    if((x[jj]-binx[ii])<0.03&&(x[jj]-binx[ii])>0){
		tmpY[nn]=Ratio3[jj];
		tmpYerr[nn]=Rerr3[jj];
		tmpX[nn]=x[jj];
		nn++;
	    }
	}
	if(nn==0)continue;
	Double_t var=0.0,x_weight=0.0,R_weight=0.0;
        for(int kk=0;kk<nn;kk++){
            x_weight+=tmpX[kk]/(tmpYerr[kk]*tmpYerr[kk]);
	    var+=1.0/(tmpYerr[kk]*tmpYerr[kk]);
 	    R_weight+=tmpY[kk]/(tmpYerr[kk]*tmpYerr[kk]);
	}       
	x_final[ii]=x_weight/var;
	Ratio_final[ii]=R_weight/var;
	Rerr_final[ii]=sqrt(1.0/var);
 	outfile<<x_final[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<endl;
        gDp->SetPoint(ii,x_final[ii],Ratio_final[ii]);
        gDp->SetPointError(ii,0,Rerr_final[ii]);
    }

    outfile.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    gDp->SetMarkerStyle(8);
    gDp->SetMarkerColor(4);
    gDp->Draw("AP");
    gDp->SetTitle("D/p;xbj;");

}
