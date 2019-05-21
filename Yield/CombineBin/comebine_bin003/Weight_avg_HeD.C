#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_HeD()
{
    Double_t x[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    Double_t Pos_err[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString filename;
    filename="bin003/Ratio_HeD.dat";
    int totalN = ReadFile(filename,x,Ratio,Rerr,RadCor,kin);
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}
    
    cout<<"Number of points:  "<<totalN<<endl;
    ofstream outfile1;
    outfile1.open("CorrRatio_HeD.dat");
    for(int ii=0;ii<totalN;ii++){

	/* Endcup correction */
	Double_t tmp_ECC=1.0-TMath::Exp(ECCA_HeD*x[ii]+ECCB_HeD);
	Ratio1[ii]=Ratio[ii]*tmp_ECC;
	Rerr1[ii]=Rerr[ii]*tmp_ECC;

        /* positron correction */
	Double_t tmp_pD2=1.0-TMath::Exp(pA_D2*x[ii]+pB_D2);
	Double_t tmp_pHe3=1.0-TMath::Exp(pA_He3*x[ii]+pB_He3);
	Ratio2[ii]=Ratio1[ii]*(tmp_pHe3/tmp_pD2);
	Rerr2[ii]=Rerr1[ii]*(tmp_pHe3/tmp_pD2);

	/* radiative correction */
	Ratio3[ii]=Ratio2[ii]*RadCor[ii];	
	Rerr3[ii]=Rerr2[ii]*RadCor[ii];

        Double_t pHe3_Var=exp(2.0*(pA_He3*x[ii]+pB_He3))*(pow(x[ii],2)*pHe3_VA+pHe3_VB+2.0*x[ii]*pHe3_COV_AB);
        Double_t pD2_Var=exp(2.0*(pA_D2*x[ii]+pB_D2))*(pow(x[ii],2)*pD2_VA+pD2_VB+2.0*x[ii]*pD2_COV_AB);
        Pos_err[ii]=sqrt(pHe3_Var/(tmp_pHe3*tmp_pHe3)+pD2_Var/(tmp_pD2*tmp_pD2))*Ratio3[ii]; //positron absolute error on ratio   
	outfile1<<x[ii]<<"  "<<Ratio3[ii]<<"  "<<Rerr3[ii]<<"  "<<kin[ii]<<endl;
    }     
    outfile1.close();
    Double_t binx[26]={0.0};
    for(int ii=0;ii<26;ii++) binx[ii]=0.15+ii*0.03;

    Double_t x_final[26]={0.0},Ratio_final[26]={0.0},Rerr_final[26]={0.0};
    Double_t Rerr_pos[26]={0.0};

    TGraphErrors *gHeD=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/bin003/HeD_final.dat");
    outfile<<"x     Ratio     Ratio_err      relative_err"<<endl;

    ofstream outfile2;
    outfile2.open("ERROR/HeD_error.dat");
    outfile2<<"x   positron_err  relative_err"<<endl;

    int mm=0;
    for(int ii=0;ii<26;ii++){
	int nn=0;
        Double_t tmpY[5]={0.0},tmpYerr[5]={0.0},tmpX[5]={0.0};;
        Double_t tmpEpos[5]={0.0};
	for(int jj=0;jj<totalN;jj++){
	    if(x[jj]==0)continue;
	    if((x[jj]-binx[ii])<0.03&&(x[jj]-binx[ii])>0){
		tmpY[nn]=Ratio3[jj];
		tmpYerr[nn]=Rerr3[jj];
		tmpX[nn]=x[jj];
                tmpEpos[nn]=Pos_err[jj];
		nn++;
	    }
	}
	if(nn==0)continue;
	Double_t var=0.0,x_weight=0.0,R_weight=0.0;
	Double_t Epos_weight=0.0;
        for(int kk=0;kk<nn;kk++){
            x_weight+=tmpX[kk]/(tmpYerr[kk]*tmpYerr[kk]);
	    var+=1.0/(tmpYerr[kk]*tmpYerr[kk]);
 	    R_weight+=tmpY[kk]/(tmpYerr[kk]*tmpYerr[kk]);
            Epos_weight+=pow(tmpEpos[kk],2)/pow(tmpYerr[kk],4);
	}       
	x_final[ii]=x_weight/var;
	Ratio_final[ii]=R_weight/var;
	Rerr_final[ii]=sqrt(1.0/var);
        Rerr_pos[ii]=sqrt(Epos_weight)/var;
        outfile2<<x_final[ii]<<"  "<<Rerr_pos[ii]<<"  "<<Rerr_pos[ii]/Ratio_final[ii]<<endl;
 	outfile<<x_final[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<endl;
        gHeD->SetPoint(mm,x_final[ii],Ratio_final[ii]);
        gHeD->SetPointError(mm,0,Rerr_final[ii]);
	mm++;
    }
    outfile2.close();
    outfile.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    gHeD->SetMarkerStyle(8);
    gHeD->SetMarkerColor(4);
    gHeD->Draw("AP");
    gHeD->SetTitle("He3/D2;xbj;");

}
