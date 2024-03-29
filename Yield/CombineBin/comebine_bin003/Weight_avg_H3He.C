#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_H3He()
{
    Double_t x[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0},Ratio4[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0},Rerr4[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    Double_t Pos_err[MAXNUM]={0.0};
    Double_t ECC_err[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString filename;
    filename="bin003/Ratio_H3He.dat";
    int totalN = ReadFile(filename,x,Ratio,Rerr,RadCor,kin);
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}
   
    ifstream infile1;
    infile1.open("TriDedcay.dat");
    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    Double_t totalQ[12]={0.0},totalQfH[12]={0.0};
    while(tmp.ReadLine(infile1)){
	  if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          totalQ[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          totalQfH[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    infile1.close();

    cout<<"Number of points:  "<<totalN<<endl;
    for(int ii=0;ii<totalN;ii++){

	/* Endcup correction */

	Double_t tmp_ECC=1.0+TMath::Exp(ECCA_H3He*x[ii]+ECCB_H3He);
	Ratio1[ii]=Ratio[ii]*tmp_ECC;
	Rerr1[ii]=Rerr[ii]*tmp_ECC;

        /* positron correction */
	Double_t tmp_pHe3=1.0-TMath::Exp(pA_He3*x[ii]+pB_He3);
	Double_t tmp_pH3=1.0-TMath::Exp(pA_H3*x[ii]+pB_H3);
	Ratio2[ii]=Ratio1[ii]*(tmp_pH3/tmp_pHe3);
	Rerr2[ii]=Rerr1[ii]*(tmp_pH3/tmp_pHe3);

	/* radiative correction */
	Ratio3[ii]=Ratio2[ii]*RadCor[ii];	
	Rerr3[ii]=Rerr2[ii]*RadCor[ii];

	/* tritium decay correction */
        int KKin=-1;
        if(kin[ii]<6)KKin=kin[ii];
	if(kin[ii]==7)KKin=6; 
	if(kin[ii]==9)KKin=7; 
	if(kin[ii]==11)KKin=8; 
	if(kin[ii]==13)KKin=9; 
	if(kin[ii]==15)KKin=10;
	if(kin[ii]==16)KKin=11;

   	Ratio4[ii]=Ratio3[ii]*totalQ[KKin]/(totalQ[KKin]-totalQfH[KKin])-totalQfH[KKin]/(totalQ[KKin]-totalQfH[KKin]);
	Rerr4[ii]=(totalQ[KKin]/(totalQ[KKin]-totalQfH[KKin]))*Rerr3[ii];

	/*  positron absolute error */
        Double_t pHe3_Var=exp(2.0*(pA_He3*x[ii]+pB_He3))*(pow(x[ii],2)*pHe3_VA+pHe3_VB+2.0*x[ii]*pHe3_COV_AB);
        Double_t pH3_Var=exp(2.0*(pA_H3*x[ii]+pB_H3))*(pow(x[ii],2)*pH3_VA+pH3_VB+2.0*x[ii]*pH3_COV_AB);
        Pos_err[ii]=sqrt(pHe3_Var/(tmp_pHe3*tmp_pHe3)+pH3_Var/(tmp_pH3*tmp_pH3))*Ratio4[ii]; //positron absolute error on ratio

	/*  End Cap absolute error */
        ECC_err[ii]=exp(ECCA_H3He*x[ii]+ECCB_H3He)*sqrt(x[ii]*x[ii]*ECCVA_H3He+ECCVB_H3He+2.0*ECC_CovH3He*x[ii]);

    }     

    Double_t binx[26]={0.0};
    for(int ii=0;ii<26;ii++) binx[ii]=0.15+ii*0.03;

    Double_t x_final[26]={0.0},Ratio_final[26]={0.0},Rerr_final[26]={0.0};
    Double_t Rerr_pos[26]={0.0},Rerr_ECC[26]={0.0};

    TGraphErrors *gH3He=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/bin003/H3He_final.dat");
    outfile<<"x     Ratio     Ratio_err      relative_err"<<endl;

    ofstream outfile1;
    outfile1.open("ERROR/H3He_error.dat");
    outfile1<<"x   e+_err    e+_rel_err     ECC_err     ECC_rel_err"<<endl;

    int mm=0;
    for(int ii=0;ii<26;ii++){
	int nn=0;
        Double_t tmpY[5]={0.0},tmpYerr[5]={0.0},tmpX[5]={0.0};;
        Double_t tmpEpos[5]={0.0},tmpE_ECC[5]={0.0};
	for(int jj=0;jj<totalN;jj++){
	    if(x[jj]==0)continue;
	    if((x[jj]-binx[ii])<0.03&&(x[jj]-binx[ii])>0){
		tmpY[nn]=Ratio4[jj];
		tmpYerr[nn]=Rerr4[jj];
		tmpX[nn]=x[jj];
                tmpEpos[nn]=Pos_err[jj];
                tmpE_ECC[nn]=ECC_err[jj];
		nn++;
	    }
	}
	if(nn==0)continue;
	Double_t var=0.0,x_weight=0.0,R_weight=0.0;
        Double_t Epos_weight=0.0,E_ECCweight=0.0;
        for(int kk=0;kk<nn;kk++){
            x_weight+=tmpX[kk]/(tmpYerr[kk]*tmpYerr[kk]);
	    var+=1.0/(tmpYerr[kk]*tmpYerr[kk]);
 	    R_weight+=tmpY[kk]/(tmpYerr[kk]*tmpYerr[kk]);
            Epos_weight+=pow(tmpEpos[kk],2)/pow(tmpYerr[kk],4);
            E_ECCweight+=pow(tmpE_ECC[kk],2)/pow(tmpYerr[kk],4);
	}       
	x_final[ii]=x_weight/var;
	Ratio_final[ii]=R_weight/var;
	Rerr_final[ii]=sqrt(1.0/var);
        Rerr_pos[ii]=sqrt(Epos_weight)/var;
        Rerr_ECC[ii]=sqrt(E_ECCweight)/var;
        outfile1<<x_final[ii]<<"  "<<Rerr_pos[ii]<<"  "<<Rerr_pos[ii]/Ratio_final[ii]<<"  "
                <<Rerr_ECC[ii]<<"  "<<Rerr_ECC[ii]/Ratio_final[ii]<<endl;
 	outfile<<x_final[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<endl;
        gH3He->SetPoint(mm,x_final[ii],Ratio_final[ii]);
        gH3He->SetPointError(mm,0,Rerr_final[ii]);
	mm++;
    }
    outfile1.close();
    outfile.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    gH3He->SetMarkerStyle(8);
    gH3He->SetMarkerColor(4);
    gH3He->Draw("AP");
    gH3He->SetTitle("H3/He3;xbj;");

}
