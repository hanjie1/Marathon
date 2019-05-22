#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_H3He()
{
    Double_t x[MAXNUM]={0.0},Q2[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0},Ratio4[MAXNUM]={0.0},Ratio5[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0},Rerr4[MAXNUM]={0.0},Rerr5[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    Double_t Pos_err[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString filename;
    filename="Xbj_sort_H3He.dat";
    int totalN = ReadFile(filename,x,Q2,Ratio,Rerr,RadCor,kin);
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

    ifstream infile2;
    infile2.open("BinCenter/BCfactor_H3He.dat");
    from=0;
    int mm=0;
    Double_t Xbc[MAXNUM]={0.0},BCfactor[MAXNUM]={0.0};
    while(tmp.ReadLine(infile2)){
          tmp.Tokenize(content,from," ");
          Xbc[mm]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          BCfactor[mm]=atof(content.Data());
          from=0;
          mm++;
     }
    infile2.close();

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

        if(abs(Xbc[ii]-x[ii])>0.001){cout<<"Something wrong with BC factor!"<<endl;continue;}
        Ratio5[ii]=Ratio4[ii]/BCfactor[ii];
        Rerr5[ii]=Rerr4[ii]/BCfactor[ii];

        Double_t pHe3_Var=exp(2.0*(pA_He3*x[ii]+pB_He3))*(pow(x[ii],2)*pHe3_VA+pHe3_VB+2.0*x[ii]*pHe3_COV_AB);
        Double_t pH3_Var=exp(2.0*(pA_H3*x[ii]+pB_H3))*(pow(x[ii],2)*pH3_VA+pH3_VB+2.0*x[ii]*pH3_COV_AB);
        Pos_err[ii]=sqrt(pHe3_Var/(tmp_pHe3*tmp_pHe3)+pH3_Var/(tmp_pH3*tmp_pH3))*Ratio5[ii]; //positron absolute error on ratio
    }     

    Double_t Ratio_final[19]={0.0},Rerr_final[19]={0.0};
    Double_t Rerr_pos[19]={0.0};

    TGraphErrors *gH3He=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/newbin/H3He_final.dat");
    outfile<<"x     Ratio     Ratio_err    relative_err"<<endl;

    ofstream outfile1;
    outfile1.open("ERROR/H3He_error.dat");
    outfile1<<"x   positron_err  relative_err"<<endl;

    nn=0;
    for(int ii=0;ii<19;ii++){
        int tmpN=nn+nBin[ii];
        Double_t var=0.0;
        Double_t tmpR=0.0;
        Double_t Epos_weight=0.0;
        for(int jj=nn;jj<tmpN;jj++){
          if(Ratio5[jj]==0)continue;
          var=var+1.0/(Rerr5[jj]*Rerr5[jj]);
          tmpR=tmpR+Ratio5[jj]/(Rerr5[jj]*Rerr5[jj]);
          Epos_weight+=pow(Pos_err[jj],2)/pow(Rerr5[jj],4);
          nn++;
        }
        if(var==0.0)continue;
        Ratio_final[ii]=tmpR/var;
        Rerr_final[ii]=1.0/sqrt(var);
        Rerr_pos[ii]=sqrt(Epos_weight)/var;
    }

    for(int ii=0;ii<19;ii++){
        if(Ratio_final[ii]==0)continue;
        gH3He->SetPoint(ii,X_center[ii],Ratio_final[ii]);
        gH3He->SetPointError(ii,0,Rerr_final[ii]);
        outfile<<X_center[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<endl;
        outfile1<<X_center[ii]<<"  "<<Rerr_pos[ii]<<"  "<<Rerr_pos[ii]/Ratio_final[ii]<<endl;
    }
    outfile.close();
    outfile1.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    gH3He->SetMarkerStyle(8);
    gH3He->SetMarkerColor(4);
    gH3He->Draw("AP");
    gH3He->SetTitle("H3/He3;xbj;");

    TCanvas *c2=new TCanvas("c2","c2",1500,1500);
    TGraphErrors *gH3He_kin=new TGraphErrors();
    nn=0;
    for(int ii=0;ii<19;ii++){
        int tmpN=nn+nBin[ii];
        for(int jj=nn;jj<tmpN;jj++){
          if(Ratio4[jj]==0)continue;
          gH3He_kin->SetPoint(jj,X_center[ii],Ratio4[jj]);
          gH3He_kin->SetPointError(jj,0,Rerr4[jj]);
          nn++;
        }
    }
    gH3He_kin->SetMarkerStyle(8);
    gH3He_kin->SetMarkerColor(4);
    gH3He_kin->Draw("AP");

}
