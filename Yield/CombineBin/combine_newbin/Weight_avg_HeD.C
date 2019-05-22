#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_HeD()
{
    Double_t x[MAXNUM]={0.0},Q2[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0},Ratio4[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0},Rerr4[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    Double_t Pos_err[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString filename;
    filename="Xbj_sort_HeD.dat";
    int totalN = ReadFile(filename,x,Q2,Ratio,Rerr,RadCor,kin);
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}
   
    ifstream infile1;
    infile1.open("BinCenter/BCfactor_HeD.dat");
    Ssiz_t from=0;
    TString content,tmp;
    int mm=0;
    Double_t Xbc[MAXNUM]={0.0},BCfactor[MAXNUM]={0.0};
    while(tmp.ReadLine(infile1)){
          tmp.Tokenize(content,from," ");
          Xbc[mm]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          BCfactor[mm]=atof(content.Data());
          from=0;
          mm++;
     }
    infile1.close();

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

	outfile1<<x[ii]<<"  "<<Ratio3[ii]<<"  "<<Rerr3[ii]<<"  "<<kin[ii]<<endl;

        if(abs(Xbc[ii]-x[ii])>0.001){cout<<"Something wrong with BC factor!"<<endl;continue;}
        Ratio4[ii]=Ratio3[ii]/BCfactor[ii];
        Rerr4[ii]=Rerr3[ii]/BCfactor[ii];

        Double_t pHe3_Var=exp(2.0*(pA_He3*x[ii]+pB_He3))*(pow(x[ii],2)*pHe3_VA+pHe3_VB+2.0*x[ii]*pHe3_COV_AB);
        Double_t pD2_Var=exp(2.0*(pA_D2*x[ii]+pB_D2))*(pow(x[ii],2)*pD2_VA+pD2_VB+2.0*x[ii]*pD2_COV_AB);
        Pos_err[ii]=sqrt(pHe3_Var/(tmp_pHe3*tmp_pHe3)+pD2_Var/(tmp_pD2*tmp_pD2))*Ratio4[ii]; //positron absolute error on ratio
    }     
    outfile1.close();

    Double_t Ratio_final[19]={0.0},Rerr_final[19]={0.0};
    Double_t Rerr_pos[19]={0.0};

    TGraphErrors *gHeD=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/newbin/HeD_final.dat");
    outfile<<"x     Ratio     Ratio_err    relative_err"<<endl;

    ofstream outfile2;
    outfile2.open("ERROR/HeD_error.dat");
    outfile2<<"x   positron_err  relative_err"<<endl;

    int nn=0;
    for(int ii=0;ii<19;ii++){
        int tmpN=nn+nBin[ii];
        Double_t var=0.0;
        Double_t tmpR=0.0;
        Double_t Epos_weight=0.0;
        for(int jj=nn;jj<tmpN;jj++){
          if(Ratio4[jj]==0)continue;
          var=var+1.0/(Rerr4[jj]*Rerr4[jj]);
          tmpR=tmpR+Ratio4[jj]/(Rerr4[jj]*Rerr4[jj]);
          Epos_weight+=pow(Pos_err[jj],2)/pow(Rerr4[jj],4);
          nn++;
        }
        if(var==0.0)continue;
        Ratio_final[ii]=tmpR/var;
        Rerr_final[ii]=1.0/sqrt(var);
        Rerr_pos[ii]=sqrt(Epos_weight)/var;
    }

    for(int ii=0;ii<19;ii++){
        if(Ratio_final[ii]==0)continue;
        gHeD->SetPoint(ii,X_center[ii],Ratio_final[ii]);
        gHeD->SetPointError(ii,0,Rerr_final[ii]);
        outfile<<X_center[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<endl;
        outfile2<<X_center[ii]<<"  "<<Rerr_pos[ii]<<"  "<<Rerr_pos[ii]/Ratio_final[ii]<<endl;
    }
    outfile.close();
    outfile2.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    gHeD->SetMarkerStyle(8);
    gHeD->SetMarkerColor(4);
    gHeD->Draw("AP");
    gHeD->SetTitle("He3/D2;xbj;");

    TCanvas *c2=new TCanvas("c2","c2",1500,1500);
    TGraphErrors *gHeD_kin=new TGraphErrors();
    nn=0;
    for(int ii=0;ii<19;ii++){
        int tmpN=nn+nBin[ii];
        for(int jj=nn;jj<tmpN;jj++){
          if(Ratio4[jj]==0)continue;
          gHeD_kin->SetPoint(jj,X_center[ii],Ratio4[jj]);
          gHeD_kin->SetPointError(jj,0,Rerr4[jj]);
          nn++;
        }
    }
    gHeD_kin->SetMarkerStyle(8);
    gHeD_kin->SetMarkerColor(4);
    gHeD_kin->Draw("AP");


}
