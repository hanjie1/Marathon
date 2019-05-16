#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_Dp()
{
    Double_t x[MAXNUM]={0.0},Q2[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0},Ratio4[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0},Rerr4[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString filename;
    filename="Xbj_sort_Dp.dat";
    int totalN = ReadFile(filename,x,Q2,Ratio,Rerr,RadCor,kin);
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}
   
    ifstream infile1;
    infile1.open("BinCenter/BCfactor_Dp.dat");
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

	if(abs(Xbc[ii]-x[ii])>0.001){cout<<"Something wrong with BC factor!"<<endl;continue;}
	Ratio4[ii]=Ratio3[ii]/BCfactor[ii];
	Rerr4[ii]=Rerr3[ii]/BCfactor[ii];
    }     

    Double_t Ratio_final[8]={0.0},Rerr_final[8]={0.0};

    TGraphErrors *gDp=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/newbin/Dp_final.dat");
    outfile<<"x     Ratio     Ratio_err    relative_err"<<endl;
    
    int nn=0;
    for(int ii=0;ii<8;ii++){
        int tmpN=nn+nBin_Dp[ii];
	Double_t var=0.0;
	Double_t tmpR=0.0;
        for(int jj=nn;jj<tmpN;jj++){
	  if(Ratio4[jj]==0)continue;
	  var=var+1.0/(Rerr4[jj]*Rerr4[jj]);
	  tmpR=tmpR+Ratio4[jj]/(Rerr4[jj]*Rerr4[jj]);
	  nn++;
        }
	if(var==0.0)continue;
	Ratio_final[ii]=tmpR/var;
	Rerr_final[ii]=1.0/sqrt(var);
    } 

   
    for(int ii=0;ii<8;ii++){
	if(Ratio_final[ii]==0)continue;
        gDp->SetPoint(ii,X_center_Dp[ii],Ratio_final[ii]);
        gDp->SetPointError(ii,0,Rerr_final[ii]);
        outfile<<X_center_Dp[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<endl;
    }


    outfile.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    gDp->SetMarkerStyle(8);
    gDp->SetMarkerColor(4);
    gDp->Draw("AP");
    gDp->SetTitle("D/p;xbj;");

    TCanvas *c2=new TCanvas("c2","c2",1500,1500);
    TGraphErrors *gDp_kin=new TGraphErrors();
    nn=0;
    for(int ii=0;ii<8;ii++){
        int tmpN=nn+nBin_Dp[ii];
        for(int jj=nn;jj<tmpN;jj++){
          if(Ratio4[jj]==0)continue;
          gDp_kin->SetPoint(jj,X_center_Dp[ii],Ratio4[jj]);
          gDp_kin->SetPointError(jj,0,Rerr4[jj]);
          nn++;
        }
    }
    gDp_kin->SetMarkerStyle(8);
    gDp_kin->SetMarkerColor(4);
    gDp_kin->Draw("AP"); 


}
