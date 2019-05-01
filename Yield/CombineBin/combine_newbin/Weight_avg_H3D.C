#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_H3D()
{
    Double_t x[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0},Ratio4[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0},Rerr4[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString filename;
    filename="bin003/Ratio_H3D.dat";
    int totalN = ReadFile(filename,x,Ratio,Rerr,RadCor,kin);
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}
   
    ifstream infile1;
    infile1.open("TriDedcay.dat");
    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    Double_t totalQ[11]={0.0},totalQfH[11]={0.0};
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
    infile2.open("CorrRatio_HeD.dat");
    Double_t x_HeD[MAXNUM]={0.0},Ratio_HeD[MAXNUM]={0.0},Rerr_HeD[MAXNUM]={0.0};
    int kin_HeD[MAXNUM]={0};
    from=0;
    nn=0;
    while(tmp.ReadLine(infile2)){
          tmp.Tokenize(content,from,"  ");
          x_HeD[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          Ratio_HeD[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          Rerr_HeD[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          kin_HeD[nn]=atoi(content.Data());
          from=0;
          nn++;
     }
    infile2.close();

    cout<<"Number of points:  "<<totalN<<endl;
    for(int ii=0;ii<totalN;ii++){

	/* Endcup correction */

	Double_t tmp_ECC=(1.0+TMath::Exp(ECCA_H3He*x[ii]+ECCB_H3He))*(1.0-TMath::Exp(ECCA_HeD*x[ii]+ECCB_HeD));
	Ratio1[ii]=Ratio[ii]*tmp_ECC;
	Rerr1[ii]=Rerr[ii]*tmp_ECC;

        /* positron correction */
	Double_t tmp_pD2=1.0-TMath::Exp(pA_D2*x[ii]+pB_D2);
	Double_t tmp_pH3=1.0-TMath::Exp(pA_H3*x[ii]+pB_H3);
	Ratio2[ii]=Ratio1[ii]*(tmp_pH3/tmp_pD2);
	Rerr2[ii]=Rerr1[ii]*(tmp_pH3/tmp_pD2);

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

	if(abs(x[ii]-x_HeD[ii])>0.001){cout<<"point "<<ii<<" has issue!"<<endl;continue;}
   	Ratio4[ii]=Ratio3[ii]*totalQ[KKin]/(totalQ[KKin]-totalQfH[KKin])-Ratio_HeD[ii]*totalQfH[KKin]/(totalQ[KKin]-totalQfH[KKin]);
	Rerr4[ii]=sqrt(pow(totalQ[KKin]/(totalQ[KKin]-totalQfH[KKin]),2)*Rerr3[ii]*Rerr3[ii]
                 +pow(totalQfH[KKin]/(totalQ[KKin]-totalQfH[KKin]),2)*Rerr_HeD[ii]*Rerr_HeD[ii]); 
    }     

    Double_t binx[26]={0.0};
    for(int ii=0;ii<26;ii++) binx[ii]=0.15+ii*0.03;

    Double_t x_final[26]={0.0},Ratio_final[26]={0.0},Rerr_final[26]={0.0};

    TGraphErrors *gH3D=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/bin003/H3D_final.dat");
    outfile<<"x     Ratio     Ratio_err      relative_err"<<endl;
    for(int ii=0;ii<26;ii++){
	int nn=0;
        Double_t tmpY[5]={0.0},tmpYerr[5]={0.0},tmpX[5]={0.0};;
	for(int jj=0;jj<totalN;jj++){
	    if(x[jj]==0)continue;
	    if((x[jj]-binx[ii])<0.03&&(x[jj]-binx[ii])>0){
		tmpY[nn]=Ratio4[jj];
		tmpYerr[nn]=Rerr4[jj];
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
 	outfile<<x_final[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<endl;
        gH3D->SetPoint(ii,x_final[ii],Ratio_final[ii]*2.0/3.0);
        gH3D->SetPointError(ii,0,Rerr_final[ii]);
    }

    outfile.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    gH3D->SetMarkerStyle(8);
    gH3D->SetMarkerColor(4);
    gH3D->Draw("AP");
    gH3D->SetTitle("H3/D2;xbj;");

}
