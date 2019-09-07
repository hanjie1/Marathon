#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_H3D()
{
    Double_t x[MAXNUM]={0.0},Q2[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0},Ratio4[MAXNUM]={0.0},Ratio5[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0},Rerr4[MAXNUM]={0.0},Rerr5[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    Double_t CoulCor[MAXNUM]={0.0};
    Double_t Pos_err[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString filename;
    filename="combine_newbin/Xbj_sort_H3D.dat";
    int totalN = ReadFile(filename,x,Q2,Ratio,Rerr,RadCor,CoulCor,kin);
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

    ifstream infile3;
    infile3.open("BinCenter/H3D_BCfac.dat");
    from=0;
    int mm=0;
    Double_t Xbc[MAXNUM]={0.0},BCfactor[MAXNUM]={0.0};
    while(tmp.ReadLine(infile3)){
          tmp.Tokenize(content,from," ");
          Xbc[mm]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          BCfactor[mm]=atof(content.Data());
          from=0;
          mm++;
     }
    infile3.close();

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

	/* tritium decay correction */
        int KKin=-1;
        if(kin[ii]<6)KKin=kin[ii];
	if(kin[ii]==7)KKin=6; 
	if(kin[ii]==9)KKin=7; 
	if(kin[ii]==11)KKin=8; 
	if(kin[ii]==13)KKin=9; 
	if(kin[ii]==15)KKin=10;
	if(kin[ii]==16)KKin=11;

	if(abs(x[ii]-x_HeD[ii])>0.001){cout<<"point "<<ii<<" has issue!"<<endl;continue;}
   	Ratio3[ii]=Ratio2[ii]*totalQ[KKin]/(totalQ[KKin]-totalQfH[KKin])-Ratio_HeD[ii]*totalQfH[KKin]/(totalQ[KKin]-totalQfH[KKin]);
	Rerr3[ii]=sqrt(pow(totalQ[KKin]/(totalQ[KKin]-totalQfH[KKin]),2)*Rerr2[ii]*Rerr2[ii]
                 +pow(totalQfH[KKin]/(totalQ[KKin]-totalQfH[KKin]),2)*Rerr_HeD[ii]*Rerr_HeD[ii]); 

	/* radiative correction + Coulomb Correction*/
	Ratio4[ii]=Ratio3[ii]*RadCor[ii]*CoulCor[ii];	
	Rerr4[ii]=Rerr3[ii]*RadCor[ii]*CoulCor[ii];

        if(abs(Xbc[ii]-x[ii])>0.001){cout<<"Something wrong with BC factor!"<<endl;continue;}
        Ratio5[ii]=Ratio4[ii]*BCfactor[ii];
        Rerr5[ii]=Rerr4[ii]*BCfactor[ii];

        Double_t pH3_Var=exp(2.0*(pA_H3*x[ii]+pB_H3))*(pow(x[ii],2)*pH3_VA+pH3_VB+2.0*x[ii]*pH3_COV_AB);
        Double_t pD2_Var=exp(2.0*(pA_D2*x[ii]+pB_D2))*(pow(x[ii],2)*pD2_VA+pD2_VB+2.0*x[ii]*pD2_COV_AB);
        Pos_err[ii]=sqrt(pH3_Var/(tmp_pH3*tmp_pH3)+pD2_Var/(tmp_pD2*tmp_pD2))*Ratio5[ii]; //positron absolute error on ratio

	/* estimate H3 decay error */
	Double_t He_conta=totalQfH[KKin]/(totalQ[KKin]-totalQfH[KKin])/Ratio3[ii];
cout<<He_conta<<endl;
    }     

    Double_t Ratio_final[19]={0.0},Rerr_final[19]={0.0};
    Double_t Rerr_pos[19]={0.0};

    TGraphErrors *gH3D=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/newbin/H3D_final.dat");
    outfile<<"x     Ratio     Ratio_err    relative_err"<<endl;

    ofstream outfile1;
    outfile1.open("ERROR/H3D_error.dat");
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
        gH3D->SetPoint(ii,X_center[ii],Ratio_final[ii]);
        gH3D->SetPointError(ii,0,Rerr_final[ii]);
        outfile<<X_center[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<endl;
        outfile1<<X_center[ii]<<"  "<<Rerr_pos[ii]<<"  "<<Rerr_pos[ii]/Ratio_final[ii]<<endl;
    }
    outfile.close();
    outfile1.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    gH3D->SetMarkerStyle(8);
    gH3D->SetMarkerColor(4);
    gH3D->Draw("AP");
    gH3D->SetTitle("H3/D2;xbj;");
/*
    TCanvas *c2=new TCanvas("c2","c2",1500,1500);
    TGraphErrors *gH3D_kin=new TGraphErrors();
    nn=0;
    for(int ii=0;ii<19;ii++){
        int tmpN=nn+nBin[ii];
        for(int jj=nn;jj<tmpN;jj++){
          if(Ratio4[jj]==0)continue;
          gH3D_kin->SetPoint(jj,X_center[ii],Ratio4[jj]);
          gH3D_kin->SetPointError(jj,0,Rerr4[jj]);
          nn++;
        }
    }
    gH3D_kin->SetMarkerStyle(8);
    gH3D_kin->SetMarkerColor(4);
    gH3D_kin->Draw("AP");
*/

   int color[12]={1,2,3,4,6,7,8,9,46,30,12,38};
    TGraphErrors *gH3D_kin[12];
    for(int ii=0;ii<12;ii++){
        gH3D_kin[ii]=new TGraphErrors();
        gH3D_kin[ii]->SetMarkerStyle(8);
        gH3D_kin[ii]->SetMarkerColor(color[ii]);
        gH3D_kin[ii]->SetMarkerSize(1.7);
    }

    int KKin[12]={0,1,2,3,4,5,7,9,11,13,15,16};
    int num[12]={0};
    for(int ii=0;ii<MAXNUM;ii++){
        if(x[ii]==0)continue;
        for(int kk=0;kk<12;kk++){
            if(kin[ii]!=KKin[kk])continue;
            gH3D_kin[kk]->SetPoint(num[kk],x[ii],Ratio4[ii]);
            gH3D_kin[kk]->SetPointError(num[kk],0,Rerr4[ii]);
            num[kk]++;
        }
    }

   TCanvas *c2=new TCanvas("c2","c2",1500,1200);
   TMultiGraph *mg1=new TMultiGraph();
   for(int ii=0;ii<12;ii++)
        mg1->Add(gH3D_kin[ii]);
   mg1->Draw("AP");
   mg1->SetTitle(";Bjorken x;#sigma({}^{3}H)/#sigma({}^{2}H)");

   auto leg1=new TLegend(0.7,0.65,0.9,0.9);
   leg1->SetNColumns(3);
   for(int ii=0;ii<12;ii++)
      leg1->AddEntry(gH3D_kin[ii],Form("kin%d",KKin[ii]),"P");

   leg1->Draw();
   mg1->GetYaxis()->SetRangeUser(1.22,1.5);

   c2->Print("Plots/H3D_kin.pdf");


}
