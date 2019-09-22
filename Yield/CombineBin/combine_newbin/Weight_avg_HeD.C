#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_HeD()
{
    Double_t x[MAXNUM]={0.0},Q2[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0},Ratio4[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0},Rerr4[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    Double_t CoulCor[MAXNUM]={0.0};
    Double_t Pos_err[MAXNUM]={0.0};
    Double_t ECC_err[MAXNUM]={0.0};
    Double_t relTar=1.20/100.0; //relative uncertainty from target thickness uncertainty
    Double_t R_Boil[12]={0.311/100.0,0.312/100.0,0.322/100.0,0.306/100.,0.402/100.,0.47/100.0,0.411/100.,0.397/100.,0.39/100.,0.368/100.,0.36/100.,0.356/100.}; //rel uncertainty from boiling on each kin
    Double_t relBoil[MAXNUM]={0.0};//rel uncer from boiling
    Double_t relACC=0.2/100.0; //rel uncertainty from acceptance cut
    Double_t relECC=0.3/100.0; //rel uncer from End cap correction
    Double_t relRC=0.25/100.0;  //rel uncer from radiative correction
    Double_t R_BCC=0.025/100.0;  //x>=0.75 rel uncer from bin centering on individual bin
    Double_t relBCC[MAXNUM]={0.0}; //rel uncer from bin centering
    Double_t rel_totSys[MAXNUM]={0.0}; //rel total sys err

    int kin[MAXNUM]={0};

    TString filename;
    filename="combine_newbin/Xbj_sort_HeD.dat";
    int totalN = ReadFile(filename,x,Q2,Ratio,Rerr,RadCor,CoulCor,kin);
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}
   
    ifstream infile1;
    infile1.open("BinCenter/HeD_BCfac.dat");
    Ssiz_t from=0;
    TString content,tmp;
    int mm=0;
    Double_t Xbc[MAXNUM]={0.0},BCfactor[MAXNUM]={0.0};
    while(tmp.ReadLine(infile1)){
          tmp.Tokenize(content,from," ");
          Xbc[mm]=atof(content.Data());
          tmp.Tokenize(content,from," ");
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
//cout<<tmp_pHe3/tmp_pD2<<endl;
	/* radiative correction + Coulom correction*/
	Ratio3[ii]=Ratio2[ii]*RadCor[ii]*CoulCor[ii];	
	Rerr3[ii]=Rerr2[ii]*RadCor[ii]*CoulCor[ii];

	outfile1<<x[ii]<<"  "<<Ratio2[ii]<<"  "<<Rerr2[ii]<<"  "<<kin[ii]<<endl;

        if(abs(Xbc[ii]-x[ii])>0.001){cout<<"Something wrong with BC factor!"<<endl;continue;}
        Ratio4[ii]=Ratio3[ii]*BCfactor[ii];
        Rerr4[ii]=Rerr3[ii]*BCfactor[ii];

	/* positron absolute error */
        Double_t pHe3_Var=exp(2.0*(pA_He3*x[ii]+pB_He3))*(pow(x[ii],2)*pHe3_VA+pHe3_VB+2.0*x[ii]*pHe3_COV_AB);
        Double_t pD2_Var=exp(2.0*(pA_D2*x[ii]+pB_D2))*(pow(x[ii],2)*pD2_VA+pD2_VB+2.0*x[ii]*pD2_COV_AB);
        Pos_err[ii]=sqrt(pHe3_Var/(tmp_pHe3*tmp_pHe3)+pD2_Var/(tmp_pD2*tmp_pD2))*Ratio4[ii]; //positron absolute error on ratio

	/* End cap absolute error */
        ECC_err[ii]=exp(ECCA_HeD*x[ii]+ECCB_HeD)*sqrt(x[ii]*x[ii]*ECCVA_HeD+ECCVB_HeD+2.0*ECC_CovHeD*x[ii]);

    }     
    outfile1.close();

    Double_t Ratio_final[19]={0.0},Rerr_final[19]={0.0};
    Double_t Rerr_pos[19]={0.0},Rerr_ECC[19]={0.0};

    TGraphErrors *gHeD=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/newbin/HeD_final.dat");
    outfile<<"x    Q2     Ratio     stat_err    sys_err    rel_stat_err    rel_sys_err    tot    rel_tot    Norm(rel)"<<endl;

    ofstream outfile2;
    outfile2.open("../Results/newbin/HeD_final_long.dat");
    outfile2<<"x  Q2   R   stat(rel)   ACC(rel)   boil(rel)   EC(rel)  RC(rel)   BC(rel)   sys(rel)   Norm(rel)"<<endl;

    int nn=0;
    for(int ii=0;ii<19;ii++){
        int tmpN=nn+nBin[ii];
        Double_t var=0.0;
        Double_t tmpR=0.0;
        Double_t Epos_weight=0.0,E_ECCweight=0.0,E_boil=0.0,E_BCC=0.0;
        for(int jj=nn;jj<tmpN;jj++){
          if(Ratio4[jj]==0)continue;
          Double_t wi=1.0/(Rerr4[jj]*Rerr4[jj]);
          var=var+wi;
          tmpR=tmpR+Ratio4[jj]*wi;
          Epos_weight+=pow(Pos_err[jj],2)/pow(Rerr4[jj],4);
          E_ECCweight+=pow(ECC_err[jj],2)/pow(Rerr4[jj],4);

          E_boil+=pow(wi*R_Boil[kin[jj]]*Ratio4[jj],2);
	  if(x[jj]>=0.75){
             E_BCC+=pow(wi*R_BCC*Ratio4[jj],2);
	  }
          nn++;
        }
        if(var==0.0)continue;
        Ratio_final[ii]=tmpR/var;
        Rerr_final[ii]=1.0/sqrt(var);
        relBoil[ii]=sqrt(E_boil)/var/Ratio_final[ii];
        relBCC[ii]=sqrt(E_BCC)/var/Ratio_final[ii];
	
        rel_totSys[ii]=relACC*relACC+relBoil[ii]*relBoil[ii]+relECC*relECC+relRC*relRC+relBCC[ii]*relBCC[ii];
        rel_totSys[ii]=sqrt(rel_totSys[ii]);

        Rerr_pos[ii]=sqrt(Epos_weight)/var;
        Rerr_ECC[ii]=sqrt(E_ECCweight)/var;
    }

    ofstream outfile3;
    outfile3.open("forThesis.dat");
    for(int ii=0;ii<19;ii++){
        if(Ratio_final[ii]==0)continue;
        gHeD->SetPoint(ii,X_center[ii],Ratio_final[ii]);
        Double_t totalE=sqrt(pow(Rerr_final[ii]/Ratio_final[ii],2)+rel_totSys[ii]*rel_totSys[ii])*Ratio_final[ii];
        gHeD->SetPointError(ii,0,totalE);
        outfile<<fixed<<setprecision(4);
        outfile<<X_center[ii]<<"  "<<Q2[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<rel_totSys[ii]*Ratio_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<"   "<<rel_totSys[ii]<<"   "<<totalE<<"   "<<totalE/Ratio_final[ii]<<"  "<<relTar<<endl;
        outfile2<<setprecision(4);
        outfile2<<X_center[ii]<<"  "<<Q2[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<"  "<<relACC<<"  "<<relBoil[ii]<<"  "<<relECC<<"  "<<relRC<<"  "<<relBCC[ii]<<"  "<<rel_totSys[ii]<<"  "<<relTar<<endl;
	outfile3<<fixed<<setprecision(2)<<X_center[ii]<<" & "<<setprecision(2)<<Q2[ii]<<" & "<<setprecision(4)<<Ratio_final[ii]<<" & "<<Rerr_final[ii]/Ratio_final[ii]<<" & "<<rel_totSys[ii]<<" \\\\"<<endl;
	outfile3<<"\\hline"<<endl;

    }
    outfile.close();
    outfile2.close();
    outfile3.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1200);
    gHeD->SetMarkerStyle(8);
    gHeD->SetMarkerColor(4);
    gHeD->SetMarkerSize(2);
    gHeD->Draw("AP");
    gHeD->SetTitle(";Bjorken x;#sigma({}^{3}He)/#sigma({}^{2}H)");
    gHeD->GetYaxis()->SetRangeUser(1.5,1.75);
    c1->Print("Plots/HeD_final.pdf");

/*
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
*/
   int color[12]={1,2,3,4,6,7,8,9,46,30,12,38};
    TGraphErrors *gHeD_kin[12];
    for(int ii=0;ii<12;ii++){
        gHeD_kin[ii]=new TGraphErrors();
        gHeD_kin[ii]->SetMarkerStyle(8);
        gHeD_kin[ii]->SetMarkerColor(color[ii]);
        gHeD_kin[ii]->SetMarkerSize(1.7);
    }

    int KKin[12]={0,1,2,3,4,5,7,9,11,13,15,16};
    int num[12]={0};
    for(int ii=0;ii<MAXNUM;ii++){
        if(x[ii]==0)continue;
        for(int kk=0;kk<12;kk++){
            if(kin[ii]!=KKin[kk])continue;
            gHeD_kin[kk]->SetPoint(num[kk],x[ii],Ratio3[ii]);
            gHeD_kin[kk]->SetPointError(num[kk],0,Rerr3[ii]);
            num[kk]++;
        }
    }

   TCanvas *c2=new TCanvas("c2","c2",1500,1200);
   TMultiGraph *mg1=new TMultiGraph();
   for(int ii=0;ii<12;ii++)
        mg1->Add(gHeD_kin[ii]);
   mg1->Draw("AP");
   mg1->SetTitle(";Bjorken x;#sigma({}^{3}He)/#sigma({}^{2}H)");

   auto leg1=new TLegend(0.15,0.65,0.35,0.9);
   leg1->SetNColumns(3);
   for(int ii=0;ii<12;ii++)
      leg1->AddEntry(gHeD_kin[ii],Form("kin%d",KKin[ii]),"P");

   leg1->Draw();
   mg1->GetYaxis()->SetRangeUser(1.5,1.75);

   c2->Print("Plots/HeD_kin.pdf");




}
