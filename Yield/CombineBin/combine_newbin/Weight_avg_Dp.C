#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_Dp()
{
    Double_t x[MAXNUM]={0.0},Q2[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0},Q2_bin[MAXNUM]={0.0};;
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0},Ratio4[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0},Rerr4[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    Double_t CoulCor[MAXNUM]={0.0};
    Double_t ECCor[MAXNUM]={0.0};
    Double_t Pos_err[MAXNUM]={0.0};
    Double_t ECC_err[MAXNUM]={0.0};
    Double_t relTar=0.79/100.0; //relative uncertainty from target thickness uncertainty
    Double_t R_Boil[5]={0.327/100.0,0.33/100.0,0.335/100.0,0.324/100.,0.406/100.}; //rel uncertainty from boiling on each kin
    Double_t relBoil[MAXNUM]={0.0};//rel uncer from boiling
    Double_t relACC=0.2/100.0; //rel uncertainty from acceptance cut
    Double_t relECC=0.3/100.0; //rel uncer from End cap correction
    Double_t relRC=0.26/100.0;  //rel uncer from radiative correction
    Double_t relBCC=0.08/100.0;  //rel uncer from bin centering 
    Double_t rel_totSys[MAXNUM]={0.0}; //rel total sys err
    int kin[MAXNUM]={0};

    TString filename;
    filename="combine_newbin/Xbj_sort_Dp.dat";
    int totalN = ReadFile(filename,x,Q2,Ratio,Rerr,RadCor,CoulCor,kin);
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}
   
    ifstream infile1;
    infile1.open("BinCenter/Dp_BCfac.dat");
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
    for(int ii=0;ii<totalN;ii++){

	/* Endcup correction */
	ECCor[ii]=1.0+TMath::Exp(ECCA_Dp*x[ii]+ECCB_Dp);
	Ratio1[ii]=Ratio[ii]*ECCor[ii];
	Rerr1[ii]=Rerr[ii]*ECCor[ii];

        /* positron correction */
	Double_t tmp_pH1=1.0-TMath::Exp(pA_H1*x[ii]+pB_H1);
	Double_t tmp_pD2=1.0-TMath::Exp(pA_D2*x[ii]+pB_D2);
	Ratio2[ii]=Ratio1[ii]*(tmp_pD2/tmp_pH1);
	Rerr2[ii]=Rerr1[ii]*(tmp_pD2/tmp_pH1);
//cout<<tmp_pD2/tmp_pH1<<endl;
	/* radiative correction + Coulomb correction*/
	Ratio3[ii]=Ratio2[ii]*RadCor[ii]*CoulCor[ii];	
	Rerr3[ii]=Rerr2[ii]*RadCor[ii]*CoulCor[ii];

	if(abs(Xbc[ii]-x[ii])>0.001){cout<<"Something wrong with BC factor!"<<endl;continue;}
	Ratio4[ii]=Ratio3[ii]*BCfactor[ii];
	Rerr4[ii]=Rerr3[ii]*BCfactor[ii];

	/* positron absolute error */
        Double_t pH1_Var=exp(2.0*(pA_H1*x[ii]+pB_H1))*(pow(x[ii],2)*pH1_VA+pH1_VB+2.0*x[ii]*pH1_COV_AB);
        Double_t pD2_Var=exp(2.0*(pA_D2*x[ii]+pB_D2))*(pow(x[ii],2)*pD2_VA+pD2_VB+2.0*x[ii]*pD2_COV_AB);
        Pos_err[ii]=sqrt(pH1_Var/(tmp_pH1*tmp_pH1)+pD2_Var/(tmp_pD2*tmp_pD2))*Ratio4[ii]; //positron absolute error on ratio

	/* End cap absolute error */
        ECC_err[ii]=exp(ECCA_Dp*x[ii]+ECCB_Dp)*sqrt(x[ii]*x[ii]*ECCVA_Dp+ECCVB_Dp+2.0*ECC_CovDp*x[ii]);
    }     

    Double_t Ratio_final[8]={0.0},Rerr_final[8]={0.0};
    Double_t Rerr_pos[8]={0.0},Rerr_ECC[8]={0.0};

    TGraphErrors *gDp=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/newbin/Dp_final.dat");
    outfile<<"x    Q2     Ratio     stat_err    sys_err    rel_stat_err    rel_sys_err    tot    rel_tot    Norm(rel)"<<endl;
   
    ofstream outfile1;
    outfile1.open("../Results/newbin/Dp_final_long.dat");
    outfile1<<"x  Q2   R   stat(rel)    ACC(rel)   boil(rel)   EC(rel)  RC(rel)   BC(rel)   sys(rel)    Norm(rel)"<<endl;
 
    int nn=0;
    for(int ii=0;ii<8;ii++){
        int tmpN=nn+nBin_Dp[ii];
	Double_t var=0.0;
	Double_t tmpR=0.0;
	Double_t tmpQ2=0.0;
        Double_t Epos_weight=0.0,E_ECCweight=0.0, E_boil=0.0;
        for(int jj=nn;jj<tmpN;jj++){
	  if(Ratio4[jj]==0)continue;
	  Double_t wi=1.0/(Rerr4[jj]*Rerr4[jj]);
	  var=var+wi;
	  tmpR=tmpR+Ratio4[jj]*wi;
	  tmpQ2=tmpQ2+Q2[jj]*wi;
	//  E_ECCweight+=pow(wi*Ratio4[jj]*0.3/100.0,2);
	  E_boil+=pow(wi*R_Boil[kin[jj]]*Ratio4[jj],2);
          Epos_weight+=pow(Pos_err[jj],2)/pow(Rerr4[jj],4);
	  nn++;
        }
	if(var==0.0)continue;
	Ratio_final[ii]=tmpR/var;
	Q2_bin[ii]=tmpQ2/var;
	Rerr_final[ii]=1.0/sqrt(var);
	relBoil[ii]=sqrt(E_boil)/var/Ratio_final[ii];	
	Rerr_pos[ii]=sqrt(Epos_weight)/var;
//	Rerr_ECC[ii]=sqrt(E_ECCweight)/var;
	rel_totSys[ii]=relACC*relACC+relBoil[ii]*relBoil[ii]+relECC*relECC+relRC*relRC+relBCC*relBCC;
	rel_totSys[ii]=sqrt(rel_totSys[ii]);
    }

    ofstream outfile3;
    outfile3.open("forThesis.dat");
    for(int ii=0;ii<8;ii++){
	if(Ratio_final[ii]==0)continue;
        gDp->SetPoint(ii,X_center_Dp[ii],Ratio_final[ii]);
	Double_t totalE=sqrt(pow(Rerr_final[ii]/Ratio_final[ii],2)+rel_totSys[ii]*rel_totSys[ii])*Ratio_final[ii];
        gDp->SetPointError(ii,0,totalE);
	outfile<<fixed<<setprecision(4);
        outfile<<X_center_Dp[ii]<<"  "<<Q2_bin[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<rel_totSys[ii]*Ratio_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<"   "<<rel_totSys[ii]<<"   "<<totalE<<"   "<<totalE/Ratio_final[ii]<<"  "<<relTar<<endl;
	outfile1<<setprecision(4);
        outfile1<<X_center_Dp[ii]<<"  "<<Q2_bin[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<"  "<<"  "<<relACC<<"  "<<relBoil[ii]<<"  "<<relECC<<"  "<<relRC<<"  "<<relBCC<<"  "<<rel_totSys[ii]<<"   "<<relTar<<endl;
        outfile3<<fixed<<setprecision(2)<<X_center_Dp[ii]<<" & "<<setprecision(2)<<Q2_bin[ii]<<" & "<<setprecision(4)<<Ratio_final[ii]<<" & "<<Rerr_final[ii]/Ratio_final[ii]<<" & "<<rel_totSys[ii]<<" \\\\"<<endl;
        outfile3<<"\\hline"<<endl;

    }

    outfile1.close();
    outfile3.close();
    outfile.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1200);
    gDp->SetMarkerStyle(8);
    gDp->SetMarkerColor(4);
    gDp->SetMarkerSize(2);
    gDp->Draw("AP");
    gDp->SetTitle(";Bjorken x;#sigma({}^{2}H)/#sigma({}^{1}H)");
    c1->Print("Plots/Dp_final.pdf");
/*
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
*/
    int color[5]={1,2,8,4,6};    
    TGraphErrors *gDp_kin[5];
    for(int ii=0;ii<5;ii++){
	gDp_kin[ii]=new TGraphErrors();
   	gDp_kin[ii]->SetMarkerStyle(8);
   	gDp_kin[ii]->SetMarkerColor(color[ii]);
   	gDp_kin[ii]->SetMarkerSize(1.7);
    }

    int num[5]={0}; 
    for(int ii=0;ii<MAXNUM;ii++){
	if(x[ii]==0)continue;
	for(int kk=0;kk<5;kk++){
	    if(kin[ii]!=kk)continue;
            gDp_kin[kk]->SetPoint(num[kk],x[ii],Ratio3[ii]);
            gDp_kin[kk]->SetPointError(num[kk],0,Rerr3[ii]);
            num[kk]++;	    
	}
    }

   TCanvas *c2=new TCanvas("c2","c2",1500,1200);
   TMultiGraph *mg1=new TMultiGraph();
   for(int ii=0;ii<5;ii++)
        mg1->Add(gDp_kin[ii]);
   mg1->Draw("AP");
   mg1->SetTitle(";Bjorken x;#sigma({}^{2}H)/#sigma({}^{1}H)");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   //leg1->SetNColumns(2);
   for(int ii=0;ii<5;ii++)
      leg1->AddEntry(gDp_kin[ii],Form("kin%d",ii),"P");
   
   leg1->Draw();

   c2->Print("Plots/Dp_kin.pdf");


}
