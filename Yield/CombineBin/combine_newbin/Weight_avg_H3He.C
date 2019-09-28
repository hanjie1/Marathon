#include <TMath.h>
#include "CorrFactor.h"
#include "ReadFile.h"
void Weight_avg_H3He()
{
    Double_t x[MAXNUM]={0.0},Q2[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0};
    Double_t Ratio1[MAXNUM]={0.0},Ratio2[MAXNUM]={0.0},Ratio3[MAXNUM]={0.0},Ratio4[MAXNUM]={0.0},Ratio5[MAXNUM]={0.0};  
    Double_t Rerr1[MAXNUM]={0.0},Rerr2[MAXNUM]={0.0},Rerr3[MAXNUM]={0.0},Rerr4[MAXNUM]={0.0},Rerr5[MAXNUM]={0.0};  
    Double_t RadCor[MAXNUM]={0.0};
    Double_t CoulCor[MAXNUM]={0.0};
    Double_t Pos_err[MAXNUM]={0.0};
    Double_t ECC_err[MAXNUM]={0.0};
    Double_t TriDecay_err[MAXNUM]={0.0};
    Double_t relTar=1.44/100.0; //relative uncertainty from target thickness uncertainty
    Double_t R_Boil[12]={0.299/100.0,0.283/100.0,0.286/100.0,0.277/100.,0.37/100.,0.418/100.0,0.355/100.,0.361/100.,0.329/100.,0.328/100.,0.321/100.,0.317/100.}; //rel uncertainty from boiling on each kin
    Double_t relBoil[MAXNUM]={0.0};//rel uncer from boiling
    Double_t relACC=0.2/100.0; //rel uncertainty from acceptance cut
    Double_t relECC=0.3/100.0; //rel uncer from End cap correction
    Double_t relRC=0.35/100.0;  //rel uncer from radiative correction
    Double_t relBCC=0.0; //rel uncer from bin centering
    Double_t relH3dec=0.0;  //rel uncer from H3 decay
    Double_t rel_totSys[MAXNUM]={0.0}; //rel total sys err

    int kin[MAXNUM]={0};

    TString filename;
    filename="combine_newbin/Xbj_sort_H3He.dat";
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
    infile2.open("BinCenter/H3He_BCfac.dat");
    from=0;
    int mm=0;
    Double_t Xbc[MAXNUM]={0.0},BCfactor[MAXNUM]={0.0};
    while(tmp.ReadLine(infile2)){
          tmp.Tokenize(content,from," ");
          Xbc[mm]=atof(content.Data());
          tmp.Tokenize(content,from," ");
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

	/* tritium decay correction */
        int KKin=-1;
        if(kin[ii]<6)KKin=kin[ii];
	if(kin[ii]==7)KKin=6; 
	if(kin[ii]==9)KKin=7; 
	if(kin[ii]==11)KKin=8; 
	if(kin[ii]==13)KKin=9; 
	if(kin[ii]==15)KKin=10;
	if(kin[ii]==16)KKin=11;

   	Ratio3[ii]=Ratio2[ii]*totalQ[KKin]/(totalQ[KKin]-totalQfH[KKin])-totalQfH[KKin]/(totalQ[KKin]-totalQfH[KKin]);
	Rerr3[ii]=(totalQ[KKin]/(totalQ[KKin]-totalQfH[KKin]))*Rerr2[ii];

	/* radiative correction + Coulomb correction */
	Ratio4[ii]=Ratio3[ii]*RadCor[ii]*CoulCor[ii];	
	Rerr4[ii]=Rerr3[ii]*RadCor[ii]*CoulCor[ii];


        if(abs(Xbc[ii]-x[ii])>0.001){cout<<"Something wrong with BC factor!"<<endl;continue;}
        Ratio5[ii]=Ratio4[ii]*BCfactor[ii];
        Rerr5[ii]=Rerr4[ii]*BCfactor[ii];

	/* positron absolute error */
        Double_t pHe3_Var=exp(2.0*(pA_He3*x[ii]+pB_He3))*(pow(x[ii],2)*pHe3_VA+pHe3_VB+2.0*x[ii]*pHe3_COV_AB);
        Double_t pH3_Var=exp(2.0*(pA_H3*x[ii]+pB_H3))*(pow(x[ii],2)*pH3_VA+pH3_VB+2.0*x[ii]*pH3_COV_AB);
        Pos_err[ii]=sqrt(pHe3_Var/(tmp_pHe3*tmp_pHe3)+pH3_Var/(tmp_pH3*tmp_pH3))*Ratio5[ii]; //positron absolute error on ratio

	/* End cap absolute error */
        ECC_err[ii]=exp(ECCA_H3He*x[ii]+ECCB_H3He)*sqrt(x[ii]*x[ii]*ECCVA_H3He+ECCVB_H3He+2.0*ECC_CovH3He*x[ii]);	

	/* Tritium decay relative error */
	Double_t tmp_Ratio=Ratio2[ii]*totalQ[KKin]/(totalQ[KKin]-totalQfH[KKin]*(1+0.00174))-totalQfH[KKin]*(1+0.00174)/(totalQ[KKin]-totalQfH[KKin]*(1+0.00174));
	TriDecay_err[ii]=abs(tmp_Ratio-Ratio3[ii])/Ratio3[ii];
    }     

    Double_t Ratio_final[19]={0.0},Rerr_final[19]={0.0};
    Double_t Rerr_pos[19]={0.0},Rerr_ECC[19]={0.0},Rerr_H3decay[19]={0.0};

    TGraphErrors *gH3He=new TGraphErrors();
    ofstream outfile;
    outfile.open("../Results/newbin/H3He_final.dat");
    outfile<<"x    Q2     Ratio     stat_err    sys_err    rel_stat_err    rel_sys_err    tot    rel_tot   Norm(rel)"<<endl;

    ofstream outfile2;
    outfile2.open("../Results/newbin/H3He_final_long.dat");
    outfile2<<"x  Q2   R   stat(rel)   ACC(rel)   boil(rel)   EC(rel)  RC(rel)   BC(rel)   H3dec(rel)   sys(rel)   Norm(rel)"<<endl;


    nn=0;
    for(int ii=0;ii<19;ii++){
        int tmpN=nn+nBin[ii];
        Double_t var=0.0;
        Double_t tmpR=0.0;
        Double_t Epos_weight=0.0,E_ECCweight=0.0, E_H3decay=0.0,E_boil=0.0;
        for(int jj=nn;jj<tmpN;jj++){
          if(Ratio5[jj]==0)continue;
          int KKin=-1;
          if(kin[jj]<6)KKin=kin[jj];
          if(kin[jj]==7)KKin=6;
          if(kin[jj]==9)KKin=7;
          if(kin[jj]==11)KKin=8;
          if(kin[jj]==13)KKin=9;
          if(kin[jj]==15)KKin=10;
          if(kin[jj]==16)KKin=11;

          Double_t wi=1.0/(Rerr5[jj]*Rerr5[jj]);
          var=var+wi;
          tmpR=tmpR+Ratio5[jj]*wi;
          E_boil+=pow(wi*R_Boil[KKin]*Ratio5[jj],2);

          Epos_weight+=pow(Pos_err[jj],2)/pow(Rerr5[jj],4);
          E_ECCweight+=pow(ECC_err[jj],2)/pow(Rerr5[jj],4);
          E_H3decay+=pow(TriDecay_err[jj]*Ratio[5],2)/pow(Rerr5[jj],4);
          nn++;
        }
        if(var==0.0)continue;
        Ratio_final[ii]=tmpR/var;
        Rerr_final[ii]=1.0/sqrt(var);
        Rerr_pos[ii]=sqrt(Epos_weight)/var;
        Rerr_ECC[ii]=sqrt(E_ECCweight)/var;
        Rerr_H3decay[ii]=sqrt(E_H3decay)/var;

        relBoil[ii]=sqrt(E_boil)/var/Ratio_final[ii];

        rel_totSys[ii]=relACC*relACC+relBoil[ii]*relBoil[ii]+relECC*relECC+relRC*relRC+relBCC*relBCC+relH3dec*relH3dec;
        rel_totSys[ii]=sqrt(rel_totSys[ii]);
    }

    ofstream outfile3;
    outfile3.open("forThesis.dat");
    for(int ii=0;ii<19;ii++){
        if(Ratio_final[ii]==0)continue;
        gH3He->SetPoint(ii,X_center[ii],Ratio_final[ii]);
        Double_t totalE=sqrt(pow(Rerr_final[ii]/Ratio_final[ii],2)+rel_totSys[ii]*rel_totSys[ii])*Ratio_final[ii];
        gH3He->SetPointError(ii,0,totalE);
        outfile<<setprecision(4);
        outfile<<X_center[ii]<<"  "<<Q2[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]<<"  "<<rel_totSys[ii]*Ratio_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<"   "<<rel_totSys[ii]<<"   "<<totalE<<"   "<<totalE/Ratio_final[ii]<<"  "<<relTar<<endl;
        outfile2<<setprecision(4);
        outfile2<<X_center[ii]<<"  "<<Q2[ii]<<"  "<<Ratio_final[ii]<<"  "<<Rerr_final[ii]/Ratio_final[ii]<<"  "<<relACC<<"  "<<relBoil[ii]<<"  "<<relECC<<"  "<<relRC<<"  "<<relBCC<<"  "<<"  "<<relH3dec<<"  "<<rel_totSys[ii]<<"  "<<relTar<<endl;
        outfile3<<fixed<<setprecision(2)<<X_center[ii]<<" & "<<setprecision(2)<<Q2[ii]<<" & "<<setprecision(4)<<Ratio_final[ii]<<" & "<<Rerr_final[ii]/Ratio_final[ii]<<" & "<<rel_totSys[ii]<<" \\\\"<<endl;
        outfile3<<"\\hline"<<endl;

    }
    outfile.close();
    outfile2.close();
    outfile3.close();

    TCanvas *c1=new TCanvas("c1","c1",1500,1200);
    gH3He->SetMarkerStyle(8);
    gH3He->SetMarkerColor(4);
    gH3He->Draw("AP");
    gH3He->SetMarkerSize(2);
    gH3He->SetTitle(";Bjorken x;#sigma({}^{3}H)/#sigma({}^{3}He)");
    c1->Print("Plots/H3He_final.pdf");


/*
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
*/
   int color[12]={1,2,3,4,6,7,8,9,46,30,12,38};
    TGraphErrors *gH3He_kin[12];
    for(int ii=0;ii<12;ii++){
        gH3He_kin[ii]=new TGraphErrors();
        gH3He_kin[ii]->SetMarkerStyle(8);
        gH3He_kin[ii]->SetMarkerColor(color[ii]);
        gH3He_kin[ii]->SetMarkerSize(1.7);
    }

    int KKin[12]={0,1,2,3,4,5,7,9,11,13,15,16};
    int num[12]={0};
    for(int ii=0;ii<MAXNUM;ii++){
        if(x[ii]==0)continue;
        for(int kk=0;kk<12;kk++){
            if(kin[ii]!=KKin[kk])continue;
            gH3He_kin[kk]->SetPoint(num[kk],x[ii],Ratio4[ii]);
            gH3He_kin[kk]->SetPointError(num[kk],0,Rerr4[ii]);
            num[kk]++;
        }
    }

   TCanvas *c2=new TCanvas("c2","c2",1500,1200);
   TMultiGraph *mg1=new TMultiGraph();
   for(int ii=0;ii<12;ii++)
        mg1->Add(gH3He_kin[ii]);
   mg1->Draw("AP");
   mg1->SetTitle(";Bjorken x;#sigma({}^{3}H)/#sigma({}^{3}He)");

   auto leg1=new TLegend(0.7,0.65,0.9,0.9);
   leg1->SetNColumns(3);
   for(int ii=0;ii<12;ii++)
      leg1->AddEntry(gH3He_kin[ii],Form("kin%d",KKin[ii]),"P");

   leg1->Draw();
   //mg1->GetYaxis()->SetRangeUser(1.22,1.5);

   c2->Print("Plots/H3He_kin.pdf");



}
