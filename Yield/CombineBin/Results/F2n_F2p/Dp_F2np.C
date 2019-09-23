#include "ReadFile.h"
#include "libEMCR.h"
#include "libNP.h"
void Dp_F2np(){
     Double_t x[MAXBIN]={0.0},Ratio[MAXBIN]={0.0},Rerr[MAXBIN]={0.0},F2np[MAXBIN]={0.0},F2np_err[MAXBIN]={0.0};
     Double_t x1[MAXBIN]={0.0},Q2[MAXBIN]={0.0},F2np1[MAXBIN]={0.0},F2np_err1[MAXBIN]={0.0};
     Double_t x_KP[7][110]={0.0},F2n_KP[7][110]={0.0},F2p_KP[7][110]={0.0};
     Double_t x_CJ[MAXBIN]={0.0},RD[MAXBIN]={0.0};
     
     TString Rfile="newbin/Dp_final.dat";
     int nbin=ReadYield(Rfile,x,Ratio,Rerr);
     Rfile="Model/NP_CJ.out";
     int nbin1=ReadModel(Rfile,x1,Q2,F2np1,F2np_err1);
     Rfile="Model/RD_CJ.out";
     int nbin_CJ=ReadCJ(Rfile,x_CJ,RD);
     Rfile="Model/F2dis_os0tm0ht0mec0_Dav18_He3Salme";
     int nbin2=ReadKP(Rfile,x_KP[0],F2p_KP[0],F2n_KP[0]);
     Rfile="Model/F2dis_os0tm0ht0mec1_Dav18_He3Salme";
     nbin2=ReadKP(Rfile,x_KP[1],F2p_KP[1],F2n_KP[1]);
     Rfile="Model/F2dis_os0tm1ht1mec1_Dav18_He3Salme";
     nbin2=ReadKP(Rfile,x_KP[2],F2p_KP[2],F2n_KP[2]);
     Rfile="Model/F2dis_os1tm0ht0mec1_Dav18_He3Salme";
     nbin2=ReadKP(Rfile,x_KP[3],F2p_KP[3],F2n_KP[3]);
     Rfile="Model/F2dis_os1tm1ht0mec1_Dav18_He3Salme";
     nbin2=ReadKP(Rfile,x_KP[4],F2p_KP[4],F2n_KP[4]);
     Rfile="Model/F2dis_os1tm1ht1mec0_Dav18_He3Salme";
     nbin2=ReadKP(Rfile,x_KP[5],F2p_KP[5],F2n_KP[5]);
     Rfile="Model/F2dis_os1tm1ht1mec1_Dav18_He3Salme";
     nbin2=ReadKP(Rfile,x_KP[6],F2p_KP[6],F2n_KP[6]);

     TGraphErrors *gDp=new TGraphErrors();
     TGraphErrors *gDp_CJ=new TGraphErrors(); //N/P extract with CJ R
     TGraphErrors *gDp1=new TGraphErrors(); //CJ N/P
     TGraphErrors *gDp2=new TGraphErrors(); //NMC N/P
     TGraph *gDp_KP[7]; //KP N/P

     for(int ii=0;ii<nbin;ii++){
	Double_t tmpR=D2_EMC(x[ii]);
	F2np[ii]=Ratio[ii]/tmpR-1.0;
	F2np_err[ii]=Rerr[ii]/tmpR;
	gDp->SetPoint(ii,x[ii],F2np[ii]);
	gDp->SetPointError(ii,0.0,F2np_err[ii]);
     }

     for(int ii=0;ii<nbin_CJ;ii++){
	Double_t tmpNP=Ratio[ii]/RD[ii]-1.0;
	Double_t tmpNPerr=Rerr[ii]/RD[ii];
	gDp_CJ->SetPoint(ii,x[ii],tmpNP);
	gDp_CJ->SetPointError(ii,0.0,tmpNPerr);
     }

     for(int ii=0;ii<nbin1;ii++){
	gDp1->SetPoint(ii,x1[ii],F2np1[ii]);
	gDp1->SetPointError(ii,0.0,F2np_err1[ii]);
	Double_t tmpNP=NMC_NP(x1[ii],Q2[ii]);
//	Double_t tmpNP=NMC_NP(x1[ii],14.0*x1[ii]);
	gDp2->SetPoint(ii,x1[ii],tmpNP);
	gDp2->SetPointError(ii,0.0,0.0);
     }

     for(int ii=0;ii<7;ii++){
	 gDp_KP[ii]=new TGraph();
	 int nn=0;
	 for(int jj=0;jj<nbin2;jj++){
	     if(x_KP[ii][jj]<0.15 || x_KP[ii][jj]>0.4)continue;
             gDp_KP[ii]->SetPoint(nn,x_KP[ii][jj],F2n_KP[ii][jj]/F2p_KP[ii][jj]);
	     nn++;
	 }
         gDp_KP[ii]->SetLineStyle(2);
         gDp_KP[ii]->SetLineColor(ii+30);
         gDp_KP[ii]->SetLineWidth(2);
     }

     TF1 *f1=new TF1("f1","SLAC_NP(x)",0.17,0.4);

     TCanvas *c1=new TCanvas("c1","c1",1500,1200);
     TMultiGraph *mg=new TMultiGraph();
     gDp->SetMarkerStyle(8);
     gDp->SetMarkerColor(2);
     gDp_CJ->SetMarkerStyle(8);
     gDp_CJ->SetMarkerColor(4);
     gDp1->SetLineStyle(8);
     gDp1->SetLineColor(4);
     gDp1->SetLineWidth(2);
     gDp2->SetLineStyle(8);
     gDp2->SetLineColor(8);
     gDp2->SetLineWidth(2);
     f1->SetLineStyle(1);
     f1->SetLineColor(4);
     mg->Add(gDp,"P");
     mg->Add(gDp_CJ,"P");
     mg->Add(gDp1,"L");
     mg->Add(gDp2,"L");
     mg->Add(gDp_KP[6],"L");
     mg->Draw("A");
     mg->SetTitle(";Bjorken x;F_{2}^{n}/F_{2}^{p}");
     f1->Draw("same");
     mg->GetYaxis()->SetRangeUser(0.5,0.9);

   auto leg1=new TLegend(0.55,0.75,0.9,0.9);
   leg1->SetNColumns(2);
   leg1->AddEntry(gDp,"#scale[1]{MARATHON (KP)}","P");
   leg1->AddEntry(gDp_CJ,"#scale[1]{MARATHON (CJ)}","P");
   leg1->AddEntry(gDp1,"#scale[1]{CJ15}","L");
   leg1->AddEntry(gDp2,"#scale[1]{NMC}","L");
   leg1->AddEntry(gDp_KP[6],"#scale[1]{KP}","L");
   leg1->AddEntry(f1,"#scale[1]{SLAC}","L");
   leg1->Draw();

     c1->Print("NP_Dp.pdf");
}
