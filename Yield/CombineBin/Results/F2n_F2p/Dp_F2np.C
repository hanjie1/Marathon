#include "ReadFile.h"
#include "libEMCR.h"
#include "libNP.h"
void Dp_F2np(){
     Double_t x[MAXBIN]={0.0},Ratio[MAXBIN]={0.0},Rerr[MAXBIN]={0.0},F2np[MAXBIN]={0.0},F2np_err[MAXBIN]={0.0};
     Double_t x1[MAXBIN]={0.0},Q2[MAXBIN]={0.0},F2np1[MAXBIN]={0.0},F2np_err1[MAXBIN]={0.0};
     Double_t x_KP[7][110]={0.0},F2n_KP[7][110]={0.0},F2p_KP[7][110]={0.0};
     Double_t x_CJ[MAXBIN]={0.0},RD[MAXBIN]={0.0};
     Double_t x_CJ1[15]={0.0},F2n_CJ[15]={0.0},F2p_CJ[15]={0.0},F2nE_CJ[15]={0.0},F2pE_CJ[15]={0.0};
     Double_t x2[MAXBIN]={0.0},Ratio2[MAXBIN]={0.0},Rerr2[MAXBIN]={0.0},F2np2[MAXBIN]={0.0},F2np2_err[MAXBIN]={0.0};
     Double_t x_H3He[MAXBIN]={0.0},Ratio_H3He[MAXBIN]={0.0},Rerr_H3He[MAXBIN]={0.0},F2np_H3He[MAXBIN]={0.0},F2np_err_H3He[MAXBIN]={0.0};
     
     TString Rfile="newbin/Dp_final.dat";
     int nbin=ReadYield(Rfile,x,Ratio,Rerr);
     Rfile="newbin/H3He_final.dat";
     int nbin_H3He=ReadYield(Rfile,x_H3He,Ratio_H3He,Rerr_H3He);
     Rfile="Model/NP_CJ.out";
     int nbin1=ReadModel(Rfile,x1,Q2,F2np1,F2np_err1);
     Rfile="Model/RD_CJ.out";
     int nbin_CJ=ReadCJ(Rfile,x_CJ,RD);
     Rfile="Model/F2dp_Whitlow.out";
     int nbin_W=ReadModel(Rfile,x2,Q2,Ratio2,Rerr2);
     Rfile="Model/CJ_all.csv";
     int nbin_CJ1=ReadCJall(Rfile,x_CJ1,F2p_CJ,F2pE_CJ,F2n_CJ,F2nE_CJ);


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
     TGraphErrors *gH3He=new TGraphErrors();  // n/p from H3/He
     TGraphErrors *gDp_CJ=new TGraphErrors(); //N/P extract with CJ R
     TGraphErrors *gDp1=new TGraphErrors(); //CJ N/P
     TGraphErrors *gDp2=new TGraphErrors(); //NMC N/P
     TGraphErrors *gDp3=new TGraphErrors(); // N/P from Whitlow d/p with KP R
     TGraphErrors *gDp4=new TGraphErrors(); // N/P from Whitlow d/p with CJ R
     TGraph *gDp_KP[7]; //KP N/P

     for(int ii=0;ii<nbin;ii++){
	Double_t tmpR=D2_EMC(x[ii]);
	F2np[ii]=Ratio[ii]/tmpR-1.0;
	F2np_err[ii]=Rerr[ii]/tmpR;
	gDp->SetPoint(ii,x[ii],F2np[ii]);
	gDp->SetPointError(ii,0.0,F2np_err[ii]);
     }

     for(int ii=0;ii<nbin_H3He;ii++){
	if(x_H3He[ii]>0.4)continue;
	Double_t tmpR31=H3_EMC(x_H3He[ii]);
	Double_t tmpR32=He_EMC(x_H3He[ii]);
	Double_t SR=tmpR31/tmpR32;

        F2np_H3He[ii]=(2.0*Ratio_H3He[ii]-SR)/(2.0*SR-Ratio_H3He[ii]);
        F2np_err_H3He[ii]=abs(3.0*SR/((2.0*SR-Ratio_H3He[ii])*(2.0*SR-Ratio_H3He[ii])))*Rerr_H3He[ii];

	gH3He->SetPoint(ii,x_H3He[ii],F2np_H3He[ii]);
	gH3He->SetPointError(ii,0.0,F2np_err_H3He[ii]);
     }

     for(int ii=0;ii<nbin_W;ii++){
	Double_t tmpR=D2_EMC(x2[ii]);
	F2np2[ii]=Ratio2[ii]/tmpR-1.0;
	F2np2_err[ii]=Rerr2[ii]/tmpR;
	gDp3->SetPoint(ii,x2[ii],F2np2[ii]);
	gDp3->SetPointError(ii,0.0,F2np2_err[ii]);
     }

     for(int ii=0;ii<nbin_CJ;ii++){
	Double_t tmpNP=Ratio[ii]/RD[ii]-1.0;
	Double_t tmpNPerr=Rerr[ii]/RD[ii];
	gDp_CJ->SetPoint(ii,x[ii],tmpNP);
	gDp_CJ->SetPointError(ii,0.0,tmpNPerr);
     }

     for(int ii=0;ii<nbin_CJ;ii++){
	Double_t tmpNP=Ratio2[ii]/RD[ii]-1.0;
	Double_t tmpNPerr=Rerr2[ii]/RD[ii];
	gDp4->SetPoint(ii,x2[ii],tmpNP);
	gDp4->SetPointError(ii,0.0,tmpNPerr);
     }

    for(int ii=0;ii<nbin_CJ;ii++){
        gDp1->SetPoint(ii,x_CJ1[ii],F2n_CJ[ii]/F2p_CJ[ii]);
        Double_t tmpErr=F2n_CJ[ii]/F2p_CJ[ii]*sqrt(pow(F2pE_CJ[ii]/F2p_CJ[ii],2)+pow(F2nE_CJ[ii]/F2n_CJ[ii],2));
        gDp1->SetPointError(ii,0,tmpErr);
    }


     for(int ii=0;ii<nbin1;ii++){
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
     gDp->SetMarkerSize(1.6);
     gH3He->SetMarkerStyle(8);
     gH3He->SetMarkerColor(4);
     gH3He->SetMarkerSize(1.6);
     gDp_CJ->SetMarkerStyle(8);
     gDp_CJ->SetMarkerColor(4);
//     gDp3->SetMarkerStyle(8);
//     gDp3->SetMarkerColor(8);
//     gDp4->SetMarkerStyle(8);
//     gDp4->SetMarkerColor(6);
//     gDp4->SetMarkerSize(1.6);
     gDp3->SetFillStyle(3003);
     gDp3->SetFillColor(38);
     gDp4->SetFillStyle(3004);
     gDp4->SetFillColor(46);



     gDp1->SetFillStyle(3001);
     gDp1->SetFillColor(40);
     gDp2->SetLineStyle(8);
     gDp2->SetLineColor(2);
     gDp2->SetLineWidth(2);
//     f1->SetLineStyle(1);
//     f1->SetLineColor(4);
     mg->Add(gDp3,"E3");
     mg->Add(gDp4,"E3");
     mg->Add(gDp1,"E3");
     mg->Add(gDp,"P");
     mg->Add(gH3He,"P");
//     mg->Add(gDp_CJ,"P");
     mg->Add(gDp2,"L");
     mg->Add(gDp_KP[6],"L");
     mg->Draw("A");
     mg->SetTitle(";Bjorken x;F_{2}^{n}/F_{2}^{p}");
//     f1->Draw("same");
     mg->GetYaxis()->SetRangeUser(0.5,0.85);

   auto leg1=new TLegend(0.5,0.75,0.9,0.9);
   leg1->SetNColumns(2);
   leg1->AddEntry(gDp,"#scale[1]{MARATHON d/p (KP)}","P");
   leg1->AddEntry(gH3He,"#scale[1]{MARATHON H3/He (KP)}","P");
   leg1->AddEntry(gDp4,"#scale[1]{Whitlow d/p (CJ)}","F");
//   leg1->AddEntry(gDp_CJ,"#scale[1]{MARATHON (CJ)}","P");
   leg1->AddEntry(gDp1,"#scale[1]{CJ15}","F");
   leg1->AddEntry(gDp2,"#scale[1]{NMC}","L");
   leg1->AddEntry(gDp_KP[6],"#scale[1]{KP}","L");
   leg1->AddEntry(gDp3,"#scale[1]{Whitlow (KP)}","F");
   leg1->Draw();

     c1->Print("NP_Dp.pdf");
}
