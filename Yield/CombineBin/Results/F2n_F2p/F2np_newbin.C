#include "ReadFile.h"
#include "libEMCR.h"

int F2np_newbin(){
   Double_t x[19]={0.0},H3D[19]={0.0},H3D_err[19]={0.0};
   Double_t HeD[19]={0.0},HeD_err[19]={0.0};
   Double_t H3He[19]={0.0},H3He_err[19]={0.0};
   Double_t Dp_x[8]={0.0},Dp[8]={0.0},Dp_err[8]={0.0};
   Double_t x_KP[110]={0.0},F2n_KP[110]={0.0},F2p_KP[110]={0.0};
   Double_t x_CJ[15]={0.0},F2n_CJ[15]={0.0},F2p_CJ[15]={0.0},F2nE_CJ[15]={0.0},F2pE_CJ[15]={0.0};

   TString Rfile="newbin/Dp_final.dat";
   int nbin1=ReadYield(Rfile,Dp_x,Dp,Dp_err);
   Rfile="newbin/HeD_final.dat";
   int nbin2=ReadYield(Rfile,x,HeD,HeD_err);
   Rfile="newbin/H3D_final.dat";
   nbin2=ReadYield(Rfile,x,H3D,H3D_err);
   Rfile="newbin/H3He_final.dat";
   nbin2=ReadYield(Rfile,x,H3He,H3He_err);

   Rfile="Model/F2dis_os1tm1ht1mec1_Dav18_He3Salme";
   int nbin_KP=ReadKP(Rfile,x_KP,F2p_KP,F2n_KP);
   Rfile="Model/CJ_all.csv";
   int nbin_CJ=ReadCJall(Rfile,x_CJ,F2p_CJ,F2pE_CJ,F2n_CJ,F2nE_CJ);

    Double_t H3He_np[19]={0.0},H3He_npErr[19]={0.0};
    Double_t H3D_np[19]={0.0},H3D_npErr[19]={0.0};
    Double_t HeD_np[19]={0.0},HeD_npErr[19]={0.0};
    ofstream outfile1;
    outfile1.open("results/F2np.dat");
    outfile1<<"x    H3He   H3He_err    HeD    HeD_err    H3D    H3D_err "<<endl;
    for(int ii=0;ii<19;ii++){
	Double_t He_R=He_EMC(x[ii]);
	Double_t H3_R=H3_EMC(x[ii]);
	Double_t D2_R=D2_EMC(x[ii]);

	Double_t SR=H3_R/He_R;
	Double_t SR_HeD=He_R/D2_R;
	Double_t SR_H3D=H3_R/D2_R;

	H3He_np[ii]=(2.0*H3He[ii]-SR)/(2.0*SR-H3He[ii]);
	H3He_npErr[ii]=abs(3.0*SR/((2.0*SR-H3He[ii])*(2.0*SR-H3He[ii])))*H3He_err[ii];

	HeD_np[ii]=(HeD[ii]-2.0*SR_HeD)/(SR_HeD-HeD[ii]);
	HeD_npErr[ii]=abs(SR_HeD/((SR_HeD-HeD[ii])*(SR_HeD-HeD[ii])))*HeD_err[ii];

	H3D_np[ii]=(H3D[ii]-SR_H3D)/(2.0*SR_H3D-H3D[ii]);
	H3D_npErr[ii]=abs(SR_H3D/((2.0*SR_H3D-H3D[ii])*(2.0*SR_H3D-H3D[ii])))*H3D_err[ii];
	outfile1<<x[ii]<<"  "<<H3He_np[ii]<<"  "<<H3He_npErr[ii]<<"  "<<HeD_np[ii]<<"  "<<HeD_npErr[ii]<<"  "
		<<H3D_np[ii]<<"  "<<H3D_npErr[ii]<<endl;
    }
    outfile1.close();

    ofstream outfile2;
    outfile2.open("results/F2np_Dp.dat");
    outfile2<<"x    Dp    Dp_err "<<endl;

    Double_t Dp_np[8]={0.0},Dp_npErr[8]={0.0},D2_R[8]={0.0};
    for(int ii=0;ii<8;ii++){
	D2_R[ii]=D2_EMC(Dp_x[ii]);
	Dp_np[ii]=Dp[ii]/D2_R[ii]-1.0;
	Dp_npErr[ii]=1.0/D2_R[ii]*Dp_err[ii];
	outfile2<<Dp_x[ii]<<"  "<<Dp_np[ii]<<"  "<<Dp_npErr[ii]<<endl;
    } 
    outfile2.close();

    TGraphErrors *gH3He=new TGraphErrors(19,x,H3He_np,0,H3He_npErr);
    TGraphErrors *gH3D=new TGraphErrors(19,x,H3D_np,0,H3D_npErr);
    TGraphErrors *gHeD=new TGraphErrors(19,x,HeD_np,0,HeD_npErr);
    TGraphErrors *gDp=new TGraphErrors(8,Dp_x,Dp_np,0,Dp_npErr);

    TGraphErrors *gKP=new TGraphErrors();
    for(int ii=0;ii<nbin_KP;ii++){
	if(x_KP[ii]>0.85)continue;
	gKP->SetPoint(ii,x_KP[ii],F2n_KP[ii]/F2p_KP[ii]);
    }

    TGraphErrors *gCJ=new TGraphErrors();
    for(int ii=0;ii<nbin_CJ;ii++){
	gCJ->SetPoint(ii,x_CJ[ii],F2n_CJ[ii]/F2p_CJ[ii]);
	Double_t tmpErr=F2n_CJ[ii]/F2p_CJ[ii]*sqrt(pow(F2pE_CJ[ii]/F2p_CJ[ii],2)+pow(F2nE_CJ[ii]/F2n_CJ[ii],2));
	gCJ->SetPointError(ii,0,tmpErr);
    }

   TGraphErrors *Dp_norm=new TGraphErrors(1);
   Dp_norm->SetPoint(0,Dp_x[3],0.4);
   Dp_norm->SetPointError(0,0,0.0079*Dp[3]/D2_R[3]);

   gStyle->SetEndErrorSize(4);

    TCanvas *c1=new TCanvas("c1","c1",1500,1200);
    TMultiGraph *mg1=new TMultiGraph();
    gDp->SetMarkerColor(4);
    gDp->SetMarkerStyle(24);
    gDp->SetMarkerSize(2);
    gDp->SetLineColor(4);
    Dp_norm->SetLineColor(4);
    Dp_norm->SetLineWidth(2);
    
    gH3He->SetMarkerColor(2);
    gH3He->SetMarkerStyle(26);
    gH3He->SetMarkerSize(2);
    gH3He->SetLineColor(2);
    gHeD->SetMarkerColor(kViolet-1);
    gHeD->SetMarkerStyle(32);
    gHeD->SetMarkerSize(2);
    gHeD->SetLineColor(kViolet-1);
    gH3D->SetMarkerColor(1);
    gH3D->SetMarkerStyle(25);
    gH3D->SetMarkerSize(2);
    gH3D->SetLineColor(1);
    gKP->SetLineColor(8);
    gKP->SetLineStyle(1);
    gKP->SetLineWidth(2);
    gCJ->SetFillColor(29);
    gCJ->SetFillStyle(3001);
    mg1->Add(gCJ,"E3");
    mg1->Add(gDp,"P");
    mg1->Add(gH3He,"P");
    mg1->Add(gKP,"L");
    mg1->Add(Dp_norm,"L");
//    mg1->Add(gH3D);
//    mg1->Add(gHeD);
    mg1->Draw("A"); 
    mg1->SetTitle(";Bjorken x;F_{2}^{n} / F_{2}^{p}");
    mg1->GetYaxis()->SetLabelOffset(0.0005);
    mg1->GetYaxis()->SetRangeUser(0.3,0.85);
 
   auto leg1=new TLegend(0.65,0.65,0.85,0.85);
   //leg1->SetNColumns(2);
   leg1->AddEntry(gDp,"#scale[0.8]{F_{2}^{^{2}H} / F_{2}^{^{1}H}}","P");
   leg1->AddEntry(gH3He,"#scale[0.8]{F_{2}^{^{3}H} / F_{2}^{^{3}He}}","P");
   leg1->AddEntry(gKP,"#scale[0.7]{KP model}","L");
   leg1->AddEntry(gCJ,"#scale[0.7]{CJ model}","F");
   leg1->SetMargin(0.4);
   leg1->Draw();

   TLatex latex;
   latex.SetTextSize(0.025);
   latex.DrawLatex(Dp_x[3]+0.04,0.4,"Norm. uncer. from F_{2}^{^{2}H} / F_{2}^{^{1}H}");

   c1->Print("F2np_H3He.pdf");

    TCanvas *c2=new TCanvas("c2","c2",1500,1200);
    TMultiGraph *mg2=new TMultiGraph();
    mg2->Add(gDp,"P");
    mg2->Add(gH3He,"P");
//    mg2->Add(gKP,"L");
    mg2->Add(gH3D,"P");
    mg2->Add(gHeD,"P");
    mg2->Draw("A"); 
    mg2->SetTitle(";Bjorken x;F_{2}^{n} / F_{2}^{p}");
    mg2->GetYaxis()->SetLabelOffset(0.0005);
    mg2->GetXaxis()->SetRangeUser(0.2,0.4);
    mg2->GetYaxis()->SetRangeUser(0.4,0.95);
 
   auto leg2=new TLegend(0.7,0.7,0.9,0.9);
   //leg2->SetNColumns(2);
   leg2->AddEntry(gDp,"#scale[0.8]{F_{2}^{^{2}H} / F_{2}^{^{1}H}}","P");
   leg2->AddEntry(gH3He,"#scale[0.8]{F_{2}^{^{3}H} / F_{2}^{^{3}He}}","P");
   leg2->AddEntry(gH3D,"#scale[0.8]{F_{2}^{^{3}H} / F_{2}^{^{2}H}}","P");
   leg2->AddEntry(gHeD,"#scale[0.8]{F_{2}^{^{3}He} / F_{2}^{^{2}H}}","P");
//   leg2->AddEntry(gKP,"#scale[0.7]{KP model}","L");
//   leg2->AddEntry(gCJ,"#scale[0.7]{CJ model}","F");
   leg2->SetMargin(0.4);
   leg2->Draw();
   c2->Print("F2np_all.pdf");
 
    return 0;
}
