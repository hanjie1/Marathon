#include "ReadFile.h"
#include "libEMCR.h"

int H3He_np(){
   Double_t x[19]={0.0},Q2[19]={0.0};
   Double_t H3He[19]={0.0},H3He_err[19]={0.0},H3He_ST[19]={0.0},H3He_SY[19]={0.0};
   Double_t x_KP[110]={0.0},F2n_KP[110]={0.0},F2p_KP[110]={0.0};

   TString Rfile="newbin/H3He_final.dat";
   int nbin=ReadYieldFinal(Rfile,x,Q2,H3He,H3He_ST,H3He_SY);
   Rfile="Model/F2dis_os1tm1ht1mec1_Dav18_He3Salme";
   int nbin_KP=ReadKP(Rfile,x_KP,F2p_KP,F2n_KP);

   Double_t Nc=0.976;

    ofstream outfile;
    outfile.open("results/F2np_final.dat");
    outfile<<"x   Q2   n/p    total_err     stat_err     sys_err"<<endl;;

    ofstream outfile1;
    outfile1.open("ForThesis.dat");

    Double_t H3He_np[19]={0.0},H3He_npErr[19]={0.0};
    Double_t np_ST[19]={0.0},np_SY[19]={0.0};
    for(int ii=0;ii<19;ii++){
	Double_t He_R=He_EMC(x[ii]);
	Double_t H3_R=H3_EMC(x[ii]);
	Double_t SR=H3_R/He_R;

	Double_t tmpH3He=H3He[ii]*Nc;
	Double_t tmpST=H3He_ST[ii]*Nc;
	Double_t tmpSY=H3He_SY[ii]*Nc;
	Double_t tmpH3He_err=sqrt(tmpST*tmpST+tmpSY*tmpSY);
 	
	H3He_np[ii]=(2.0*tmpH3He-SR)/(2.0*SR-tmpH3He);
	H3He_npErr[ii]=abs(3.0*SR/((2.0*SR-tmpH3He)*(2.0*SR-tmpH3He)))*tmpH3He_err;
	np_ST[ii]=abs(3.0*SR/((2.0*SR-tmpH3He)*(2.0*SR-tmpH3He)))*tmpST;
	np_SY[ii]=abs(3.0*SR/((2.0*SR-tmpH3He)*(2.0*SR-tmpH3He)))*tmpSY;
	outfile<<x[ii]<<"  "<<Q2[ii]<<"  "<<H3He_np[ii]<<"  "<<H3He_npErr[ii]<<"  "<<np_ST[ii]<<"  "<<np_SY[ii]<<endl;
	outfile1<<fixed<<setprecision(2)<<x[ii]<<" & "<<Q2[ii]<<" & "<<setprecision(3)<<H3He_np[ii]<<" & "<<np_ST[ii]<<" &  "<<np_SY[ii]<<"  \\\\"<<endl;
	outfile1<<"\\hline"<<endl;
    }
    outfile.close();
    outfile1.close();
    TGraphErrors *gH3He=new TGraphErrors(19,x,H3He_np,0,H3He_npErr);

    TGraphErrors *gKP=new TGraphErrors();
    for(int ii=0;ii<nbin_KP;ii++){
	if(x_KP[ii]>0.85)continue;
	gKP->SetPoint(ii,x_KP[ii],F2n_KP[ii]/F2p_KP[ii]);
    }

   gStyle->SetEndErrorSize(4);

    TCanvas *c1=new TCanvas("c1","c1",1500,1200);
    TMultiGraph *mg1=new TMultiGraph();
    gH3He->SetMarkerColor(2);
    gH3He->SetMarkerStyle(8);
    gH3He->SetMarkerSize(2);
    gH3He->SetLineColor(2);
    gKP->SetLineColor(4);
    gKP->SetLineStyle(1);
    gKP->SetLineWidth(2);
    mg1->Add(gH3He,"P");
    mg1->Add(gKP,"L");
    mg1->Draw("A"); 
    mg1->SetTitle(";Bjorken x;F_{2}^{n} / F_{2}^{p}");
    mg1->GetYaxis()->SetLabelOffset(0.0005);
    mg1->GetYaxis()->SetRangeUser(0.3,0.85);
 
   auto leg1=new TLegend(0.65,0.7,0.85,0.85);
   leg1->AddEntry(gH3He,"#scale[1]{MARATHON}","P");
   leg1->AddEntry(gKP,"#scale[1]{KP model}","L");
   leg1->SetMargin(0.3);
   leg1->Draw();

   c1->Print("F2np_final.pdf");

    return 0;
}
