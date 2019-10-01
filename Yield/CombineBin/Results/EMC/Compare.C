#include "ReadFile.h"
void Compare(){
     Double_t x[19]={0.0},H3D[19]={0.0},H3D_err[19]={0.0},H3D_ST[19]={0.0},H3D_SY[19]={0.0};
     Double_t HeD[19]={0.0},HeD_err[19]={0.0},HeD_ST[19]={0.0},HeD_SY[19]={0.0};

     TString Rfile;
     Rfile="He_EMC_iso.dat";
     ReadISO(Rfile,x,HeD,HeD_ST,HeD_SY);
     Rfile="H3_EMC_iso.dat";
     ReadISO(Rfile,x,H3D,H3D_ST,H3D_SY);

     TGraphErrors *gHeD=new TGraphErrors();
     TGraphErrors *gH3D=new TGraphErrors();

     for(int ii=0;ii<19;ii++){
	gHeD->SetPoint(ii,x[ii],HeD[ii]);
	HeD_err[ii]=sqrt(HeD_ST[ii]*HeD_ST[ii]+HeD_SY[ii]*HeD_SY[ii]);
	gHeD->SetPointError(ii,0,HeD_err[ii]);
	
	gH3D->SetPoint(ii,x[ii],H3D[ii]);
	H3D_err[ii]=sqrt(H3D_ST[ii]*H3D_ST[ii]+H3D_SY[ii]*H3D_SY[ii]);
	gH3D->SetPointError(ii,0,H3D_err[ii]);
     }

     TCanvas *c1=new TCanvas("c1","c1",1500,1200);
     TMultiGraph *mg=new TMultiGraph();
     gH3D->SetMarkerStyle(8);
     gH3D->SetMarkerColor(2);
     gH3D->SetMarkerSize(2);
     gH3D->SetLineColor(2);
     gHeD->SetMarkerStyle(21);
     gHeD->SetMarkerColor(4);
     gHeD->SetMarkerSize(2);
     gHeD->SetLineColor(4);
     mg->Add(gH3D);
     mg->Add(gHeD);
     mg->Draw("AP");
     
     mg->SetTitle(";Bjorken x;");
     mg->GetYaxis()->SetRangeUser(0.90,1.1);
     mg->GetXaxis()->SetLimits(0.05,0.9);

     TLine *l1=new TLine(0.05,1,0.9,1);
     l1->SetLineColor(1);
     l1->SetLineStyle(7);
     l1->Draw("same");

   auto leg1=new TLegend(0.65,0.7,0.89,0.88);
   leg1->AddEntry(gHeD,"#scale[0.6]{(F_{2}^{^{3}He}/F_{2}^{^{2}H})_{iso}}","P");
   leg1->AddEntry(gH3D,"#scale[0.6]{(F_{2}^{^{3}H}/F_{2}^{^{2}H})_{iso}}","P");
   leg1->SetMargin(0.3);
   leg1->Draw();
   c1->Print("EMC_compare.pdf");

}
