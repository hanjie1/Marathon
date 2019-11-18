#include "ReadFile.h"
#include "libEMCISO.h"
void H3_EMC(){
     Double_t x[19]={0.0},Q2[19]={0.0},F2np[19]={0.0},F2np_err[19]={0.0},F2np_ST[19]={0.0},F2np_SY[19]={0.0};
     Double_t H3D[19]={0.0},H3D_err[19]={0.0},H3D_ST[19]={0.0},H3D_SY[19]={0.0};
     Double_t H3D_iso[19]={0.0},H3D_err_iso[19]={0.0},H3D_iso_ST[19]={0.0},H3D_iso_SY[19]={0.0};
     Double_t pH3D[19]={0.0},pH3D_err[19]={0.0};

     Double_t Z=1.0,A=3.0;

     TString Rfile;
     Rfile="F2n_F2p/results/F2np_final.dat";
     int nbin1=ReadNP(Rfile,x,F2np,F2np_err,F2np_ST,F2np_SY);
     Rfile="newbin/H3D_final.dat";
     int nbin2=ReadYield(Rfile,x,Q2,H3D,H3D_err,H3D_ST,H3D_SY);

     ofstream outfile;
     outfile.open("H3_EMC_iso.dat");
     for(int ii=0;ii<19;ii++){
	 Double_t tmpIso=(1.0+F2np[ii])/(Z+(A-Z)*F2np[ii]);
  	 pH3D[ii]=H3D[ii]*2.0/A;	
  	 pH3D_err[ii]=H3D_err[ii]*2.0/A;	

	 H3D_iso[ii]=H3D[ii]*tmpIso;
	 H3D_err_iso[ii]=H3D_iso[ii]*sqrt(pow(H3D_err[ii]/H3D[ii],2)
			+pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_err[ii],2));
	 H3D_iso_ST[ii]=H3D_iso[ii]*sqrt(pow(H3D_ST[ii]/H3D[ii],2)
                        +pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_ST[ii],2));
	 H3D_iso_SY[ii]=H3D_iso[ii]*sqrt(pow(H3D_SY[ii]/H3D[ii],2)
                        +pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_SY[ii],2));
	 outfile<<fixed<<setprecision(2)<<x[ii]<<" & "<<Q2[ii]<<" & "<<setprecision(3)<<H3D_iso[ii]<<" & "<<H3D_iso_ST[ii]<<" & "<<H3D_iso_SY[ii]<<"  \\\\"<<endl;
	 outfile<<"\\hline"<<endl; 
     } 
     outfile.close();

     TGraphErrors *gH3D=new TGraphErrors(19,x,pH3D,0,pH3D_err);
     TGraphErrors *gH3D_iso=new TGraphErrors(19,x,H3D_iso,0,H3D_err_iso);

     auto f1_KP=new TF1("f1_KP","He_ISO(x)",0.16,0.85);
     auto f1_SLAC=new TF1("f1_SLAC","SLAC_EMC(x)",0.16,0.85);
     auto f1_SLAC_den=new TF1("f1_SLAC_den","SLAC_EMC_Den(x,3.0,1.0)",0.16,0.85);


     TLine *l1=new TLine(0.05,1,0.9,1);
     l1->SetLineColor(1);
     l1->SetLineStyle(7);

   gStyle->SetEndErrorSize(4);
     TCanvas *c1=new TCanvas("c1","c1",1500,1200);
     TMultiGraph *mg=new TMultiGraph();
     gH3D->SetMarkerStyle(8);
     gH3D->SetMarkerColor(4);
     gH3D->SetMarkerSize(1.5);
     gH3D->SetLineColor(4);
     gH3D_iso->SetMarkerStyle(8);
     gH3D_iso->SetMarkerColor(2);
     gH3D_iso->SetMarkerSize(1.5);
     gH3D_iso->SetLineColor(2);
     mg->Add(gH3D,"P");
     mg->Add(gH3D_iso,"P");
     mg->Draw("A");
     mg->SetTitle(";Bjorken x;(F_{2}^{^{3}H}/F_{2}^{^{2}H})_{iso}");
     mg->GetXaxis()->SetLimits(0.05,0.9);
     mg->GetYaxis()->SetRangeUser(0.8,1.15);

     f1_KP->SetLineColor(12);
     f1_KP->SetLineStyle(1);
     f1_KP->Draw("same");

     f1_SLAC->SetLineColor(12);
     f1_SLAC->SetLineStyle(3);
     f1_SLAC->Draw("same");

     f1_SLAC_den->SetLineColor(12);
     f1_SLAC_den->SetLineStyle(7);
     f1_SLAC_den->Draw("same");

     l1->Draw("same");
     
   auto leg1=new TLegend(0.55,0.73,0.88,0.88);
   leg1->SetNColumns(2);
   leg1->AddEntry(gH3D_iso,"#scale[1]{MARATHON}","P");
   leg1->AddEntry(gH3D,"#scale[1]{MARATHON no Iso. Cor.}","P");
   leg1->AddEntry(f1_KP,"#scale[1]{KP model}","L");
   leg1->AddEntry(f1_SLAC,"#scale[1]{SLAC}","L");
   leg1->AddEntry(f1_SLAC_den,"#scale[1]{SLAC density model}","L");
   //leg1->SetMargin(0.4);
   leg1->Draw();

     c1->Print("H3_EMC.pdf");
}
