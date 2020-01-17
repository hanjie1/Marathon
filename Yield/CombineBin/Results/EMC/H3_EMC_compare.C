#include "ReadFile.h"
#include "libEMCISO.h"
void H3_EMC_compare(){
     Double_t x[19]={0.0},Q2[19]={0.0},F2np[19]={0.0},F2np_err[19]={0.0},F2np_ST[19]={0.0},F2np_SY[19]={0.0};
     Double_t H3D[19]={0.0},H3D_err[19]={0.0},H3D_ST[19]={0.0},H3D_SY[19]={0.0};
     Double_t H3D_iso[19]={0.0},H3D_err_iso[19]={0.0},H3D_iso_ST[19]={0.0},H3D_iso_SY[19]={0.0};
     Double_t pH3D[19]={0.0},pH3D_err[19]={0.0};
     Double_t pH3D_norm[19]={0.0},pH3D_err_norm[19]={0.0};
     Double_t H3D_KP[19]={0.0},H3D_KPerr[19]={0.0},H3D_SLAC[19]={0.0},H3D_SLAC_err[19]={0.0};

     Double_t Z=1.0,A=3.0;
     Double_t Nc=0.0;

     TString Rfile;
     Rfile="F2n_F2p/results/F2np_final.dat";
     int nbin1=ReadNP(Rfile,x,F2np,F2np_err,F2np_ST,F2np_SY);
     Rfile="newbin/H3D_final.dat";
     int nbin2=ReadYield(Rfile,x,Q2,H3D,H3D_err,H3D_ST,H3D_SY);

     for(int ii=0;ii<19;ii++){
	 Double_t tmpIso=(1.0+F2np[ii])/(Z+(A-Z)*F2np[ii]);
  	 pH3D[ii]=H3D[ii]*2.0/A;	
  	 pH3D_err[ii]=H3D_err[ii]*2.0/A;	
  	 pH3D_norm[ii]=pH3D[ii]*(1.+Nc);	
  	 pH3D_err_norm[ii]=pH3D_err[ii]*(1.+Nc);	

	 H3D_iso[ii]=H3D[ii]*tmpIso*(1.0+Nc);
	 H3D_err_iso[ii]=H3D_iso[ii]*sqrt(pow(H3D_err[ii]/H3D[ii],2)
			+pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_err[ii],2));
	 H3D_iso_ST[ii]=H3D_iso[ii]*sqrt(pow(H3D_ST[ii]/H3D[ii],2)
                        +pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_ST[ii],2));
	 H3D_iso_SY[ii]=H3D_iso[ii]*sqrt(pow(H3D_SY[ii]/H3D[ii],2)
                        +pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_SY[ii],2));

	 Double_t eKP_EMC=He_ISO(x[ii]);
         H3D_KP[ii]=H3D_iso[ii]/eKP_EMC;   
	 H3D_KPerr[ii]=H3D_err_iso[ii]/eKP_EMC;

	 Double_t eSLAC_EMC=SLAC_EMC(x[ii]);
         H3D_SLAC[ii]=H3D_iso[ii]/eSLAC_EMC;    
         H3D_SLAC_err[ii]=H3D_err_iso[ii]/eSLAC_EMC;    
     } 

     TGraphErrors *gH3D_KP=new TGraphErrors(19,x,H3D_KP,0,H3D_KPerr);
     TGraphErrors *gH3D_SLAC=new TGraphErrors(19,x,H3D_SLAC,0,H3D_SLAC_err);

   gStyle->SetEndErrorSize(4);
     TCanvas *c1=new TCanvas("c1","c1",1500,1200);
     TMultiGraph *mg=new TMultiGraph();
     gH3D_KP->SetMarkerStyle(8);
     gH3D_KP->SetMarkerColor(4);
     gH3D_KP->SetMarkerSize(1.5);
     gH3D_KP->SetLineColor(4);
     gH3D_SLAC->SetMarkerStyle(8);
     gH3D_SLAC->SetMarkerColor(2);
     gH3D_SLAC->SetMarkerSize(1.5);
     gH3D_SLAC->SetLineColor(2);
     mg->Add(gH3D_SLAC,"P");
     mg->Add(gH3D_KP,"P");
     mg->Draw("A");
//     mg->SetTitle(";Bjorken x;(F_{2}^{^{3}He}/F_{2}^{^{2}H})_{iso}");
//     mg->GetYaxis()->SetRangeUser(0.85,1.25);
//     mg->GetXaxis()->SetRangeUser(0,0.9);

   auto leg1=new TLegend(0.15,0.75,0.7,0.9);
   leg1->SetNColumns(3);
   leg1->AddEntry(gH3D_KP,"KP","P");
   leg1->AddEntry(gH3D_SLAC,"SLAC","P");
   leg1->SetMargin(0.3);
   leg1->Draw();

}
