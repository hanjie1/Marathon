#include "ReadFile.h"
#include "libEMCISO.h"
void He3_EMC_compare(){
     Double_t x[19]={0.0},Q2[19]={0.0},F2np[19]={0.0},F2np_err[19]={0.0},F2np_ST[19]={0.0},F2np_SY[19]={0.0};
     Double_t HeD[19]={0.0},HeD_err[19]={0.0},HeD_ST[19]={0.0},HeD_SY[19]={0.0};
     Double_t HeD_iso[19]={0.0},HeD_err_iso[19]={0.0},HeD_iso_ST[19]={0.0},HeD_iso_SY[19]={0.0};
     Double_t pHeD[19]={0.0},pHeD_err[19]={0.0};
     Double_t pHeD_norm[19]={0.0},pHeD_err_norm[19]={0.0};
     Double_t HeD_KP[19]={0.0},HeD_KPerr[19]={0.0},HeD_SLAC[19]={0.0},HeD_SLAC_err[19]={0.0};

     Double_t Z=2.0,A=3.0;
     Double_t Nc=0.024;

     TString Rfile;
     Rfile="F2n_F2p/results/F2np_final.dat";
     int nbin1=ReadNP(Rfile,x,F2np,F2np_err,F2np_ST,F2np_SY);
     Rfile="newbin/HeD_final.dat";
     int nbin2=ReadYield(Rfile,x,Q2,HeD,HeD_err,HeD_ST,HeD_SY);

     for(int ii=0;ii<19;ii++){
	 Double_t tmpIso=(1.0+F2np[ii])/(Z+(A-Z)*F2np[ii]);
  	 pHeD[ii]=HeD[ii]*2.0/A;	
  	 pHeD_err[ii]=HeD_err[ii]*2.0/A;	
  	 pHeD_norm[ii]=pHeD[ii]*(1.+Nc);	
  	 pHeD_err_norm[ii]=pHeD_err[ii]*(1.+Nc);	

	 HeD_iso[ii]=HeD[ii]*tmpIso*(1.0+Nc);
	 HeD_err_iso[ii]=HeD_iso[ii]*sqrt(pow(HeD_err[ii]/HeD[ii],2)
			+pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_err[ii],2));
	 HeD_iso_ST[ii]=HeD_iso[ii]*sqrt(pow(HeD_ST[ii]/HeD[ii],2)
                        +pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_ST[ii],2));
	 HeD_iso_SY[ii]=HeD_iso[ii]*sqrt(pow(HeD_SY[ii]/HeD[ii],2)
                        +pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_SY[ii],2));

	 Double_t eKP_EMC=He_ISO(x[ii]);
         HeD_KP[ii]=HeD_iso[ii]/eKP_EMC;   
	 HeD_KPerr[ii]=HeD_err_iso[ii]/eKP_EMC;

	 Double_t eSLAC_EMC=SLAC_EMC(x[ii]);
         HeD_SLAC[ii]=HeD_iso[ii]/eSLAC_EMC;    
         HeD_SLAC_err[ii]=HeD_err_iso[ii]/eSLAC_EMC;    
     } 

     TGraphErrors *gHeD_KP=new TGraphErrors(19,x,HeD_KP,0,HeD_KPerr);
     TGraphErrors *gHeD_SLAC=new TGraphErrors(19,x,HeD_SLAC,0,HeD_SLAC_err);

   gStyle->SetEndErrorSize(4);
     TCanvas *c1=new TCanvas("c1","c1",1500,1200);
     TMultiGraph *mg=new TMultiGraph();
     gHeD_KP->SetMarkerStyle(8);
     gHeD_KP->SetMarkerColor(4);
     gHeD_KP->SetMarkerSize(1.5);
     gHeD_KP->SetLineColor(4);
     gHeD_SLAC->SetMarkerStyle(8);
     gHeD_SLAC->SetMarkerColor(2);
     gHeD_SLAC->SetMarkerSize(1.5);
     gHeD_SLAC->SetLineColor(2);
     mg->Add(gHeD_SLAC,"P");
     mg->Add(gHeD_KP,"P");
     mg->Draw("A");
//     mg->SetTitle(";Bjorken x;(F_{2}^{^{3}He}/F_{2}^{^{2}H})_{iso}");
//     mg->GetYaxis()->SetRangeUser(0.85,1.25);
//     mg->GetXaxis()->SetRangeUser(0,0.9);

   auto leg1=new TLegend(0.15,0.75,0.7,0.9);
   leg1->SetNColumns(3);
   leg1->AddEntry(gHeD_KP,"KP","P");
   leg1->AddEntry(gHeD_SLAC,"SLAC","P");
   leg1->SetMargin(0.3);
   leg1->Draw();

}
