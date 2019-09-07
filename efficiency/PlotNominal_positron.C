#include "SetCut.h"
#include "TF1.h"

Double_t background(Double_t *x, Double_t *par){
	return par[0]+par[1]*TMath::Exp(par[2]*x[0]);
}

Double_t peak(Double_t *x, Double_t *par){
	return par[2]*TMath::Gaus(x[0],par[0],par[1]);
}

Double_t Fitfunction(Double_t *x, Double_t *par){
        return background(x,par) + peak(x,&par[3]);
}

void PlotNominal_positron()
{
     TCut CK = "L.cer.asum_c>1500";
     TCut VZ = Form("L.tr.vz>%f && L.tr.vz<%f",vz_min[1],vz_max[1]);
  
     TString rootpath="/lustre19/expphy/cache/halla/triton/prod/marathon/pass2/kin1p";
     TChain *T=new TChain("T");
     T->Add(Form("%s/tritium_2437*",rootpath.Data()));

     TH1F *hEP=new TH1F("hEP","positron E/P distribution",100,0.,1.3);

     T->Draw("(L.prl1.e+L.prl2.e)/(1000.0*L.gold.p)>>hEP",trigger2+CK+TRK+ACC+beta+VZ);

     TF1 *fitFuc=new TF1("fitFuc",Fitfunction,0,1.3,6);

     fitFuc->SetParameters(1,1,1,1,1,1);

     hEP->Fit("fitFuc");


}


