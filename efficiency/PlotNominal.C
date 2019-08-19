#include "SetCut.h"
void PlotNominal()
{
     TCut CK = "L.cer.asum_c>1500";
     TCut VZ = Form("L.tr.vz>%f && L.tr.vz<%f",vz_min[1],vz_max[1]);
  
     TString rootpath="/lustre19/expphy/cache/halla/triton/prod/marathon/pass2/kin1";
     TChain *T=new TChain("T");
     T->Add(Form("%s/tritium_1213*",rootpath.Data()));

     TChain *T1=new TChain("T");
     T1->Add(Form("%s/tritium_1209*",rootpath.Data()));

     TH1F *htar_H3=new TH1F("htar_H3","target vertex z",400,-0.2,0.2);
     TH1F *htar_EM=new TH1F("htar_EM","target vertex z",400,-0.2,0.2);

     T->Draw("L.tr.vz>>htar_H3",trigger2+CK+Ep+TRK+ACC+beta);
     T1->Draw("L.tr.vz>>htar_EM",trigger2+CK+Ep+TRK+ACC+beta,"same");

     htar_H3->GetXaxis()->SetRangeUser(-0.15,-0.1);
     int max_bin = htar_H3->GetMaximumBin();
     htar_H3->GetXaxis()->SetRangeUser(-0.2,0.2);
     Double_t ptar_H3 = htar_H3->GetBinContent(max_bin);

     htar_EM->GetXaxis()->SetRangeUser(-0.15,-0.1);
     max_bin = htar_EM->GetMaximumBin();
     htar_EM->GetXaxis()->SetRangeUser(-0.2,0.2);
     Double_t ptar_EM = htar_EM->GetBinContent(max_bin);

     htar_H3->Scale(ptar_EM/ptar_H3);
     htar_EM->SetLineColor(2);

     auto leg1=new TLegend(0.6,0.6,0.85,0.85);
     leg1->AddEntry(htar_H3,"{}^3H target","1L");
     leg1->AddEntry(htar_EM,"empty cell","L");
     //leg1->Draw();

}


