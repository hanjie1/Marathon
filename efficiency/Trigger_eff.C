#include "GetTrees.h"
#include "SetCut.h"
struct EFF_VAR{
      Double_t T1_eff;
      Double_t T1_err;
      Double_t T2_eff;
      Double_t T2_err;
};

void trigger(int run_number, int kin, EFF_VAR &TT){
     TString TreeName="T";
     TChain *T=GetTree(run_number,kin,TreeName);

     Double_t NT1=0.0,NT2=0.0,NT3=0.0,NT2_2=0.0;

     if(kin<16){
        TCut EE="(L.prl1.e+L.prl2.e)/3100.0>0.8 && (L.prl1.e+L.prl2.e)/3100.0<1.2";
        TCut PS="(L.prl1.e/L.prl2.e>0.5) && (L.prl1.e/L.prl2.e<4)";
        TCut CK = "L.cer.asum_c>2200";

        NT3=T->GetEntries(EE+CK+trigger3);
        NT2=T->GetEntries(EE+CK+trigger2);

        NT1=T->GetEntries(EE+trigger1+beta+ACC+PS);
        NT2_2=T->GetEntries(EE+trigger2+beta+ACC+PS);
     }

     if(kin==16){
        TCut EE="(R.ps.e+R.sh.e)/2900.0>0.8 && (R.ps.e+R.sh.e)/2900.0<1.2";
        TCut CK_R= "R.cer.asum_c>3000";
        TCut PS="(R.ps.e/R.sh.e>0.2) && (R.ps.e/R.sh.e<0.8)";

        NT3=T->GetEntries(EE+CK_R+trigger6+PS);
        NT2=T->GetEntries(EE+CK_R+trigger5+PS);

        NT1=T->GetEntries(EE+trigger4+beta_R+ACC_R+PS);
        NT2_2=T->GetEntries(EE+trigger5+beta_R+ACC_R+PS);
     }

     Double_t eff1=NT2/NT3;
     Double_t eff1_err=sqrt(1.0/NT2-1.0/NT3)*eff1;

     Double_t eff2=NT2_2/NT1;
     Double_t eff2_err=sqrt(1.0/NT2_2-1.0/NT1)*eff2;

     TT.T1_eff=eff1;
     TT.T1_err=eff1_err;

     TT.T2_eff=eff2;
     TT.T2_err=eff2_err;

     cout<<kin<<": "<<eff1<<"  "<<eff1_err<<"  "<<eff2<<"  "<<eff2_err<<endl;

     return;
}

void Trigger_eff(){
     int kin[12]={0,1,2,3,4,5,7,9,11,13,15,16};
 
     int H1_nrun[5]={2479,1213,1240,1286,2582};
     int D2_nrun[12]={2493,1214,1242,1288,2579,1350,1371,1405,1493,1639,1862,91280};
     int He_nrun[12]={2531,1219,1236,1283,2576,1353,1373,1408,1485,1632,1868,91365};
     int H3_nrun[12]={2503,1224,1245,1290,2577,1355,1377,1413,1489,1635,1880,91369};
     
     TGraphErrors *gH1_T1=new TGraphErrors(5);
     TGraphErrors *gD2_T1=new TGraphErrors(12);
     TGraphErrors *gHe_T1=new TGraphErrors(12);
     TGraphErrors *gH3_T1=new TGraphErrors(12);

     TGraphErrors *gH1_T2=new TGraphErrors(5);
     TGraphErrors *gD2_T2=new TGraphErrors(12);
     TGraphErrors *gHe_T2=new TGraphErrors(12);
     TGraphErrors *gH3_T2=new TGraphErrors(12);

     EFF_VAR H1_EFF[5];
     EFF_VAR D2_EFF[12];
     EFF_VAR He_EFF[12];
     EFF_VAR H3_EFF[12];

     for(int ii=0;ii<12;ii++){
         Double_t KKin=1.0*kin[ii];
         if(ii<5){
	    EFF_VAR H1_single;
            trigger(H1_nrun[ii],kin[ii],H1_single);
	    H1_EFF[ii]=H1_single;
            gH1_T1->SetPoint(ii,KKin,H1_single.T1_eff);
            gH1_T1->SetPointError(ii,0.0,H1_single.T1_err);
            gH1_T2->SetPoint(ii,KKin,H1_single.T2_eff);
            gH1_T2->SetPointError(ii,0.0,H1_single.T2_err);
         }

	 EFF_VAR D2_single;
         trigger(D2_nrun[ii],kin[ii],D2_single);
	 D2_EFF[ii]=D2_single;
         gD2_T1->SetPoint(ii,KKin,D2_single.T1_eff);
         gD2_T1->SetPointError(ii,0.0,D2_single.T1_err);
         gD2_T2->SetPoint(ii,KKin,D2_single.T2_eff);
         gD2_T2->SetPointError(ii,0.0,D2_single.T2_err);

	 EFF_VAR He_single;
         trigger(He_nrun[ii],kin[ii],He_single);
	 He_EFF[ii]=He_single;
         gHe_T1->SetPoint(ii,KKin,He_single.T1_eff);
         gHe_T1->SetPointError(ii,0.0,He_single.T1_err);
         gHe_T2->SetPoint(ii,KKin,He_single.T2_eff);
         gHe_T2->SetPointError(ii,0.0,He_single.T2_err);

	 EFF_VAR H3_single;
         trigger(H3_nrun[ii],kin[ii],H3_single);
	 H3_EFF[ii]=H3_single;
         gH3_T1->SetPoint(ii,KKin,H3_single.T1_eff);
         gH3_T1->SetPointError(ii,0.0,H3_single.T1_err);
         gH3_T2->SetPoint(ii,KKin,H3_single.T2_eff);
         gH3_T2->SetPointError(ii,0.0,H3_single.T2_err);
     }

        gH1_T1->SetMarkerStyle(8);
        gH1_T1->SetMarkerColor(1);
        gD2_T1->SetMarkerStyle(8);
        gD2_T1->SetMarkerColor(3);
        gHe_T1->SetMarkerStyle(8);
        gHe_T1->SetMarkerColor(4);
        gH3_T1->SetMarkerStyle(8);
        gH3_T1->SetMarkerColor(2);

        TCanvas *c1=new TCanvas("c1");
        TMultiGraph *mg=new TMultiGraph();
        mg->Add(gH1_T1);
        mg->Add(gD2_T1);
        mg->Add(gHe_T1);
        mg->Add(gH3_T1);
        mg->Draw("AP");
        mg->GetXaxis()->SetTitle("kinematic");
        mg->GetYaxis()->SetTitle("trigger efficiency");

        auto leg1=new TLegend(0.7,0.6,0.85,0.85);
        leg1->AddEntry(gH1_T1,"H1","P");
        leg1->AddEntry(gD2_T1,"D2","P");
        leg1->AddEntry(gHe_T1,"He3","P");
        leg1->AddEntry(gH3_T1,"H3","P");
	leg1->Draw();

        gH1_T2->SetMarkerStyle(8);
        gH1_T2->SetMarkerColor(1);
        gD2_T2->SetMarkerStyle(8);
        gD2_T2->SetMarkerColor(3);
        gHe_T2->SetMarkerStyle(8);
        gHe_T2->SetMarkerColor(4);
        gH3_T2->SetMarkerStyle(8);
        gH3_T2->SetMarkerColor(2);

        TCanvas *c2=new TCanvas("c2");
        TMultiGraph *mg1=new TMultiGraph();
        mg1->Add(gH1_T2);
        mg1->Add(gD2_T2);
        mg1->Add(gHe_T2);
        mg1->Add(gH3_T2);
        mg1->Draw("AP");
        mg1->GetXaxis()->SetTitle("kinematic");
        mg1->GetYaxis()->SetTitle("trigger efficiency");

        auto leg2=new TLegend(0.7,0.6,0.85,0.85);
        leg2->AddEntry(gH1_T2,"H1","P");
        leg2->AddEntry(gD2_T2,"D2","P");
        leg2->AddEntry(gHe_T2,"He3","P");
        leg2->AddEntry(gH3_T2,"H3","P");
	leg1->Draw();

        c1->Print("trigger.pdf[");
        c1->Print("trigger.pdf");
        c2->Print("trigger.pdf");
        c2->Print("trigger.pdf]");
     

}
