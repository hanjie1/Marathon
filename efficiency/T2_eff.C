#include "GetChain.h"
#include "SetCut.h"
#include "GetRunList.h"
struct EFF_VAR{
      Double_t T1_eff;
      Double_t T1_err;
      Double_t T2_eff;
      Double_t T2_err;
};

void trigger(int run_number[],int size, int kin, EFF_VAR &TT){
     TString TreeName="T";
     TChain *T=GetChain(run_number,size,kin,TreeName);

     Double_t NT1=0.0,NT2=0.0,NT3=0.0,NT2_2=0.0;

     if(kin<16){
        TCut EE="(L.prl1.e+L.prl2.e)/3100.0>0.8 && (L.prl1.e+L.prl2.e)/3100.0<1.2";
        TCut PS="L.prl1.e>(-16.0/15.0*L.prl2.e+1600.0)";
        TCut CK = "L.cer.asum_c>2200 && L.cer.asum_c<8000";

        NT3=T->GetEntries(EE+CK+trigger3);
        NT2=T->GetEntries(EE+CK+trigger2);

        NT1=T->GetEntries(EE+trigger1+beta+ACC+PS+TRK);
        NT2_2=T->GetEntries(EE+trigger2+beta+ACC+PS+TRK);
     }

     if(kin==16){
        TCut EE="(R.ps.e+R.sh.e)/2900.0>0.8 && (R.ps.e+R.sh.e)/2900.0<1.2";
        TCut CK_R= "R.cer.asum_c>3000 && R.cer.asum_c<9000";
        TCut PS="R.ps.e>400 && R.sh.e>1000";

        NT3=T->GetEntries(EE+CK_R+trigger6+PS);
        NT2=T->GetEntries(EE+CK_R+trigger5+PS);

        NT1=T->GetEntries(EE+trigger4+beta_R+ACC_R+PS+TRK_R);
        NT2_2=T->GetEntries(EE+trigger5+beta_R+ACC_R+PS+TRK_R);
     }

     Double_t eff1=NT2/NT3;
     Double_t eff1_err=sqrt(1.0/NT2-1.0/NT3)*eff1;

     Double_t eff2=NT2_2/NT1;
     Double_t eff2_err=sqrt(1.0/NT2_2-1.0/NT1)*eff2;

     TT.T1_eff=eff1;
     TT.T1_err=eff1_err;

     TT.T2_eff=eff2*eff1;
     TT.T2_err=sqrt(eff1*eff1*eff2_err*eff2_err+eff2*eff2*eff1_err*eff1_err);

     cout<<kin<<": "<<eff1<<"  "<<eff1_err<<"  "<<eff2<<"  "<<eff2_err<<endl;

     return;
}

void Trigger_eff(){
     int kin[4]={1,2,3,4};
 
     int C_nrun[4]={1207,1233,1279,2568};
   
     TGraphErrors *gC_T1=new TGraphErrors();
     TGraphErrors *gC_T2=new TGraphErrors();

     EFF_VAR H1_EFF[4];

     int nnp=0;
     int nnrun[1]={0};
     for(int ii=0;ii<4;ii++){
         Double_t KKin=1.0*kin[ii];
	 EFF_VAR C_single;
	 nnrun[0]=C_nrun[ii];
         trigger(nnrun,1,kin[ii],C_single);
	 C_EFF[ii]=C_single;
            gC_T1->SetPoint(nnp,KKin,C_single.T1_eff);
            gC_T1->SetPointError(nnp,0.0,C_single.T1_err);
            gC_T2->SetPoint(nnp,KKin,C_single.T2_eff);
            gC_T2->SetPointError(nnp,0.0,C_single.T2_err);
	 nnp++;
     }

     Double_t R_Dp[5]={0.0},R_HeD[12]={0.0},R_H3D[12]={0.0},R_H3He[12]={0.0};
     Double_t RDp_err[5]={0.0},RHeD_err[12]={0.0},RH3D_err[12]={0.0},RH3He_err[12]={0.0};
     TGraphErrors *gDp_T2=new TGraphErrors();
     TGraphErrors *gHeD_T2=new TGraphErrors();
     TGraphErrors *gH3D_T2=new TGraphErrors();
     TGraphErrors *gH3He_T2=new TGraphErrors();

     nnp=0;
     for(int ii=8;ii<11;ii++){
	if(ii<5){
	   R_Dp[ii]=D2_EFF[ii].T2_eff/H1_EFF[ii].T2_eff;
	   RDp_err[ii]=sqrt(D2_EFF[ii].T2_err*D2_EFF[ii].T2_err/(H1_EFF[ii].T2_eff*H1_EFF[ii].T2_eff)
			+pow(D2_EFF[ii].T2_eff,2)/pow(H1_EFF[ii].T2_eff,4)*pow(H1_EFF[ii].T2_err,2));
	   gDp_T2->SetPoint(nnp,kin[ii],R_Dp[ii]);
	   gDp_T2->SetPointError(nnp,0.0,RDp_err[ii]);
	}

	R_HeD[ii]=He_EFF[ii].T2_eff/D2_EFF[ii].T2_eff;
	RHeD_err[ii]=sqrt(He_EFF[ii].T2_err*He_EFF[ii].T2_err/(D2_EFF[ii].T2_eff*D2_EFF[ii].T2_eff)
	           +pow(He_EFF[ii].T2_eff,2)/pow(D2_EFF[ii].T2_eff,4)*pow(D2_EFF[ii].T2_err,2));
	gHeD_T2->SetPoint(nnp,kin[ii],R_HeD[ii]);
	gHeD_T2->SetPointError(nnp,0.0,RHeD_err[ii]);

	R_H3D[ii]=H3_EFF[ii].T2_eff/D2_EFF[ii].T2_eff;
	RH3D_err[ii]=sqrt(H3_EFF[ii].T2_err*H3_EFF[ii].T2_err/(D2_EFF[ii].T2_eff*D2_EFF[ii].T2_eff)
	           +pow(H3_EFF[ii].T2_eff,2)/pow(D2_EFF[ii].T2_eff,4)*pow(D2_EFF[ii].T2_err,2));
	gH3D_T2->SetPoint(nnp,kin[ii],R_H3D[ii]);
	gH3D_T2->SetPointError(nnp,0.0,RH3D_err[ii]);

	R_H3He[ii]=H3_EFF[ii].T2_eff/He_EFF[ii].T2_eff;
	RH3He_err[ii]=sqrt(H3_EFF[ii].T2_err*H3_EFF[ii].T2_err/(He_EFF[ii].T2_eff*He_EFF[ii].T2_eff)
	           +pow(H3_EFF[ii].T2_eff,2)/pow(He_EFF[ii].T2_eff,4)*pow(He_EFF[ii].T2_err,2));
	gH3He_T2->SetPoint(nnp,kin[ii],R_H3He[ii]);
	gH3He_T2->SetPointError(nnp,0.0,RH3He_err[ii]);
	nnp++;
     }

        gH1_T1->SetMarkerStyle(8);
        gH1_T1->SetMarkerColor(1);
        gD2_T1->SetMarkerStyle(8);
        gD2_T1->SetMarkerColor(3);
        gHe_T1->SetMarkerStyle(8);
        gHe_T1->SetMarkerColor(4);
        gH3_T1->SetMarkerStyle(8);
        gH3_T1->SetMarkerColor(2);

        gDp_T2->SetMarkerStyle(8);
        gDp_T2->SetMarkerColor(1);
        gHeD_T2->SetMarkerStyle(8);
        gHeD_T2->SetMarkerColor(3);
        gH3D_T2->SetMarkerStyle(8);
        gH3D_T2->SetMarkerColor(4);
        gH3He_T2->SetMarkerStyle(8);
        gH3He_T2->SetMarkerColor(2);


        TCanvas *c1=new TCanvas("c1","c1",600,800);
        TMultiGraph *mg=new TMultiGraph();
        mg->Add(gH1_T1);
        mg->Add(gD2_T1);
        mg->Add(gHe_T1);
        mg->Add(gH3_T1);
        mg->Draw("AP");
        mg->GetXaxis()->SetTitle("kinematic");
        mg->GetYaxis()->SetTitle("T1 efficiency");

        auto leg1=new TLegend(0.2,0.2,0.35,0.45);
        leg1->AddEntry(gH1_T1,"{}^{1}H","P");
        leg1->AddEntry(gD2_T1,"{}^{2}H","P");
        leg1->AddEntry(gHe_T1,"{}^{3}He","P");
        leg1->AddEntry(gH3_T1,"{}^{3}H","P");
	leg1->Draw();

        gH1_T2->SetMarkerStyle(8);
        gH1_T2->SetMarkerColor(1);
        gD2_T2->SetMarkerStyle(8);
        gD2_T2->SetMarkerColor(3);
        gHe_T2->SetMarkerStyle(8);
        gHe_T2->SetMarkerColor(4);
        gH3_T2->SetMarkerStyle(8);
        gH3_T2->SetMarkerColor(2);

        TCanvas *c2=new TCanvas("c2","c2",600,800);
        TMultiGraph *mg1=new TMultiGraph();
        mg1->Add(gH1_T2);
        mg1->Add(gD2_T2);
        mg1->Add(gHe_T2);
        mg1->Add(gH3_T2);
        mg1->Draw("AP");
        mg1->GetXaxis()->SetTitle("kinematic");
        mg1->GetYaxis()->SetTitle("T2 efficiency");

        auto leg2=new TLegend(0.2,0.2,0.35,0.45);
        leg2->AddEntry(gH1_T2,"{}^{1}H","P");
        leg2->AddEntry(gD2_T2,"{}^{2}H","P");
        leg2->AddEntry(gHe_T2,"{}^{3}He","P");
        leg2->AddEntry(gH3_T2,"{}^{3}H","P");
	leg2->Draw();

        TCanvas *c3=new TCanvas("c3","c3",600,800);
        TMultiGraph *mg3=new TMultiGraph();
        mg3->Add(gDp_T2);
        mg3->Add(gHeD_T2);
        mg3->Add(gH3D_T2);
        mg3->Add(gH3He_T2);
        mg3->Draw("AP");
        mg3->GetXaxis()->SetTitle("kinematic");
        mg3->GetYaxis()->SetTitle("T2 efficiency ratio");

        auto leg3=new TLegend(0.2,0.2,0.35,0.45);
        leg3->AddEntry(gDp_T2,"{}^{2}H/{}^{1}H","P");
        leg3->AddEntry(gHeD_T2,"{}^{3}He/{}^{2}H","P");
        leg3->AddEntry(gH3D_T2,"{}^{3}H/{}^{2}H","P");
        leg3->AddEntry(gH3He_T2,"{}^{3}H/{}^{3}He","P");
	leg3->Draw();
/*

        c1->Print("trigger.pdf[");
        c1->Print("trigger.pdf");
        c2->Print("trigger.pdf");
        c3->Print("trigger.pdf");
        c3->Print("trigger.pdf]");
*/   

}
