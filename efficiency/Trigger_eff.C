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
     int kin[12]={0,1,2,3,4,5,7,9,11,13,15,16};
 
     int H1_nrun[5]={2479,1213,1240,1286,2582};
     int D2_nrun[7]={2493,1214,1242,1288,2579,1350,1371};
     int He_nrun[7]={2531,1219,1236,1283,2576,1353,1373};
     int H3_nrun[7]={2503,1224,1245,1290,2577,1355,1377};
   
     int nkin9=5,nkin11=6,nkin13=20,nkin15=60,nkin16=150;
 
     int D2_kin9[nkin9];
     int D2_kin11[nkin11];
     int D2_kin13[nkin13];
     int D2_kin15[nkin15];
     int D2_kin16[nkin16];

     GetRunList(D2_kin9,nkin9,9,"D2");
     GetRunList(D2_kin11,nkin11,11,"D2");
     GetRunList(D2_kin13,nkin13,13,"D2");
     GetRunList(D2_kin15,nkin15,15,"D2");
     GetRunList(D2_kin16,nkin16,16,"D2");

     int He_kin9[nkin9];
     int He_kin11[nkin11];
     int He_kin13[nkin13];
     int He_kin15[nkin15];
     int He_kin16[nkin16];

     GetRunList(He_kin9,nkin9,9,"He3");
     GetRunList(He_kin11,nkin11,11,"He3");
     GetRunList(He_kin13,nkin13,13,"He3");
     GetRunList(He_kin15,nkin15,15,"He3");
     GetRunList(He_kin16,nkin16,16,"He3");

     int H3_kin9[nkin9];
     int H3_kin11[nkin11];
     int H3_kin13[nkin13];
     int H3_kin15[nkin15];
     int H3_kin16[nkin16];

     GetRunList(H3_kin9,nkin9,9,"H3");
     GetRunList(H3_kin11,nkin11,11,"H3");
     GetRunList(H3_kin13,nkin13,13,"H3");
     GetRunList(H3_kin15,nkin15,15,"H3");
     GetRunList(H3_kin16,nkin16,16,"H3");
 
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

     int nnrun[1]={0};
     for(int ii=0;ii<12;ii++){
         Double_t KKin=1.0*kin[ii];
         if(ii<5){
	    EFF_VAR H1_single;
	    nnrun[0]=H1_nrun[ii];
            trigger(nnrun,1,kin[ii],H1_single);
	    H1_EFF[ii]=H1_single;
            gH1_T1->SetPoint(ii,KKin,H1_single.T1_eff);
            gH1_T1->SetPointError(ii,0.0,H1_single.T1_err);
            gH1_T2->SetPoint(ii,KKin,H1_single.T2_eff);
            gH1_T2->SetPointError(ii,0.0,H1_single.T2_err);
         }

	 EFF_VAR D2_single;
         if(ii<7){
	    nnrun[0]=D2_nrun[ii];
	    trigger(nnrun,1,kin[ii],D2_single);
	 }
         if(kin[ii]==9)trigger(D2_kin9,nkin9,kin[ii],D2_single);
         if(kin[ii]==11)trigger(D2_kin11,nkin11,kin[ii],D2_single);
         if(kin[ii]==13)trigger(D2_kin13,nkin13,kin[ii],D2_single);
         if(kin[ii]==15)trigger(D2_kin15,nkin15,kin[ii],D2_single);
         if(kin[ii]==16)trigger(D2_kin16,nkin16,kin[ii],D2_single);
	 D2_EFF[ii]=D2_single;
         gD2_T1->SetPoint(ii,KKin,D2_single.T1_eff);
         gD2_T1->SetPointError(ii,0.0,D2_single.T1_err);
         gD2_T2->SetPoint(ii,KKin,D2_single.T2_eff);
         gD2_T2->SetPointError(ii,0.0,D2_single.T2_err);

	 EFF_VAR He_single;
         if(ii<7){
	    nnrun[0]=He_nrun[ii];
	    trigger(nnrun,1,kin[ii],He_single);
	 }
         if(kin[ii]==9)trigger(He_kin9,nkin9,kin[ii],He_single);
         if(kin[ii]==11)trigger(He_kin11,nkin11,kin[ii],He_single);
         if(kin[ii]==13)trigger(He_kin13,nkin13,kin[ii],He_single);
         if(kin[ii]==15)trigger(He_kin15,nkin15,kin[ii],He_single);
         if(kin[ii]==16)trigger(He_kin16,nkin16,kin[ii],He_single);

	 He_EFF[ii]=He_single;
         gHe_T1->SetPoint(ii,KKin,He_single.T1_eff);
         gHe_T1->SetPointError(ii,0.0,He_single.T1_err);
         gHe_T2->SetPoint(ii,KKin,He_single.T2_eff);
         gHe_T2->SetPointError(ii,0.0,He_single.T2_err);

	 EFF_VAR H3_single;
         if(ii<7){
	    nnrun[0]=H3_nrun[ii];
	    trigger(nnrun,1,kin[ii],H3_single);
	 }
         if(kin[ii]==9)trigger(H3_kin9,nkin9,kin[ii],H3_single);
         if(kin[ii]==11)trigger(H3_kin11,nkin11,kin[ii],H3_single);
         if(kin[ii]==13)trigger(H3_kin13,nkin13,kin[ii],H3_single);
         if(kin[ii]==15)trigger(H3_kin15,nkin15,kin[ii],H3_single);
         if(kin[ii]==16)trigger(H3_kin16,nkin16,kin[ii],H3_single);

	 H3_EFF[ii]=H3_single;
         gH3_T1->SetPoint(ii,KKin,H3_single.T1_eff);
         gH3_T1->SetPointError(ii,0.0,H3_single.T1_err);
         gH3_T2->SetPoint(ii,KKin,H3_single.T2_eff);
         gH3_T2->SetPointError(ii,0.0,H3_single.T2_err);
     }

     Double_t R_Dp[5]={0.0},R_HeD[12]={0.0},R_H3D[12]={0.0},R_H3He[12]={0.0};
     Double_t RDp_err[5]={0.0},RHeD_err[12]={0.0},RH3D_err[12]={0.0},RH3He_err[12]={0.0};
     TGraphErrors *gDp_T2=new TGraphErrors(5);
     TGraphErrors *gHeD_T2=new TGraphErrors(12);
     TGraphErrors *gH3D_T2=new TGraphErrors(12);
     TGraphErrors *gH3He_T2=new TGraphErrors(12);

     for(int ii=0;ii<12;ii++){
	if(ii<5){
	   R_Dp[ii]=D2_EFF[ii].T2_eff/H1_EFF[ii].T2_eff;
	   RDp_err[ii]=sqrt(D2_EFF[ii].T2_err*D2_EFF[ii].T2_err/(H1_EFF[ii].T2_eff*H1_EFF[ii].T2_eff)
			+pow(D2_EFF[ii].T2_eff,2)/pow(H1_EFF[ii].T2_eff,4)*pow(H1_EFF[ii].T2_err,2));
	   gDp_T2->SetPoint(ii,kin[ii],R_Dp[ii]);
	   gDp_T2->SetPointError(ii,0.0,RDp_err[ii]);
	}

	R_HeD[ii]=He_EFF[ii].T2_eff/D2_EFF[ii].T2_eff;
	RHeD_err[ii]=sqrt(He_EFF[ii].T2_err*He_EFF[ii].T2_err/(D2_EFF[ii].T2_eff*D2_EFF[ii].T2_eff)
	           +pow(He_EFF[ii].T2_eff,2)/pow(D2_EFF[ii].T2_eff,4)*pow(D2_EFF[ii].T2_err,2));
	gHeD_T2->SetPoint(ii,kin[ii],R_HeD[ii]);
	gHeD_T2->SetPointError(ii,0.0,RHeD_err[ii]);

	R_H3D[ii]=H3_EFF[ii].T2_eff/D2_EFF[ii].T2_eff;
	RH3D_err[ii]=sqrt(H3_EFF[ii].T2_err*H3_EFF[ii].T2_err/(D2_EFF[ii].T2_eff*D2_EFF[ii].T2_eff)
	           +pow(H3_EFF[ii].T2_eff,2)/pow(D2_EFF[ii].T2_eff,4)*pow(D2_EFF[ii].T2_err,2));
	gH3D_T2->SetPoint(ii,kin[ii],R_H3D[ii]);
	gH3D_T2->SetPointError(ii,0.0,RH3D_err[ii]);

	R_H3He[ii]=H3_EFF[ii].T2_eff/He_EFF[ii].T2_eff;
	RH3He_err[ii]=sqrt(H3_EFF[ii].T2_err*H3_EFF[ii].T2_err/(He_EFF[ii].T2_eff*He_EFF[ii].T2_eff)
	           +pow(H3_EFF[ii].T2_eff,2)/pow(He_EFF[ii].T2_eff,4)*pow(He_EFF[ii].T2_err,2));
	gH3He_T2->SetPoint(ii,kin[ii],R_H3He[ii]);
	gH3He_T2->SetPointError(ii,0.0,RH3He_err[ii]);
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
        leg2->AddEntry(gH1_T2,"H1","P");
        leg2->AddEntry(gD2_T2,"D2","P");
        leg2->AddEntry(gHe_T2,"He3","P");
        leg2->AddEntry(gH3_T2,"H3","P");
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
        leg3->AddEntry(gDp_T2,"D/p","P");
        leg3->AddEntry(gHeD_T2,"He3/D2","P");
        leg3->AddEntry(gH3D_T2,"H3/D2","P");
        leg3->AddEntry(gH3He_T2,"H3/He","P");
	leg3->Draw();


        c1->Print("trigger.pdf[");
        c1->Print("trigger.pdf");
        c2->Print("trigger.pdf");
        c3->Print("trigger.pdf");
        c3->Print("trigger.pdf]");
     

}
