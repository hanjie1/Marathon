#include "SetCut.h"
#include "GetRunList.h"
#include "GetChain.h"
void piontoe(int nrun[], int size, int kin, Double_t &Ep_pet, Double_t &Ep_pec){
	TChain *T;
	T=GetChain(nrun,size,kin,"T");

        int KKin=0;  //kinematica series number
        if(kin<6)KKin=kin;
        if(kin==7)KKin=6;
        if(kin==9)KKin=7;
        if(kin==11)KKin=8;
        if(kin==13)KKin=9;
        if(kin==15)KKin=10;

	TH1F *hCK_e=new TH1F("hCK_e","electron CK distribution",500,0,15000);
	TH1F *hCK_pi=new TH1F("hCK_pi","pion CK distribution",500,0,15000);

	TH1F *hEp_e=new TH1F("hEp_e","electron E/p distribution",100,0,1.4);
	TH1F *hEp_pi=new TH1F("hEp_pi","pion E/p distribution",100,0,1.4);
	TH1F *hEp_pi_final=new TH1F("hEp_pi_final","pion E/p distribution",100,0,1.4);

	Ep_pet=0.0;
	Ep_pec=0.0;

	if(kin<16){
           TCut VZ=Form("L.tr.vz<%f && L.tr.vz>%f",vz_max[KKin],vz_min[KKin]);
	   TCut E_pi="(L.prl1.e+L.prl2.e)/(L.gold.p*1000)<0.2";
	   TCut CK_e="L.cer.asum_c>3500 && L.cer.asum_c<7000";
	   TCut CK_pi="L.cer.asum_c<200";
           TCut PS_e="L.prl1.e>(-16.0/15.0*L.prl2.e+1600.0)";


	   TCanvas *c1=new TCanvas("c1");
           T->Draw("L.cer.asum_c>>hCK_e",TRK+ACC+trigger1+VZ+beta+Ep+PS_e);
           T->Draw("L.cer.asum_c>>hCK_pi",TRK+ACC+trigger1+VZ+beta+E_pi,"same");
	   hCK_pi->SetLineColor(2);

	   int nbin1=hCK_e->FindBin(1500);
	   int nbin2=hCK_e->FindBin(15000);
	   Double_t CK_epass=hCK_e->Integral(nbin1,nbin2);
           Int_t max_bin = hCK_e->GetMaximumBin();
           Double_t nCK_e = hCK_e->GetBinContent(max_bin);

	   Double_t CK_pipass=hCK_pi->Integral(nbin1,nbin2);
           max_bin = hCK_pi->GetMaximumBin();
           Double_t nCK_pi = hCK_pi->GetBinContent(max_bin);

	   Double_t scale1=nCK_pi/nCK_e;
	   hCK_e->Scale(scale1);

	   Double_t CK_totalpi=hCK_pi->Integral();
	   Double_t CK_totale=hCK_e->Integral();
	   Double_t CK_pe=CK_totalpi/(scale1*CK_totale);
	   Double_t CK_pitoe=CK_pipass/(CK_epass*scale1);

	   c1->Update();
	   gPad->SetLogy();
	   c1->Print("PID_results/CK_kin15.png");	

//	   int nbin1=0,nbin2=0,max_bin=0;
	   TCanvas *c2=new TCanvas("c2");
           T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_pi",TRK+ACC+trigger1+VZ+beta+CK_pi);
           T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_e",TRK+ACC+trigger1+VZ+beta+CK_e,"same");
           //T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_pi",TRK+ACC+trigger1+VZ+beta+CK_pi,"same");
           T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_pi_final",TRK+ACC+trigger1+VZ+beta+CK_pi+Ep,"same");
	   hEp_pi->SetLineColor(2);

	   nbin1=hEp_e->FindBin(0.8);
	   nbin2=hEp_e->FindBin(1.2);
	   Double_t Ep_epass=hEp_e->Integral(nbin1,nbin2);
	   Double_t Ep_pipass=hEp_pi->Integral(nbin1,nbin2);

	   hEp_e->GetXaxis()->SetRangeUser(0.0,0.4);
           max_bin = hEp_e->GetMaximumBin();
	   hEp_e->GetXaxis()->SetRangeUser(0,1.5);
           Double_t nEp_e = hEp_e->GetBinContent(max_bin);

	   hEp_pi->GetXaxis()->SetRangeUser(0.0,0.4);
           max_bin = hEp_pi->GetMaximumBin();
	   hEp_pi->GetXaxis()->SetRangeUser(0.0,1.5);
           Double_t nEp_pi = hEp_pi->GetBinContent(max_bin);

	   Double_t scale2=nEp_pi/nEp_e;
	//   hEp_e->Scale(scale2);

	   Double_t Ep_totalpi=hEp_pi->Integral();
	   Double_t Ep_totale=hEp_e->Integral();
	   Ep_pet=Ep_totalpi/(scale2*Ep_totale);
	   Ep_pec=Ep_pipass/(Ep_epass*scale2);

	   c2->Update();
	   gPad->SetLogy();
	   c2->Print(Form("PID_results/Ep_%d.png",kin));	
	}

	return;
}


void Pion_contamination()
{
     int kin[12]={0,1,2,3,4,5,7,9,11,13,15,16};
     
     int H1_nrun[5]={2479,1213,1240,1286,2582};
     int D2_nrun[7]={2493,1214,1242,1288,2579,1350,1371};
     int He_nrun[7]={2531,1219,1236,1283,2576,1353,1373};
     int H3_nrun[7]={2503,1224,1245,1290,2577,1355,1377};

     int nkin9=5,nkin11=6,nkin13=20,nkin15=75,nkin16=150;

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

     Double_t H1_pet[12]={0.0},H1_pec[12]={0.0};
     Double_t D2_pet[12]={0.0},D2_pec[12]={0.0};
     Double_t He_pet[12]={0.0},He_pec[12]={0.0};
     Double_t H3_pet[12]={0.0},H3_pec[12]={0.0};

     ofstream outfile;
     outfile.open("PID_results/pitoe.txt");
     outfile<<"kin    H1_pet    H1_pec    D2_pet    D2_pec    He_pet    He_pec    H3_pet    H3_pec"<<endl;
     for(int ii=10;ii<11;ii++){
	 int nrun[1]={0};
	 Double_t pet=0.0,pec=0.0;
/*	 if(ii<5){
	    nrun[0]=H1_nrun[ii];
	    piontoe(nrun,1,kin[ii],pet,pec);
	    H1_pet[ii]=pet;
	    H1_pec[ii]=pec;
	 }
	 if(ii<7){
	    nrun[0]=D2_nrun[ii];
	    piontoe(nrun,1,kin[ii],pet,pec);
	 }
	 if(kin[ii]==9)piontoe(D2_kin9,nkin9,kin[ii],pet,pec);
	 if(kin[ii]==11)piontoe(D2_kin11,nkin11,kin[ii],pet,pec);
	 if(kin[ii]==13)piontoe(D2_kin13,nkin13,kin[ii],pet,pec);
*/	 if(kin[ii]==15)piontoe(D2_kin15,nkin15,kin[ii],pet,pec);
	 if(kin[ii]==16)piontoe(D2_kin16,nkin16,kin[ii],pet,pec);
	 D2_pet[ii]=pet;
	 D2_pec[ii]=pec;
/*
	 if(ii<7){
	    nrun[0]=He_nrun[ii];
	    piontoe(nrun,1,kin[ii],pet,pec);
	 }
	 if(kin[ii]==9)piontoe(He_kin9,nkin9,kin[ii],pet,pec);
	 if(kin[ii]==11)piontoe(He_kin11,nkin11,kin[ii],pet,pec);
	 if(kin[ii]==13)piontoe(He_kin13,nkin13,kin[ii],pet,pec);
	 if(kin[ii]==15)piontoe(He_kin15,nkin15,kin[ii],pet,pec);
	 if(kin[ii]==16)piontoe(He_kin16,nkin16,kin[ii],pet,pec);
	 He_pet[ii]=pet;
	 He_pec[ii]=pec;

	 if(ii<7){
	    nrun[0]=H3_nrun[ii];
	    piontoe(nrun,1,kin[ii],pet,pec);
	 }
	 if(kin[ii]==9)piontoe(H3_kin9,nkin9,kin[ii],pet,pec);
	 if(kin[ii]==11)piontoe(H3_kin11,nkin11,kin[ii],pet,pec);
	 if(kin[ii]==13)piontoe(H3_kin13,nkin13,kin[ii],pet,pec);
	 if(kin[ii]==15)piontoe(H3_kin15,nkin15,kin[ii],pet,pec);
	 if(kin[ii]==16)piontoe(H3_kin16,nkin16,kin[ii],pet,pec);
	 H3_pet[ii]=pet;
	 H3_pec[ii]=pec;
*/
	outfile<<kin[ii]<<"  "<<H1_pet[ii]<<"  "<<H1_pec[ii]<<"  "<<D2_pet[ii]<<"  "<<D2_pec[ii]<<"  "
		<<He_pet[ii]<<"  "<<He_pec[ii]<<"  "<<H3_pet[ii]<<"  "<<H3_pec[ii]<<endl;
     }
 
     outfile.close();

 
}
