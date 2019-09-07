#include "SetCut.h"
#include "GetRunList.h"
#include "GetChain.h"
void piontoe(int nrun[], int size, int kin, TString target, Double_t &Ep_pet, Double_t &Ep_pec){
	TChain *T;
	T=GetChain(nrun,size,kin,"T");

        int KKin=0;  //kinematica series number
        if(kin<6)KKin=kin;
        if(kin==7)KKin=6;
        if(kin==9)KKin=7;
        if(kin==11)KKin=8;
        if(kin==13)KKin=9;
        if(kin==15)KKin=10;

        delete gROOT->FindObject("hEp_e");
        delete gROOT->FindObject("hEp_pi");
        delete gROOT->FindObject("hEp_pi_final");

	TH1F *hEp_e=new TH1F("hEp_e","electron E/p distribution",100,0,1.4);
	TH1F *hEp_pi=new TH1F("hEp_pi","pion E/p distribution",100,0,1.4);
	TH1F *hEp_pi_final=new TH1F("hEp_pi_final","pion E/p distribution",100,0,1.4);

	Ep_pet=0.0;
	Ep_pec=0.0;

	TCut VZ,CK_e,CK_pi;
        TCanvas *c1=new TCanvas("c1");
	if(kin<16){
           VZ=Form("L.tr.vz<%f && L.tr.vz>%f",vz_max[KKin],vz_min[KKin]);
	   CK_e="L.cer.asum_c>1500";
	   CK_pi="L.cer.asum_c<200";

           T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_e",TRK+ACC+trigger2+VZ+beta+CK_e+W2,"HIST");
           T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_pi",TRK+ACC+trigger1+VZ+beta+CK_pi+W2,"same HIST");
	}
	else{
	   CK_e="R.cer.asum_c>2000";
	   CK_pi="R.cer.asum_c<200";

           T->Draw("(R.ps.e+R.sh.e)/(1000*R.gold.p)>>hEp_e",TRK_R+ACC_R+trigger5+VZ_R+beta_R+CK_e+W2_R,"HIST");
           T->Draw("(R.ps.e+R.sh.e)/(1000*R.gold.p)>>hEp_pi",TRK_R+ACC_R+trigger4+VZ_R+beta_R+CK_pi+W2_R,"same HIST");

	}


	   hEp_e->SetLineColor(2);
	   hEp_e->SetLineWidth(3);
	   hEp_pi->SetLineColor(4);
	   hEp_pi->SetLineWidth(3);

	   Double_t tEp_epass=hEp_e->Integral();
	   Double_t tEp_pi=hEp_e->Integral();

	   hEp_e->GetXaxis()->SetRangeUser(0.02,0.05);
           int min1 = hEp_e->GetMinimumBin();
           Double_t tmp_max1 = hEp_e->GetBinContent(min1);

	   hEp_e->GetXaxis()->SetRangeUser(0.05,0.15);
           int max2 = hEp_e->GetMaximumBin();
           Double_t tmp_max2 = hEp_e->GetBinContent(max2);
	   hEp_e->GetXaxis()->SetRangeUser(0,1.4);

	   if(tmp_max2>tmp_max1){ //has pion contamination
              if(kin<16)
		T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_pi_final",TRK+ACC+trigger1+VZ+beta+CK_pi+Ep,"same HIST");
	      else
                T->Draw("(R.ps.e+R.sh.e)/(1000*R.gold.p)>>hEp_pi_final",TRK_R+ACC_R+trigger4+VZ_R+beta_R+CK_pi+Ep_R,"same HIST");
	      hEp_pi_final->SetLineColor(4);
	      hEp_pi_final->SetLineWidth(3);
              Double_t tEp_pipass=hEp_pi_final->Integral();

	      hEp_pi->GetXaxis()->SetRangeUser(0.0,0.4);
              int max_bin = hEp_pi->GetMaximumBin();
	      hEp_pi->GetXaxis()->SetRangeUser(0.0,1.4);
              Double_t nEp_pi = hEp_pi->GetBinContent(max_bin);
	      Double_t scale=tmp_max2/nEp_pi;
	      hEp_pi->Scale(scale);
	      hEp_pi_final->Scale(scale);

	      Ep_pec=tEp_pipass*scale/tEp_epass;
	      Ep_pet=tEp_pi*scale/tEp_epass;
	      hEp_pi->SetFillColor(38);
	      hEp_pi_final->SetFillColor(9);
	   }
	   else{
              hEp_pi->GetXaxis()->SetRangeUser(0.0,0.4);
              int max_bin = hEp_pi->GetMaximumBin();
              hEp_pi->GetXaxis()->SetRangeUser(0.0,1.4);
              Double_t nEp_pi = hEp_pi->GetBinContent(max_bin);

	      if(nEp_pi>tmp_max1){
	         hEp_e->GetYaxis()->SetRangeUser(1,nEp_pi+1000);
	      }
	   }

	  
	   gPad->SetLogy();
	   c1->Print(Form("PID_results/Ep_%s_%d.png",target.Data(),kin));	
	   delete c1;	

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
     TString target;
     for(int ii=0;ii<12;ii++){
	 int nrun[1]={0};
	 Double_t pet=0.0,pec=0.0;
	 target="H1";
	 if(ii<5){
	    nrun[0]=H1_nrun[ii];
	    piontoe(nrun,1,kin[ii],target,pet,pec);
	    H1_pet[ii]=pet;
	    H1_pec[ii]=pec;
	 }
	 pet=0.0;pec=0.0;
	 target="D2";
	 if(ii<7){
	    nrun[0]=D2_nrun[ii];
	    piontoe(nrun,1,kin[ii],target,pet,pec);
	 }
	 if(kin[ii]==9)piontoe(D2_kin9,nkin9,kin[ii],target,pet,pec);
	 if(kin[ii]==11)piontoe(D2_kin11,nkin11,kin[ii],target,pet,pec);
	 if(kin[ii]==13)piontoe(D2_kin13,nkin13,kin[ii],target,pet,pec);
	 if(kin[ii]==15)piontoe(D2_kin15,nkin15,kin[ii],target,pet,pec);
	 if(kin[ii]==16)piontoe(D2_kin16,nkin16,kin[ii],target,pet,pec);
	 D2_pet[ii]=pet;
	 D2_pec[ii]=pec;
	 
	 pet=0.0;pec=0.0;
	 target="He3";
	 if(ii<7){
	    nrun[0]=He_nrun[ii];
	    piontoe(nrun,1,kin[ii],target,pet,pec);
	 }
	 if(kin[ii]==9)piontoe(He_kin9,nkin9,kin[ii],target,pet,pec);
	 if(kin[ii]==11)piontoe(He_kin11,nkin11,kin[ii],target,pet,pec);
	 if(kin[ii]==13)piontoe(He_kin13,nkin13,kin[ii],target,pet,pec);
	 if(kin[ii]==15)piontoe(He_kin15,nkin15,kin[ii],target,pet,pec);
	 if(kin[ii]==16)piontoe(He_kin16,nkin16,kin[ii],target,pet,pec);
	 He_pet[ii]=pet;
	 He_pec[ii]=pec;

	 pet=0.0;pec=0.0;
	 target="H3";
	 if(ii<7){
	    nrun[0]=H3_nrun[ii];
	    piontoe(nrun,1,kin[ii],target,pet,pec);
	 }
	 if(kin[ii]==9)piontoe(H3_kin9,nkin9,kin[ii],target,pet,pec);
	 if(kin[ii]==11)piontoe(H3_kin11,nkin11,kin[ii],target,pet,pec);
	 if(kin[ii]==13)piontoe(H3_kin13,nkin13,kin[ii],target,pet,pec);
	 if(kin[ii]==15)piontoe(H3_kin15,nkin15,kin[ii],target,pet,pec);
	 if(kin[ii]==16)piontoe(H3_kin16,nkin16,kin[ii],target,pet,pec);
	 H3_pet[ii]=pet;
	 H3_pec[ii]=pec;

	outfile<<kin[ii]<<"  "<<H1_pet[ii]<<"  "<<H1_pec[ii]<<"  "<<D2_pet[ii]<<"  "<<D2_pec[ii]<<"  "
		<<He_pet[ii]<<"  "<<He_pec[ii]<<"  "<<H3_pet[ii]<<"  "<<H3_pec[ii]<<endl;
     }
 
     outfile.close();

 
}
