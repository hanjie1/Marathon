#include "SetCut.h"
#include "GetRunList.h"
#include "GetChain.h"
void piontoe(int nrun[], int size, int kin, TString target, Double_t &Ep_pec, Double_t &pecErr){
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

	Double_t ep_cut1=0.0;
	Double_t ep_cut2=0.7;

	TCut VZ,CK_e,CK_pi;
        TCanvas *c1=new TCanvas("c1");
	if(kin<16){
           VZ=Form("L.tr.vz<%f && L.tr.vz>%f",vz_max[KKin],vz_min[KKin]);
	   CK_e="L.cer.asum_c>1500";
	   CK_pi="L.cer.asum_c<200";

           T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_e",TRK+ACC+trigger2+VZ+beta+CK_e+W2,"HIST");
           T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_pi",TRK+ACC+trigger1+VZ+beta+CK_pi+W2,"same HIST");
	   ep_cut1=0.1; 
	}
	else{
	   CK_e="R.cer.asum_c>2000";
	   CK_pi="R.cer.asum_c<200";

           T->Draw("(R.ps.e+R.sh.e)/(1000*R.gold.p)>>hEp_e",TRK_R+ACC_R+trigger5+VZ_R+beta_R+CK_e+W2_R,"HIST");
           T->Draw("(R.ps.e+R.sh.e)/(1000*R.gold.p)>>hEp_pi",TRK_R+ACC_R+trigger4+VZ_R+beta_R+CK_pi+W2_R,"same HIST");
	   ep_cut1=0.15; 
	}


	   hEp_e->SetLineColor(2);
	   hEp_e->SetLineWidth(3);
	   hEp_pi->SetLineColor(4);
	   hEp_pi->SetLineWidth(3);

           int nbin1=hEp_e->FindBin(ep_cut1);
           int nbin2=hEp_e->FindBin(ep_cut2);
           int nbin3=hEp_e->FindBin(1.4);

           Double_t N3 = hEp_e->Integral(nbin1,nbin2);
           Double_t N4 = hEp_e->Integral(nbin2,nbin3);

           nbin1=hEp_pi->FindBin(ep_cut1);
           nbin2=hEp_pi->FindBin(ep_cut2);
           nbin3=hEp_pi->FindBin(1.4);

           Double_t N1 = hEp_pi->Integral(nbin1,nbin2);
           Double_t N2 = hEp_pi->Integral(nbin2,nbin3);
         
  	   Ep_pec=(N3/N1*N2)/N4;    
	   pecErr=Ep_pec*sqrt(1.0/N1+1.0/N2+1.0/N3+1.0/N4);

	   delete c1;	

	return;
}


void Pion_contamination_Javier()
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

     Double_t H1_pecErr[12]={0.0},H1_pec[12]={0.0};
     Double_t D2_pecErr[12]={0.0},D2_pec[12]={0.0};
     Double_t He_pecErr[12]={0.0},He_pec[12]={0.0};
     Double_t H3_pecErr[12]={0.0},H3_pec[12]={0.0};

     ofstream outfile;
     outfile.open("PID_results/pitoe_Javier.txt");
     outfile<<"kin    H1_pec   D2_pec   He_pec    H3_pec"<<endl;
     TString target;
     for(int ii=0;ii<12;ii++){
	 int nrun[1]={0};
	 Double_t pecErr=0.0,pec=0.0;
	 target="H1";
	 if(ii<5){
	    nrun[0]=H1_nrun[ii];
	    piontoe(nrun,1,kin[ii],target,pec,pecErr);
	    H1_pecErr[ii]=pecErr;
	    H1_pec[ii]=pec;
	 }
	 pecErr=0.0;pec=0.0;
	 target="D2";
	 if(ii<7){
	    nrun[0]=D2_nrun[ii];
	    piontoe(nrun,1,kin[ii],target,pec,pecErr);
	 }
	 if(kin[ii]==9)piontoe(D2_kin9,nkin9,kin[ii],target,pec,pecErr);
	 if(kin[ii]==11)piontoe(D2_kin11,nkin11,kin[ii],target,pec,pecErr);
	 if(kin[ii]==13)piontoe(D2_kin13,nkin13,kin[ii],target,pec,pecErr);
	 if(kin[ii]==15)piontoe(D2_kin15,nkin15,kin[ii],target,pec,pecErr);
	 if(kin[ii]==16)piontoe(D2_kin16,nkin16,kin[ii],target,pec,pecErr);
	 D2_pecErr[ii]=pecErr;
	 D2_pec[ii]=pec;
	 
	 pecErr=0.0;pec=0.0;
	 target="He3";
	 if(ii<7){
	    nrun[0]=He_nrun[ii];
	    piontoe(nrun,1,kin[ii],target,pec,pecErr);
	 }
	 if(kin[ii]==9)piontoe(He_kin9,nkin9,kin[ii],target,pec,pecErr);
	 if(kin[ii]==11)piontoe(He_kin11,nkin11,kin[ii],target,pec,pecErr);
	 if(kin[ii]==13)piontoe(He_kin13,nkin13,kin[ii],target,pec,pecErr);
	 if(kin[ii]==15)piontoe(He_kin15,nkin15,kin[ii],target,pec,pecErr);
	 if(kin[ii]==16)piontoe(He_kin16,nkin16,kin[ii],target,pec,pecErr);
	 He_pecErr[ii]=pecErr;
	 He_pec[ii]=pec;

	 pecErr=0.0;pec=0.0;
	 target="H3";
	 if(ii<7){
	    nrun[0]=H3_nrun[ii];
	    piontoe(nrun,1,kin[ii],target,pec,pecErr);
	 }
	 if(kin[ii]==9)piontoe(H3_kin9,nkin9,kin[ii],target,pec,pecErr);
	 if(kin[ii]==11)piontoe(H3_kin11,nkin11,kin[ii],target,pec,pecErr);
	 if(kin[ii]==13)piontoe(H3_kin13,nkin13,kin[ii],target,pec,pecErr);
	 if(kin[ii]==15)piontoe(H3_kin15,nkin15,kin[ii],target,pec,pecErr);
	 if(kin[ii]==16)piontoe(H3_kin16,nkin16,kin[ii],target,pec,pecErr);
	 H3_pecErr[ii]=pecErr;
	 H3_pec[ii]=pec;

	outfile<<kin[ii]<<"  "<<H1_pec[ii]<<"  "<<H1_pecErr[ii]<<"  "<<D2_pec[ii]<<"  "<<D2_pecErr[ii]<<"  "
		<<He_pec[ii]<<"  "<<He_pecErr[ii]<<"  "<<H3_pec[ii]<<"  "<<H3_pecErr[ii]<<endl;
     }
 
     outfile.close();

 
}
