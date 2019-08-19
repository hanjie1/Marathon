#include "GetTrees.h"
#include "GetRunList.h"
struct EFF_VAR{
      Double_t Notrack_eff;
      Double_t Notrack_err;
      Double_t Onetrack_eff;
      Double_t Onetrack_err;
      Double_t Mtrack_eff;
      Double_t Mtrack_err;
};

void track(int run_number, int kin, int HRS, EFF_VAR &TT){
     TString TreeName="T";
     TChain *T=GetTree(run_number,kin,TreeName);

     //HRS=0 is LHRS, HRS=1 is RHRS
     Double_t eff_no=0.0,eff_noerr=0.0,eff_mt=0.0,eff_mterr=0.0;
     Double_t eff_1=0.0,eff_1err=0.0;
     if(HRS==0){
	TCut T2="(DL.evtypebits>>2)&1";
	TCut CK="L.cer.asum_c>2200 && L.cer.asum_c<10000";
	TCut EE="(L.prl1.e+L.prl2.e)/3100.0>0.75";
        TCut PS="L.prl1.e>(-16.0/15.0*L.prl2.e+1600.0)";
	TCut notrack="L.tr.n==0";
	TCut onetrack="L.tr.n==1";
	TCut mtrack="L.tr.n>1";
	
	Double_t ntrack=T->GetEntries(T2+CK+EE);
	Double_t ntrack0=T->GetEntries(T2+CK+EE+notrack);
	Double_t ntrack1=T->GetEntries(T2+CK+EE+onetrack);
	Double_t ntrackM=T->GetEntries(T2+CK+EE+mtrack);

	eff_no=ntrack0/ntrack;
	eff_noerr=eff_no*sqrt(1.0/ntrack0-1.0/ntrack);

	eff_1=ntrack1/ntrack;
	eff_1err=eff_1*sqrt(1.0/ntrack1-1.0/ntrack);

	eff_mt=ntrackM/ntrack;
	eff_mterr=eff_mt*sqrt(1.0/ntrackM-1.0/ntrack);
     }

     if(HRS==1){
	TCut T2="(DR.evtypebits>>5)&1";
	TCut CK="R.cer.asum_c>3000 && R.cer.asum_c<10000";
	TCut EE="(R.ps.e+R.sh.e)/2900.0>0.75";
	TCut notrack="R.tr.n==0";
	TCut onetrack="R.tr.n==1";
	TCut mtrack="R.tr.n>1";
	
	Double_t ntrack=T->GetEntries(T2+CK+EE);
	Double_t ntrack0=T->GetEntries(T2+CK+EE+notrack);
	Double_t ntrack1=T->GetEntries(T2+CK+EE+onetrack);
	Double_t ntrackM=T->GetEntries(T2+CK+EE+mtrack);

	eff_no=ntrack0/ntrack;
	eff_noerr=eff_no*sqrt(1.0/ntrack0-1.0/ntrack);

	eff_1=ntrack1/ntrack;
	eff_1err=eff_1*sqrt(1.0/ntrack1-1.0/ntrack);

	eff_mt=ntrackM/ntrack;
	eff_mterr=eff_mt*sqrt(1.0/ntrackM-1.0/ntrack);
     }

     TT.Notrack_eff=eff_no;
     TT.Notrack_err=eff_noerr;
     TT.Onetrack_eff=eff_1;
     TT.Onetrack_err=eff_1err;
     TT.Mtrack_eff=eff_mt;
     TT.Mtrack_err=eff_mterr;

     return;
}


void Track_eff(){
     TString target[4]={"H1","D2","He3","H3"};
     int kin[4]={0,1,2,1};

     int H1_nkin[4]={2476,1213,1240,90849};
     int D2_nkin[4]={2492,1214,1243,90877};
     int He_nkin[4]={2528,1210,1250,90855};
     int H3_nkin[4]={2501,1216,1246,90875};

     Double_t H1_tr0[4]={0.0},H1_tr0err[4]={0.0};
     Double_t D2_tr0[4]={0.0},D2_tr0err[4]={0.0};
     Double_t He_tr0[4]={0.0},He_tr0err[4]={0.0};
     Double_t H3_tr0[4]={0.0},H3_tr0err[4]={0.0};

     TGraphErrors *gH1=new TGraphErrors();
     TGraphErrors *gD2=new TGraphErrors();
     TGraphErrors *gHe=new TGraphErrors();
     TGraphErrors *gH3=new TGraphErrors();

     TGraphErrors *gH1_R=new TGraphErrors();
     TGraphErrors *gD2_R=new TGraphErrors();
     TGraphErrors *gHe_R=new TGraphErrors();
     TGraphErrors *gH3_R=new TGraphErrors();

     for(int ii=0;ii<4;ii++){
   	 EFF_VAR TT;
         int HRS=0;
	 if(ii==3)HRS=1;
	 track(H1_nkin[ii],kin[ii],HRS,TT);
	 H1_tr0[ii]=TT.Onetrack_eff;
	 H1_tr0err[ii]=TT.Onetrack_err;
  	 if(HRS==0){
	    gH1->SetPoint(ii,kin[ii],H1_tr0[ii]);
	    gH1->SetPointError(ii,0,H1_tr0err[ii]);
	 }
  	 else{
	    gH1_R->SetPoint(0,kin[ii],H1_tr0[ii]);
	    gH1_R->SetPointError(0,0,H1_tr0err[ii]);
	 }

	 track(D2_nkin[ii],kin[ii],HRS,TT);
	 D2_tr0[ii]=TT.Onetrack_eff;
	 D2_tr0err[ii]=TT.Onetrack_err;
  	 if(HRS==0){
	    gD2->SetPoint(ii,kin[ii],D2_tr0[ii]);
	    gD2->SetPointError(ii,0,D2_tr0err[ii]);
	 }
  	 else{
	    gD2_R->SetPoint(0,kin[ii],D2_tr0[ii]);
	    gD2_R->SetPointError(0,0,D2_tr0err[ii]);
	 }

	 track(He_nkin[ii],kin[ii],HRS,TT);
	 He_tr0[ii]=TT.Onetrack_eff;
	 He_tr0err[ii]=TT.Onetrack_err;
  	 if(HRS==0){
	    gHe->SetPoint(ii,kin[ii],He_tr0[ii]);
	    gHe->SetPointError(ii,0,He_tr0err[ii]);
	 }
  	 else{
	    gHe_R->SetPoint(0,kin[ii],He_tr0[ii]);
	    gHe_R->SetPointError(0,0,He_tr0err[ii]);
	 }


	 track(H3_nkin[ii],kin[ii],HRS,TT);
	 H3_tr0[ii]=TT.Onetrack_eff;
	 H3_tr0err[ii]=TT.Onetrack_err;
  	 if(HRS==0){
	    gH3->SetPoint(ii,kin[ii],H3_tr0[ii]);
	    gH3->SetPointError(ii,0,H3_tr0err[ii]);
	 }
  	 else{
	    gH3_R->SetPoint(0,kin[ii],H3_tr0[ii]);
	    gH3_R->SetPointError(0,0,H3_tr0err[ii]);
	 }


	 cout<<"kin"<<kin[ii]<<":  "<<H1_tr0[ii]<<"  "<<H1_tr0err[ii]<<endl;
	 cout<<"kin"<<kin[ii]<<":  "<<D2_tr0[ii]<<"  "<<D2_tr0err[ii]<<endl;
	 cout<<"kin"<<kin[ii]<<":  "<<He_tr0[ii]<<"  "<<He_tr0err[ii]<<endl;
	 cout<<"kin"<<kin[ii]<<":  "<<H3_tr0[ii]<<"  "<<H3_tr0err[ii]<<endl;
     }

     gH1->SetMarkerColor(1);
     gH1->SetMarkerStyle(8);
     gD2->SetMarkerColor(4);
     gD2->SetMarkerStyle(8);
     gHe->SetMarkerColor(8);
     gHe->SetMarkerStyle(8);
     gH3->SetMarkerColor(2);
     gH3->SetMarkerStyle(8);

     gH1_R->SetMarkerColor(1);
     gH1_R->SetMarkerStyle(22);
     gD2_R->SetMarkerColor(4);
     gD2_R->SetMarkerStyle(22);
     gHe_R->SetMarkerColor(8);
     gHe_R->SetMarkerStyle(22);
     gH3_R->SetMarkerColor(2);
     gH3_R->SetMarkerStyle(22);

     TMultiGraph *mg=new TMultiGraph();
     mg->Add(gH1);
     mg->Add(gD2);
     mg->Add(gHe);
     mg->Add(gH3);
     mg->Add(gH1_R);
     mg->Add(gD2_R);
     mg->Add(gHe_R);
     mg->Add(gH3_R);
     mg->Draw("AP"); 

     

}

