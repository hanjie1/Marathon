#include "ReadFile.h"
#include "libEMCR.h"
#include "libNP.h"
void Norm_KP(){
   Double_t xmin=0.0,xmax=1.0;
   cout<<"what is the minimum x?   ";
   cin>>xmin;
   cout<<"what is the maximum x?   ";
   cin>>xmax;

   Double_t x[19]={0.0},H3D[19]={0.0},H3D_err[19]={0.0};
   Double_t HeD[19]={0.0},HeD_err[19]={0.0};
   Double_t H3He[19]={0.0},H3He_err[19]={0.0};
   Double_t Dp_x[8]={0.0},Dp[8]={0.0},Dp_err[8]={0.0};

   TString Rfile="newbin/Dp_final.dat";
   int nbin1=ReadYield(Rfile,Dp_x,Dp,Dp_err);
   Rfile="newbin/HeD_final.dat";
   int nbin2=ReadYield(Rfile,x,HeD,HeD_err);
   Rfile="newbin/H3D_final.dat";
   nbin2=ReadYield(Rfile,x,H3D,H3D_err);
   Rfile="newbin/H3He_final.dat";
   nbin2=ReadYield(Rfile,x,H3He,H3He_err);

    Double_t H3He_np[19]={0.0},H3He_npErr[19]={0.0};
    Double_t chi2_min=100.0;
    Double_t Nc_H3He=0.0;

    TGraphErrors *gfinal_H3He;
    for(int ii=0;ii<71;ii++){
	Double_t Nc=-0.035+ii*0.001;
	Double_t tmpChi2=0.0;
	for(int jj=0;jj<19;jj++){
	    Double_t He_R=He_EMC(x[jj]);
	    Double_t H3_R=H3_EMC(x[jj]);
	    Double_t SR=H3_R/He_R;
	    Double_t tmp_H3He=(1.0-Nc)*H3He[jj];
	    Double_t tmp_H3He_err=(1.0-Nc)*H3He_err[jj];
            H3He_np[jj]=(2.0*tmp_H3He-SR)/(2.0*SR-tmp_H3He);
	    H3He_npErr[jj]=abs(3.0*SR/((2.0*SR-tmp_H3He)*(2.0*SR-tmp_H3He)))*tmp_H3He_err;
	
	    if(x[jj]<xmin || x[jj]>xmax)continue;
	    Double_t tmpKP=KP_NP(x[jj]);
	    tmpChi2=tmpChi2+pow((H3He_np[jj]-tmpKP)/H3He_npErr[jj],2);
  	 } 
	    TGraphErrors *gH3He=new TGraphErrors(19,x,H3He_np,0,H3He_npErr);
	    Double_t np_H3He=0.0;

	    if(tmpChi2<chi2_min){
		gfinal_H3He=(TGraphErrors*)gH3He->Clone();
		Nc_H3He=Nc;
		chi2_min=tmpChi2;
	    }
	
    }

    ofstream outfile1;
    outfile1.open("H3D_Norm_err.dat");
    Double_t H3D_np[19]={0.0},H3D_npErr[19]={0.0};
    chi2_min=100.0;
    Double_t Nc_H3D=0.0;
    TGraphErrors *gfinal_H3D;
    for(int ii=0;ii<71;ii++){
	Double_t Nc=-0.035+ii*0.001;
	Double_t tmpChi2=0.0;
        outfile1<<"-------- "<<Nc<<"-------"<<endl;
	for(int jj=0;jj<19;jj++){
	    Double_t D2_R=D2_EMC(x[jj]);
	    Double_t H3_R=H3_EMC(x[jj]);
	    Double_t SR=H3_R/D2_R;
	    Double_t tmp_H3D=(1.0-Nc)*H3D[jj];
	    Double_t tmp_H3D_err=(1.0-Nc)*H3D_err[jj];
            H3D_np[jj]=(tmp_H3D-SR)/(2.0*SR-tmp_H3D);
            H3D_npErr[jj]=abs(SR/((2.0*SR-tmp_H3D)*(2.0*SR-tmp_H3D)))*tmp_H3D_err;
	
	    if(x[jj]<xmin || x[jj]>xmax)continue;
	    Double_t tmpKP=KP_NP(x[jj]);
	    tmpChi2=tmpChi2+pow((H3D_np[jj]-tmpKP)/H3D_npErr[jj],2);

	    Double_t tmpER31=H3EMC_ERROR(x[jj]);
	    Double_t tmpER21=D2EMC_ERROR(x[jj]);
	    Double_t tmpErr=SR*sqrt(tmpER31*tmpER31+tmpER21*tmpER21);
	    Double_t relKP=tmp_H3D/pow((2.0*SR-tmp_H3D),2)*tmpErr/H3D_np[jj];
	    outfile1<<x[jj]<<"  "<<relKP<<"  "<<tmpErr<<"  "<<H3D_npErr[jj]/H3D_np[jj]<<endl;
  	 } 
	    TGraphErrors *gH3D=new TGraphErrors(19,x,H3D_np,0,H3D_npErr);
	    if(tmpChi2<chi2_min){
		gfinal_H3D=(TGraphErrors*)gH3D->Clone();
		Nc_H3D=Nc;
		chi2_min=tmpChi2;
	    }
    }
    outfile1.close();

    Double_t HeD_np[19]={0.0},HeD_npErr[19]={0.0};
    chi2_min=100.0;
    Double_t Nc_HeD=0.0;
    TGraphErrors *gfinal_HeD;
    for(int ii=0;ii<71;ii++){
	Double_t Nc=-0.035+ii*0.001;
	Double_t tmpChi2=0.0;
	for(int jj=0;jj<19;jj++){
	    Double_t D2_R=D2_EMC(x[jj]);
	    Double_t He_R=He_EMC(x[jj]);
	    Double_t SR=He_R/D2_R;
	    Double_t tmp_HeD=(1.0-Nc)*HeD[jj];
	    Double_t tmp_HeD_err=(1.0-Nc)*HeD_err[jj];
            HeD_np[jj]=(tmp_HeD-2.0*SR)/(SR-tmp_HeD);
            HeD_npErr[jj]=abs(SR/((SR-tmp_HeD)*(SR-tmp_HeD)))*tmp_HeD_err;

	    if(x[jj]<xmin || x[jj]>xmax)continue;
	    Double_t tmpKP=KP_NP(x[jj]);
	    tmpChi2=tmpChi2+pow((HeD_np[jj]-tmpKP)/HeD_npErr[jj],2);
  	 } 
	    TGraphErrors *gHeD=new TGraphErrors(19,x,HeD_np,0,HeD_npErr);
	    if(tmpChi2<chi2_min){
		gfinal_HeD=(TGraphErrors*)gHeD->Clone();
		Nc_HeD=Nc;
		chi2_min=tmpChi2;
	    }
    }

    ofstream outfile3;
    outfile3.open("Dp_Norm_err.dat");
    Double_t Dp_np[8]={0.0},Dp_npErr[8]={0.0};
    chi2_min=100.0;
    Double_t Nc_Dp=0.0;
    TGraphErrors *gfinal_Dp;
    for(int ii=0;ii<71;ii++){
	Double_t Nc=-0.035+ii*0.001;
	Double_t tmpChi2=0.0;
	outfile3<<"-----------"<<Nc<<"----------"<<endl;
	for(int jj=0;jj<8;jj++){
	    Double_t D2_R=D2_EMC(Dp_x[jj]);
	    Double_t tmp_Dp=(1.0-Nc)*Dp[jj];
	    Double_t tmp_Dp_err=(1.0-Nc)*Dp_err[jj];
            Dp_np[jj]=tmp_Dp/D2_R-1.0;
            Dp_npErr[jj]=abs(1.0/D2_R)*tmp_Dp_err;

	    if(Dp_x[jj]<xmin || Dp_x[jj]>xmax)continue;
	    Double_t tmpKP=KP_NP(Dp_x[jj]);
	    tmpChi2=tmpChi2+pow((Dp_np[jj]-tmpKP)/Dp_npErr[jj],2);

	    Double_t tmpER21=D2EMC_ERROR(x[jj])*D2_R;
	    Double_t relKP=tmp_Dp/pow(D2_R,2)*tmpER21/Dp_np[jj];
	    outfile3<<x[jj]<<"  "<<relKP<<"  "<<tmpER21<<"  "<<Dp_npErr[jj]/Dp_np[jj]<<endl;
  	 } 
	    TGraphErrors *gDp=new TGraphErrors(8,Dp_x,Dp_np,0,Dp_npErr);
	    if(tmpChi2<chi2_min){
		gfinal_Dp=(TGraphErrors*)gDp->Clone();
		Nc_Dp=Nc;
		chi2_min=tmpChi2;
	    }
    }
    outfile3.close();

    TGraph *gKP=new TGraph();
    for(int ii=0;ii<19;ii++){
        gKP->SetPoint(ii,x[ii],KP_NP(x[ii]));
    }

    cout<<"H3He:   "<<Nc_H3He<<endl;
    cout<<"H3D:   "<<Nc_H3D<<endl;
    cout<<"HeD:   "<<Nc_HeD<<endl;
    cout<<"Dp:   "<<Nc_Dp<<endl;

    TCanvas *c1=new TCanvas("c1","c1",1500,1000);
    TMultiGraph *mg=new TMultiGraph();
    gfinal_Dp->SetMarkerColor(4);
    gfinal_Dp->SetMarkerStyle(24);
    gfinal_Dp->SetMarkerSize(2);
    gfinal_Dp->SetLineColor(4);
    gfinal_H3He->SetMarkerColor(2);
    gfinal_H3He->SetMarkerStyle(26);
    gfinal_H3He->SetMarkerSize(2);
    gfinal_H3He->SetLineColor(2);
    gfinal_HeD->SetMarkerColor(kViolet-1);
    gfinal_HeD->SetMarkerStyle(32);
    gfinal_HeD->SetMarkerSize(2);
    gfinal_HeD->SetLineColor(kViolet-1);
    gfinal_H3D->SetMarkerColor(1);
    gfinal_H3D->SetMarkerStyle(25);
    gfinal_H3D->SetMarkerSize(2);
    gfinal_H3D->SetLineColor(1);
    gKP->SetLineStyle(1);
    gKP->SetLineColor(8);
    mg->Add(gfinal_H3He,"P");
    mg->Add(gfinal_H3D,"P");
    mg->Add(gfinal_HeD,"P");
    mg->Add(gfinal_Dp,"P");
    mg->Add(gKP,"L");
    mg->Draw("A");
    mg->SetTitle(";Bjorken x;F_{2}^{n} / F_{2}^{p}");
    mg->GetXaxis()->SetRangeUser(0.2,0.4);
    mg->GetYaxis()->SetRangeUser(0.5,0.8);

   auto leg1=new TLegend(0.65,0.65,0.88,0.88);
   leg1->AddEntry(gfinal_Dp,"#scale[0.8]{F_{2}^{^{2}H} / F_{2}^{^{1}H} (+0.6%)}","P");
   leg1->AddEntry(gfinal_H3He,"#scale[0.8]{F_{2}^{^{3}H} / F_{2}^{^{3}He} (-2.6%)}","P");
   leg1->AddEntry(gfinal_H3D,"#scale[0.8]{F_{2}^{^{3}H} / F_{2}^{^{2}H} (-0.2%)}","P");
   leg1->AddEntry(gfinal_HeD,"#scale[0.8]{F_{2}^{^{3}He} / F_{2}^{^{2}H} (+2.4%)}","P");
   leg1->AddEntry(gKP,"#scale[0.7]{KP model}","L");
   leg1->SetMargin(0.4);
   leg1->Draw();


    c1->Print("fit.pdf");
}
