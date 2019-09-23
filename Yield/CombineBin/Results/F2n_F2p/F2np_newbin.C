#include "ReadFile.h"
Double_t D2_EMC(Double_t xbj){
	 Double_t R=1.0;
	 R=0.959252+1.70218*xbj-21.3065*pow(xbj,2)+141.459*pow(xbj,3)-568.77*pow(xbj,4)+1443.77*pow(xbj,5)
	   -2338.81*pow(xbj,6)+2353.76*pow(xbj,7)-1344.9*pow(xbj,8)+334.502*pow(xbj,9);

	 return R;
}

Double_t He_EMC(Double_t x){
	Double_t R=1.0;
	R=0.949891+2.38613*x-30.205*pow(x,2)+200.737*pow(x,3)-807.611*pow(x,4)+2050.04*pow(x,5)-3319.42*pow(x,6)
	  +3338.26*pow(x,7)-1906.05*pow(x,8)+473.934*pow(x,9);	
	return R;
}

Double_t H3_EMC(Double_t x){
	Double_t R=1.0;
	R=0.939068+2.78422*x-35.2978*pow(x,2)+234.697*pow(x,3)-940.068*pow(x,4)+2368.52*pow(x,5)-3802.3*pow(x,6)
	  +3791.3*pow(x,7)-2147.29*pow(x,8)+529.779*pow(x,9);

	return R;
}

int F2np_newbin(){
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
    Double_t H3D_np[19]={0.0},H3D_npErr[19]={0.0};
    Double_t HeD_np[19]={0.0},HeD_npErr[19]={0.0};
    ofstream outfile1;
    outfile1.open("results/F2np.dat");
    outfile1<<"x    H3He   H3He_err    HeD    HeD_err    H3D    H3D_err "<<endl;
    for(int ii=0;ii<19;ii++){
	Double_t He_R=He_EMC(x[ii]);
	Double_t H3_R=H3_EMC(x[ii]);
	Double_t D2_R=D2_EMC(x[ii]);

	Double_t SR=H3_R/He_R;
	Double_t SR_HeD=He_R/D2_R;
	Double_t SR_H3D=H3_R/D2_R;

	H3He_np[ii]=(2.0*H3He[ii]-SR)/(2.0*SR-H3He[ii]);
	H3He_npErr[ii]=abs(3.0*SR/((2.0*SR-H3He[ii])*(2.0*SR-H3He[ii])))*H3He_err[ii];

	HeD_np[ii]=(HeD[ii]-2.0*SR_HeD)/(SR_HeD-HeD[ii]);
	HeD_npErr[ii]=abs(SR_HeD/((SR_HeD-HeD[ii])*(SR_HeD-HeD[ii])))*HeD_err[ii];

	H3D_np[ii]=(H3D[ii]-SR_H3D)/(2.0*SR_H3D-H3D[ii]);
	H3D_npErr[ii]=abs(SR_H3D/((2.0*SR_H3D-H3D[ii])*(2.0*SR_H3D-H3D[ii])))*H3D_err[ii];
	outfile1<<x[ii]<<"  "<<H3He_np[ii]<<"  "<<H3He_npErr[ii]<<"  "<<HeD_np[ii]<<"  "<<HeD_npErr[ii]<<"  "
		<<H3D_np[ii]<<"  "<<H3D_npErr[ii]<<endl;
    }
    outfile1.close();

    ofstream outfile2;
    outfile2.open("results/F2np_Dp.dat");
    outfile2<<"x    Dp    Dp_err "<<endl;

    Double_t Dp_np[8]={0.0},Dp_npErr[8]={0.0};
    for(int ii=0;ii<8;ii++){
	Double_t D2_R=D2_EMC(Dp_x[ii]);
	Dp_np[ii]=Dp[ii]/D2_R-1.0;
	Dp_npErr[ii]=1.0/D2_R*Dp_err[ii];
	outfile2<<Dp_x[ii]<<"  "<<Dp_np[ii]<<"  "<<Dp_npErr[ii]<<endl;
    } 
    outfile2.close();

    TGraphErrors *gH3He=new TGraphErrors(19,x,H3He_np,0,H3He_npErr);
    TGraphErrors *gH3D=new TGraphErrors(19,x,H3D_np,0,H3D_npErr);
    TGraphErrors *gHeD=new TGraphErrors(19,x,HeD_np,0,HeD_npErr);
    TGraphErrors *gDp=new TGraphErrors(8,Dp_x,Dp_np,0,Dp_npErr);

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    TMultiGraph *mg1=new TMultiGraph();
    gDp->SetMarkerColor(8);
    gDp->SetMarkerStyle(8);
    gDp->SetMarkerSize(1.7);
    gH3He->SetMarkerColor(2);
    gH3He->SetMarkerStyle(22);
    gH3He->SetMarkerSize(1.7);
    gHeD->SetMarkerColor(4);
    gHeD->SetMarkerStyle(33);
    gHeD->SetMarkerSize(1.7);
    gH3D->SetMarkerColor(1);
    gH3D->SetMarkerStyle(29);
    gH3D->SetMarkerSize(1.7);
    mg1->Add(gDp);
    mg1->Add(gH3He);
//    mg1->Add(gH3D);
//    mg1->Add(gHeD);
    mg1->Draw("AP"); 
    mg1->GetYaxis()->SetRangeUser(0.5,0.85);
  
    return 0;
}
