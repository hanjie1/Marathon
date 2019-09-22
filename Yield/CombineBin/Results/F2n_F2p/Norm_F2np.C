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

void Norm_F2np(){
   int nfun=1;
   cout<<"Which function?  ";
   cin>>nfun;
   Double_t xup=0.4,xlo=0.2;
   cout<<"what is the upper range of x?  ";
   cin>>xup;
   cout<<"what is the lower range of x?  ";
   cin>>xlo;

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

    Double_t Dp_np[8]={0.0},Dp_npErr[8]={0.0};
    for(int ii=0;ii<8;ii++){
	Double_t D2_R=D2_EMC(Dp_x[ii]);
	Dp_np[ii]=Dp[ii]/D2_R-1.0;
	Dp_npErr[ii]=1.0/D2_R*Dp_err[ii];
    } 
    TGraphErrors *gDp=new TGraphErrors(8,Dp_x,Dp_np,0,Dp_npErr);

    Double_t FDp_p[3]={0.0};
    TF1 *fDp;
    Double_t np_fix=0.0;
    if(nfun==2){
       fDp=new TF1("fDp","pol2",xlo,xup);
       gDp->Fit("fDp","QR");
       fDp=gDp->GetFunction("fDp");
       FDp_p[0]=fDp->GetParameter(0);
       FDp_p[1]=fDp->GetParameter(1);
       FDp_p[2]=fDp->GetParameter(2);
       np_fix=FDp_p[0]+0.3*FDp_p[1]+0.3*0.3*FDp_p[2];
    }

    if(nfun==1){
       fDp=new TF1("fDp","pol1",xlo,xup);
       gDp->Fit("fDp","QR");
       fDp=gDp->GetFunction("fDp");
       FDp_p[0]=fDp->GetParameter(0);
       FDp_p[1]=fDp->GetParameter(1);
       np_fix=FDp_p[0]+0.3*FDp_p[1];
    }


    Double_t H3He_np[19]={0.0},H3He_npErr[19]={0.0};
    Double_t H3He_np_final[19]={0.0},H3He_npErr_final[19]={0.0};
    TH1F *hChi_H3He=new TH1F("hChi_H3He","chi square distribution",7,0,0.4);
    Double_t dmin=1.0;
    Double_t Nc_H3He=0.0;

    TGraphErrors *gfinal;

    for(int ii=0;ii<71;ii++){
	Double_t Nc=-0.035+ii*0.001;
	for(int jj=0;jj<19;jj++){
	    Double_t He_R=He_EMC(x[jj]);
	    Double_t H3_R=H3_EMC(x[jj]);
	    Double_t SR=H3_R/He_R;
	    Double_t tmp_H3He=(1.0-Nc)*H3He[jj];
	    Double_t tmp_H3He_err=(1.0-Nc)*H3He_err[jj];
            H3He_np[jj]=(2.0*tmp_H3He-SR)/(2.0*SR-tmp_H3He);
	    H3He_npErr[jj]=abs(3.0*SR/((2.0*SR-tmp_H3He)*(2.0*SR-tmp_H3He)))*tmp_H3He_err;
  	 } 
	    TGraphErrors *gH3He=new TGraphErrors(19,x,H3He_np,0,H3He_npErr);
	    Double_t np_H3He=0.0;
	    TF1 *fH3He;
	    if(nfun==2){
               gH3He->Fit("pol2","Q","",xlo,xup);
               fH3He=gH3He->GetFunction("pol2");
               Double_t FH3He_p0=fH3He->GetParameter(0);
               Double_t FH3He_p1=fH3He->GetParameter(1);
               Double_t FH3He_p2=fH3He->GetParameter(2);
               np_H3He=FH3He_p0+0.3*FH3He_p1+0.3*0.3*FH3He_p2;
	    }
	    if(nfun==1){
               gH3He->Fit("pol1","Q","",xlo,xup);
               fH3He=gH3He->GetFunction("pol1");
               Double_t FH3He_p0=fH3He->GetParameter(0);
               Double_t FH3He_p1=fH3He->GetParameter(1);
               np_H3He=FH3He_p0+0.3*FH3He_p1;
	    }
	    Double_t chi2=fH3He->GetChisquare()/fH3He->GetNDF();
	    hChi_H3He->Fill(chi2);
	    cout<<Nc<<"  "<<abs(np_H3He-np_fix)<<endl;
	    if(abs(np_H3He-np_fix)<dmin){
		gfinal=(TGraphErrors*)gH3He->Clone();
		Nc_H3He=Nc;
		dmin=abs(np_H3He-np_fix);
	    }
	
    }
    cout<<Nc_H3He<<endl;

    TCanvas *c1=new TCanvas("c1","c1",1500,1000);
    TMultiGraph *mg=new TMultiGraph();
    gfinal->SetMarkerStyle(8);
    gfinal->SetMarkerColor(4);
    gfinal->SetMarkerSize(2);
    gDp->SetMarkerStyle(22);
    gDp->SetMarkerColor(2);
    gDp->SetMarkerSize(2);
    mg->Add(gfinal);
    mg->Add(gDp);
    mg->Draw("AP");
    c1->Print("fit.pdf");
}
