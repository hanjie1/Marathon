#include "ReadFile.h"
void FitKPerr(){
     Double_t X1[MAXBIN]={0.0},Y1[MAXBIN]={0.0};
     Double_t X2[MAXBIN]={0.0},Y2[MAXBIN]={0.0};
     Double_t X3[MAXBIN]={0.0},Y3[MAXBIN]={0.0};

     TString Rfile;
     Rfile="R21.csv";
     int nbin1=ReadYield(Rfile,X1,Y1);
     Rfile="R31.csv";
     int nbin2=ReadYield(Rfile,X2,Y2);
     Rfile="R32.csv";
     int nbin3=ReadYield(Rfile,X3,Y3);

     TGraph *gR21=new TGraph(nbin1,X1,Y1);
     TGraph *gR31=new TGraph(nbin2,X2,Y2);
     TGraph *gR32=new TGraph(nbin3,X3,Y3);

     TCanvas *c1=new TCanvas("c1","c1",1500,1200);
     gR21->Draw("AP*");
     gR21->SetTitle("R21;x;");

     TCanvas *c2=new TCanvas("c2","c2",1500,1200);
     gR31->Draw("AP*");
     gR31->SetTitle("R31;x;");

     TCanvas *c3=new TCanvas("c3","c3",1500,1200);
     gR32->Draw("AP*");
     gR32->SetTitle("R32;x;");
}
