#include "SetCut.h"
#include "GetTrees.h"

void plot_Ion(){
     int nrun,kin;
     cout<<"Input run number: "<<endl;
     cin>>nrun;
     cout<<"Input kin: "<<endl;
     cin>>kin;

     TChain *T=new TChain("T");
     T=GetTree(nrun,kin,"T");

     Double_t nentries=T->GetEntries(CK+Ep+TRK+trigger2+beta+ACC+VZ);
     cout<<nentries<<endl;
/*
     TCanvas *c1=new TCanvas("c1");
     TH1F *h1=new TH1F("h1","EKLxe.Q2-EKLx.Q2",300,0,0.002);
     T->Draw("EKLxe.Q2-EKLx.Q2>>h1",CK+Ep+TRK+trigger2+beta+ACC+VZ);

     TCanvas *c2=new TCanvas("c2");
     TH1F *h2=new TH1F("h2","EKLxe.x_bj-EKLx.x_bj",300,0.,0.0002);
     T->Draw("EKLxe.x_bj-EKLx.x_bj>>h2",CK+Ep+TRK+trigger2+beta+ACC+VZ);
*/     

     TCanvas *c3=new TCanvas("c3");
     TH1F *hNu=new TH1F("hNu","nu distribution",300,6.5,8);
     T->Draw("EKLx.nu>>hNu",CK+Ep+TRK+trigger2+beta+ACC+VZ);

     TCanvas *c4=new TCanvas("c4");
     TH1F *hAngle=new TH1F("hAngle","angle distribution",300,12,20);
     T->Draw("EKLx.angle*180/3.14>>hAngle",CK+Ep+TRK+trigger2+beta+ACC+VZ);

     TCanvas *c5=new TCanvas("c5");
     TH2F *hnu_theta=new TH2F("hnu_theta","angle distribution",200,12,20,200,7,8);
     T->Draw("EKLx.nu:EKLx.angle*180/3.14>>hnu_theta",CK+Ep+TRK+trigger2+beta+ACC+VZ,"COLZ");

}
