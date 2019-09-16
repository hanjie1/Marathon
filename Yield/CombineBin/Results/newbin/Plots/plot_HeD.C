#include <fstream>
#include "ReadFile.h"
using namespace std;

Double_t KP_EMC(Double_t x){
	Double_t ratio=1.0;
	ratio=1.02967-0.135929*x+3.92009*pow(x,2)-21.2861*pow(x,3)+64.7762*pow(x,4)
	      -129.928*pow(x,5)+169.609*pow(x,6)-127.386*pow(x,7)+41.0723*pow(x,8);
	ratio=ratio*3.0/2.0;
	return ratio;
}

int ReadHallC(Double_t x[],Double_t Ratio[],Double_t Rerr[]){
     ifstream infile;
     infile.open("/w/halla-scifs17exp/triton/hanjie/MARATHON/analysis/Yield/CombineBin/Results/Other_Data/HeD_HallC.dat");
     if(!infile.is_open())return 0;

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;

     while(tmp.ReadLine(infile)){
          if(nn==0){nn++;continue;}
          if(tmp[0]=='#')continue;
          tmp.Tokenize(content,from," ");
          x[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          Double_t tmpR=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Double_t tmpE_stat=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Double_t tmpE_sys=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Double_t tmp_Iso=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Double_t tmp_Coul=atof(content.Data());

	  Ratio[nn-1]=tmpR/(tmp_Iso*tmp_Coul)*3.0/2.0;
	  Rerr[nn-1]=sqrt(tmpE_stat*tmpE_stat+tmpE_sys*tmpE_sys)/(tmp_Iso*tmp_Coul)*3.0/2.0;

          from=0;
          nn++;
     }
    infile.close();

    nn=nn-1;
    return nn;
}

void plot_HeD()
{
   Double_t x[MAXBIN]={0.0},Ratio[MAXBIN]={0.0},Rerr[MAXBIN]={0.0};
   Double_t xC[MAXBIN]={0.0},RatioC[MAXBIN]={0.0},RerrC[MAXBIN]={0.0};
   
   int nbinC=ReadHallC(xC,RatioC,RerrC);
   auto f1=new TF1("f1","KP_EMC(x)",0.16,0.84);

   TString Rfile="newbin/HeD_final.dat";
   int nbin=ReadYield(Rfile,x,Ratio,Rerr); 

   TGraphErrors *hratio=new TGraphErrors();
   TGraphErrors *hratioC=new TGraphErrors();
   
   for(int ii=0;ii<nbin;ii++){
	hratio->SetPoint(ii,x[ii],Ratio[ii]);
	hratio->SetPointError(ii,0,Rerr[ii]);
   } 

   for(int ii=0;ii<nbinC;ii++){
	hratioC->SetPoint(ii,xC[ii],RatioC[ii]);
	hratioC->SetPointError(ii,0,RerrC[ii]);
   } 

   TCanvas *c1=new TCanvas("c1","c1",1500,1200);
   TMultiGraph *mg1=new TMultiGraph();
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(2);
   hratioC->SetMarkerStyle(8);
   hratioC->SetMarkerColor(4);
   mg1->Add(hratio);
   mg1->Add(hratioC);
   mg1->Draw("AP");
   mg1->SetTitle(";Bjorken x;#sigma({}^{3}He)/#sigma({}^{2}H)");
   mg1->GetYaxis()->SetRangeUser(1.5,1.75);

   f1->SetLineColor(4);
   f1->SetLineStyle(9);
   f1->Draw("same");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio,"MARATHON","P");
   leg1->AddEntry(hratioC,"Hall C","P");
   leg1->Draw();


   c1->Print("HeD_final.pdf");
}
