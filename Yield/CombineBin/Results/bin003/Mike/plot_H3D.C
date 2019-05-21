#include <fstream>
using namespace std;

void plot_H3D()
{
   Double_t x1[22]={0.0},Ratio1[22]={0.0},Rerr1[22]={0.0};
   Double_t x2[22]={0.0},Ratio2[22]={0.0},Rerr2[22]={0.0};

   ifstream infile1;
   infile1.open("H3D_final.dat");
   Ssiz_t from=0;
   TString content,tmp;
   int nn=0;

   while(tmp.ReadLine(infile1)){
         if(nn==0){nn++;continue;}
         tmp.Tokenize(content,from," ");
         x1[nn-1]=atof(content.Data());
         tmp.Tokenize(content,from," ");
         Ratio1[nn-1]=atof(content.Data());
         tmp.Tokenize(content,from," ");
         Rerr1[nn-1]=atof(content.Data());
         from=0;
         nn++;
    }
   infile1.close();

   ifstream infile2;
   infile2.open("H3_D2_Hanjie.csv");
   from=0;
   nn=0;
   while(tmp.ReadLine(infile2)){
         if(nn==0){nn++;continue;}
         tmp.Tokenize(content,from,",");
         x2[nn-1]=atof(content.Data());
         tmp.Tokenize(content,from,",");
         Ratio2[nn-1]=atof(content.Data());
         tmp.Tokenize(content,from," ");
         Double_t rel=atof(content.Data());
         Rerr2[nn-1]=Ratio2[nn-1]*rel;
         from=0;
         nn++;
    }
   infile2.close();


   TGraphErrors *hratio1=new TGraphErrors();
   TGraphErrors *hratio2=new TGraphErrors();
   
   for(int ii=0;ii<22;ii++){
	hratio1->SetPoint(ii,x1[ii],Ratio1[ii]);
	hratio1->SetPointError(ii,0,Rerr1[ii]);
   } 

   for(int ii=0;ii<22;ii++){
	hratio2->SetPoint(ii,x2[ii],Ratio2[ii]);
	hratio2->SetPointError(ii,0,Rerr2[ii]);
   } 



   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   TMultiGraph *mg1=new TMultiGraph();
   hratio1->SetMarkerStyle(8);
   hratio1->SetMarkerColor(2);
   hratio2->SetMarkerStyle(8);
   hratio2->SetMarkerColor(1);
   mg1->Add(hratio1);
   mg1->Add(hratio2);
   mg1->Draw("AP");
   mg1->SetTitle("H3/D2 yield ratio;xbj;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hratio1,"hanjie","P");
   leg1->AddEntry(hratio2,"mike","P");
   leg1->Draw();

}
