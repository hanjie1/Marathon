#include <fstream>
using namespace std;

void plot()
{
     ifstream file1,file2;
     file1.open("D2_run3.txt");     
     file2.open("D2_kin3_run1327.dat");     

     Double_t xbj[8]={0.0},Yield1[8]={0.0},Yield2[8]={0.0},Yerr1[8]={0.0},Yerr2[8]={0.0};

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,", ");
          tmp.Tokenize(content,from,", ");
          xbj[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,", ");
          tmp.Tokenize(content,from,", ");
          Yield1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Yerr1[nn-1]=atof(content.Data());
          from=0;
     //     cout<<Yield1[nn-1]<<"  "<<Yerr1[nn-1]<<endl;
          nn++;
     }
     file1.close();
 
     from=0;
     nn=0;
     while(tmp.ReadLine(file2)){
          tmp.Tokenize(content,from,"	");
          xbj[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"	");
          Yield2[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"	");
          Yerr2[nn]=atof(content.Data());
          from=0;
      //    cout<<Yield2[nn]<<"  "<<Yerr2[nn]<<endl;
          nn++;
     }
    file2.close();

    TGraphErrors *g1=new TGraphErrors(8,xbj,Yield1,0,Yerr1);
    TGraphErrors *g2=new TGraphErrors(8,xbj,Yield2,0,Yerr2);
  
    TMultiGraph *mg=new TMultiGraph();
    g1->SetMarkerStyle(8);
    g1->SetMarkerColor(1);
    g2->SetMarkerStyle(8);
    g2->SetMarkerColor(2);

    mg->Add(g1);
    mg->Add(g2);
    mg->Draw("AP");
    mg->SetTitle("run 1327;");

    auto leg=new TLegend(0.65,0.7,0.85,0.85);
    leg->AddEntry(g1,"Hanjie","P");
    leg->AddEntry(g2,"Tyler","P");
    leg->Draw();
}
