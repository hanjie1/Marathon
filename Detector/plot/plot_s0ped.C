void plot_s0ped()
{
     ifstream file1;
     file1.open("../OUT/Ls0_ped.dat");

     Double_t nrun[1400]={0.0};
     Double_t s0la[1400]={0.0};
     Double_t s0ra[1400]={0.0};
     Double_t flag[1400]={0.0};

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0,total=0;
     while(tmp.ReadLine(file1)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,"  ");
          nrun[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s0la[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s0ra[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          flag[nn-1]=atof(content.Data());
          nn++;
          total++;
          from=0;
     }
     file1.close();

     TGraph *hs0la=new TGraph();
     TGraph *hs0ra=new TGraph();

     int mm=0;
     for(int ii=0;ii<total;ii++){
         if(flag[ii]==1)continue;
         hs0la->SetPoint(mm,nrun[ii],s0la[ii]);
         hs0ra->SetPoint(mm,nrun[ii],s0ra[ii]);
         mm++;
     }
     TCanvas *c1=new TCanvas("c1");
     hs0la->SetMarkerStyle(8);
     hs0la->SetMarkerColor(4);
     hs0la->Draw("AP");
     hs0la->SetTitle("LHRS s0la;run number;pedestal");

     TCanvas *c2=new TCanvas("c2");
     hs0ra->SetMarkerStyle(8);
     hs0ra->SetMarkerColor(4);
     hs0ra->Draw("AP");
     hs0ra->SetTitle("LHRS s0ra;run number;pedestal");

     c1->Print("Ls0_ped.pdf[");
     c1->Print("Ls0_ped.pdf");
     c2->Print("Ls0_ped.pdf");
     c2->Print("Ls0_ped.pdf]");






}


