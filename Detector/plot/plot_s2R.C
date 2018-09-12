void plot_s2R()
{
     ifstream file1;
     file1.open("../OUT/Ls2R_peak.dat");

     Double_t nrun[1400]={0.0};
     Double_t s2R1[1400]={0.0};
     Double_t s2R2[1400]={0.0};
     Double_t s2R3[1400]={0.0};
     Double_t s2R4[1400]={0.0};
     Double_t s2R5[1400]={0.0};
     Double_t s2R6[1400]={0.0};
     Double_t s2R7[1400]={0.0};
     Double_t s2R8[1400]={0.0};
     Double_t s2R9[1400]={0.0};
     Double_t s2R10[1400]={0.0};
     Double_t s2R11[1400]={0.0};
     Double_t s2R12[1400]={0.0};
     Double_t s2R13[1400]={0.0};
     Double_t s2R14[1400]={0.0};
     Double_t s2R15[1400]={0.0};
     Double_t s2R16[1400]={0.0};
     Double_t flag[1400]={0.0};

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0,total=0;
     while(tmp.ReadLine(file1)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,"  ");
          nrun[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R2[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R3[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R4[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R5[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R6[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R7[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R8[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R9[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R10[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R11[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R12[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R13[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R14[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R15[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2R16[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          flag[nn-1]=atof(content.Data());
          nn++;
          total++;
          from=0;
     }
     file1.close();

     TGraph *hs2R0=new TGraph();
     TGraph *hs2R1=new TGraph();
     TGraph *hs2R2=new TGraph();
     TGraph *hs2R3=new TGraph();
     TGraph *hs2R4=new TGraph();
     TGraph *hs2R5=new TGraph();
     TGraph *hs2R6=new TGraph();
     TGraph *hs2R7=new TGraph();
     TGraph *hs2R8=new TGraph();
     TGraph *hs2R9=new TGraph();
     TGraph *hs2R10=new TGraph();
     TGraph *hs2R11=new TGraph();
     TGraph *hs2R12=new TGraph();
     TGraph *hs2R13=new TGraph();
     TGraph *hs2R14=new TGraph();
     TGraph *hs2R15=new TGraph();

     int mm=0;
     for(int ii=0;ii<total;ii++){
         if(flag[ii]==1)continue;
         hs2R0->SetPoint(mm,nrun[ii],s2R1[ii]);
         hs2R1->SetPoint(mm,nrun[ii],s2R2[ii]);
         hs2R2->SetPoint(mm,nrun[ii],s2R3[ii]);
         hs2R3->SetPoint(mm,nrun[ii],s2R4[ii]);
         hs2R4->SetPoint(mm,nrun[ii],s2R5[ii]);
         hs2R5->SetPoint(mm,nrun[ii],s2R6[ii]);
         hs2R6->SetPoint(mm,nrun[ii],s2R7[ii]);
         hs2R7->SetPoint(mm,nrun[ii],s2R8[ii]);
         hs2R8->SetPoint(mm,nrun[ii],s2R9[ii]);
         hs2R9->SetPoint(mm,nrun[ii],s2R10[ii]);
         hs2R10->SetPoint(mm,nrun[ii],s2R11[ii]);
         hs2R11->SetPoint(mm,nrun[ii],s2R12[ii]);
         hs2R12->SetPoint(mm,nrun[ii],s2R13[ii]);
         hs2R13->SetPoint(mm,nrun[ii],s2R14[ii]);
         hs2R14->SetPoint(mm,nrun[ii],s2R15[ii]);
         hs2R15->SetPoint(mm,nrun[ii],s2R16[ii]);
         mm++;
     }
     TCanvas *c1=new TCanvas("c1");
     hs2R0->SetMarkerStyle(8);
     hs2R0->SetMarkerColor(4);
     hs2R0->Draw("AP");
     hs2R0->SetTitle("LHRS s2R1;run number;peak");

     TCanvas *c2=new TCanvas("c2");
     hs2R1->SetMarkerStyle(8);
     hs2R1->SetMarkerColor(4);
     hs2R1->Draw("AP");
     hs2R1->SetTitle("LHRS s2R2;run number;peak");

     TCanvas *c3=new TCanvas("c3");
     hs2R2->SetMarkerStyle(8);
     hs2R2->SetMarkerColor(4);
     hs2R2->Draw("AP");
     hs2R2->SetTitle("LHRS s2R3;run number;peak");

     TCanvas *c4=new TCanvas("c4");
     hs2R3->SetMarkerStyle(8);
     hs2R3->SetMarkerColor(4);
     hs2R3->Draw("AP");
     hs2R3->SetTitle("LHRS s2R4;run number;peak");

     TCanvas *c5=new TCanvas("c5");
     hs2R4->SetMarkerStyle(8);
     hs2R4->SetMarkerColor(4);
     hs2R4->Draw("AP");
     hs2R4->SetTitle("LHRS s2R5;run number;peak");

     TCanvas *c6=new TCanvas("c6");
     hs2R5->SetMarkerStyle(8);
     hs2R5->SetMarkerColor(4);
     hs2R5->Draw("AP");
     hs2R5->SetTitle("LHRS s2R6;run number;peak");

     TCanvas *c7=new TCanvas("c7");
     hs2R6->SetMarkerStyle(8);
     hs2R6->SetMarkerColor(4);
     hs2R6->Draw("AP");
     hs2R6->SetTitle("LHRS s2R7;run number;peak");

     TCanvas *c8=new TCanvas("c8");
     hs2R7->SetMarkerStyle(8);
     hs2R7->SetMarkerColor(4);
     hs2R7->Draw("AP");
     hs2R7->SetTitle("LHRS s2R8;run number;peak");

     TCanvas *c9=new TCanvas("c9");
     hs2R8->SetMarkerStyle(8);
     hs2R8->SetMarkerColor(4);
     hs2R8->Draw("AP");
     hs2R8->SetTitle("LHRS s2R9;run number;peak");

     TCanvas *c10=new TCanvas("c10");
     hs2R9->SetMarkerStyle(8);
     hs2R9->SetMarkerColor(4);
     hs2R9->Draw("AP");
     hs2R9->SetTitle("LHRS s2R10;run number;peak");

     TCanvas *c11=new TCanvas("c11");
     hs2R10->SetMarkerStyle(8);
     hs2R10->SetMarkerColor(4);
     hs2R10->Draw("AP");
     hs2R10->SetTitle("LHRS s2R11;run number;peak");

     TCanvas *c12=new TCanvas("c12");
     hs2R11->SetMarkerStyle(8);
     hs2R11->SetMarkerColor(4);
     hs2R11->Draw("AP");
     hs2R11->SetTitle("LHRS s2R12;run number;peak");

     TCanvas *c13=new TCanvas("c13");
     hs2R12->SetMarkerStyle(8);
     hs2R12->SetMarkerColor(4);
     hs2R12->Draw("AP");
     hs2R12->SetTitle("LHRS s2R13;run number;peak");

     TCanvas *c14=new TCanvas("c14");
     hs2R13->SetMarkerStyle(8);
     hs2R13->SetMarkerColor(4);
     hs2R13->Draw("AP");
     hs2R13->SetTitle("LHRS s2R14;run number;peak");

     TCanvas *c15=new TCanvas("c15");
     hs2R14->SetMarkerStyle(8);
     hs2R14->SetMarkerColor(4);
     hs2R14->Draw("AP");
     hs2R14->SetTitle("LHRS s2R15;run number;peak");

     TCanvas *c16=new TCanvas("c16");
     hs2R15->SetMarkerStyle(8);
     hs2R15->SetMarkerColor(4);
     hs2R15->Draw("AP");
     hs2R15->SetTitle("LHRS s2R16;run number;peak");

     c1->Print("Ls2R.pdf[");
     c1->Print("Ls2R.pdf");
     c2->Print("Ls2R.pdf");
     c3->Print("Ls2R.pdf");
     c4->Print("Ls2R.pdf");
     c5->Print("Ls2R.pdf");
     c6->Print("Ls2R.pdf");
     c7->Print("Ls2R.pdf");
     c8->Print("Ls2R.pdf");
     c9->Print("Ls2R.pdf");
     c10->Print("Ls2R.pdf");
     c11->Print("Ls2R.pdf");
     c12->Print("Ls2R.pdf");
     c13->Print("Ls2R.pdf");
     c14->Print("Ls2R.pdf");
     c15->Print("Ls2R.pdf");
     c16->Print("Ls2R.pdf");
     c16->Print("Ls2R.pdf]");






}


